//! Scalar banded Gotoh affine-gap aligner for POA graphs.
//!
//! Uses Ukkonen-style band doubling combined with per-node multi-band tracking.
//! Each node may have multiple disjoint bands; bands are propagated from
//! predecessors (so reachable regions discovered on earlier nodes carry forward)
//! and merged when they are adjacent or overlap.
//!
//! # Band computation (per-node, topological order)
//! 1. Every real node that is a direct successor of the start sentinel receives an
//!    explicit seed band `[0, k.min(m)]` centred on the main diagonal.  In a linear
//!    graph this is exactly the first real node; in a bubble graph each chain-start
//!    node is seeded independently.
//! 2. Every real node derives additional bands from its real predecessors: a band
//!    `[lo_p, hi_p]` propagates as `[(lo_p+1).min(m), (hi_p+1).min(m)]`, shifting
//!    both bounds by one so the band width stays constant and the coverage tracks the
//!    alignment diagonal.  Capping at `m` prevents bands from becoming empty when the
//!    graph is longer than the query.
//! 3. All candidate intervals are sorted and merged when they overlap.
//!
//! # Memory layout
//! All bands across all nodes are stored in a single flat `Vec<Band>`, with a
//! per-node index into that slice (`node_band_range`).  The DP values (M, D, I)
//! for each band are packed into one contiguous `Vec<u32>` using prefix-sum
//! offsets.

use std::cell::RefCell;
use std::cmp::Ordering;
use std::convert::Infallible;
use std::marker::PhantomData;
use std::ops::Index;
use std::rc::Rc;

use crate::align::{
    cost_models::AlignmentCostModel,
    engine::{
        AlignOutput, AlignmentStats,
        dp::{INF, StateCol, backtrace_generic},
    },
    kernels::DPKernel,
    traits::{AlignableGraph, AlignmentEngine},
};
use crate::debug::logger::DebugLogger;
use crate::graph::{
    poa::{IndexType, POAGraph},
    traits::{GraphBase, GraphWithNodeOrdering},
};

fn log_align_stats(stats: &AlignmentStats) {
    tracing::info!(
        "Alignment stats - max k: {}, computed cells: {}, fraction of full matrix: {:.4}",
        stats.max_bandwidth,
        stats.cells_computed,
        stats.fraction_of_full_matrix,
    );
}

// ─── Engine struct ───────────────────────────────────────────────────────────

/// Reusable scratch buffers for a single `BandDoublingEngineScalar`. All fields
/// are reset on entry to `align_banded` and grow monotonically across calls so
/// repeated alignments against the same engine avoid heap churn.
#[derive(Default)]
pub(crate) struct Scratch {
    matrix: BandMatrix,
    init_data: Vec<u32>,
    node_to_rank: Vec<usize>,
    diag_buf: Vec<u32>,
}

pub struct BandDoublingEngineScalar<C, Ix>
where
    C: AlignmentCostModel,
    Ix: IndexType,
{
    costs: C,
    /// Starting bandwidth k (≥ 1).
    initial_k: usize,
    logger: Option<Rc<RefCell<DebugLogger>>>,
    scratch: RefCell<Scratch>,
    _phantom: PhantomData<Ix>,
}

impl<C, Ix> BandDoublingEngineScalar<C, Ix>
where
    C: AlignmentCostModel,
    Ix: IndexType,
{
    pub fn new(costs: C) -> Self {
        Self {
            costs,
            initial_k: 1,
            logger: None,
            scratch: RefCell::new(Scratch::default()),
            _phantom: PhantomData,
        }
    }

    /// Override the starting bandwidth (default 1).
    pub fn with_initial_k(mut self, k: usize) -> Self {
        self.initial_k = k.max(1);
        self
    }

    /// Attach a debug logger; call `DebugLogger::begin_sequence` before each `align` call.
    pub fn with_debug_logger(mut self, logger: Rc<RefCell<DebugLogger>>) -> Self {
        self.logger = Some(logger);
        self
    }
}

impl<C, Ix> AlignmentEngine<&[u8]> for BandDoublingEngineScalar<C, Ix>
where
    C: AlignmentCostModel,
    Ix: IndexType,
{
    type Graph = POAGraph<Ix>;
    type Success = AlignOutput<POAGraph<Ix>>;
    type Error = Infallible;

    fn align(
        &self,
        graph: &POAGraph<Ix>,
        query: &[u8],
    ) -> Result<AlignOutput<POAGraph<Ix>>, Infallible> {
        let m = query.len();
        let l_min = graph.l_min_real();
        let l_max = graph.l_max_real();
        let min_ge = self.costs.min_gap_extend() as usize;

        let mut k = self.initial_k.max(1);
        tracing::debug!(initial_k = k, "starting band-doubling loop");

        let forced = if m < l_min {
            l_min - m
        } else {
            m.saturating_sub(l_max)
        };

        // Total DP-cell budget (m+1) * (n_real + 1) * n_states, matching the
        // single-attempt formula in `align_banded`. Used to convert the
        // accumulated cell count into a fraction of the full matrix.
        let n_real = graph.node_count().saturating_sub(2);
        let full_matrix_cells = (m + 1)
            .saturating_mul(n_real + 1)
            .saturating_mul(C::Kernel::STATES);
        let mut total_cells_computed: usize = 0;

        let log_borrow = self.logger.as_ref().map(|rc| rc.borrow());
        let log: Option<&DebugLogger> = log_borrow.as_deref();

        let mut scratch = self.scratch.borrow_mut();

        let finalize_stats = |stats: &mut AlignmentStats, total_cells: usize| {
            stats.cells_computed = total_cells;
            stats.fraction_of_full_matrix = if full_matrix_cells == 0 {
                0.0
            } else {
                total_cells as f64 / full_matrix_cells as f64
            };
        };

        let mut iter_idx = 0usize;
        loop {
            let _k_span = tracing::info_span!("band_iter", iter = iter_idx, k).entered();
            if let Some(log) = log {
                log.begin_k_iter(k, C::Kernel::STATE_NAMES);
            }
            let attempt =
                align_banded::<C::Kernel, Ix>(&self.costs, graph, query, k, log, &mut scratch);
            // Account for the cells touched by this attempt regardless of
            // whether it reached the alignment end.
            total_cells_computed = total_cells_computed.saturating_add(scratch.matrix.data.len());

            if let Some(mut result) = attempt {
                // Ukkonen bound: with min_ge > 0, any alignment of cost s
                // contains at most s / min_ge indels, so the optimum deviates
                // from the main diagonal by at most that much; combined with
                // the length-forced indels we get required_k.
                let s = result.score as usize;
                let required_k = s
                    .checked_div(min_ge)
                    .map(|v| forced + v + 1)
                    .unwrap_or(usize::MAX);

                if k >= required_k {
                    tracing::debug!(score = result.score, k, required_k, "alignment optimal");
                    finalize_stats(&mut result.stats, total_cells_computed);
                    log_align_stats(&result.stats);
                    return Ok(result);
                }
                tracing::debug!(
                    score = result.score,
                    k,
                    required_k,
                    "band too narrow for provable optimum; doubling"
                );
            }

            let next_k = k.saturating_mul(2);
            if next_k >= m + l_max {
                // Once k covers the full matrix the DP is equivalent to canonical.
                let full_k = m + l_max;
                tracing::info!(
                    prev_k = k,
                    k = full_k,
                    "band width saturated; running full-width DP"
                );
                if let Some(log) = log {
                    log.begin_k_iter(full_k, C::Kernel::STATE_NAMES);
                }
                let mut result = align_banded::<C::Kernel, Ix>(
                    &self.costs,
                    graph,
                    query,
                    full_k,
                    log,
                    &mut scratch,
                )
                .expect("full-width band must succeed");
                total_cells_computed =
                    total_cells_computed.saturating_add(scratch.matrix.data.len());
                finalize_stats(&mut result.stats, total_cells_computed);
                log_align_stats(&result.stats);
                return Ok(result);
            }

            tracing::debug!(prev_k = k, k = next_k, "band width increased");
            k = next_k;
            iter_idx += 1;
        }
    }
}

// ─── Bands ────────────────────────────────────────────────────────────────────

/// A contiguous range of query positions associated with one graph node,
/// with indices to predecessor edges in the band DAG.
#[derive(Clone, Debug)]
pub(crate) struct Band {
    /// Topological rank of the owning node.
    pub(crate) node_rank: usize,

    /// First query position in this band (inclusive).
    pub(crate) qlo: usize,

    /// One past the last query position in this band (exclusive).
    /// The band covers the half-open interval `[qlo, qhi)`.
    pub(crate) qhi: usize,

    /// First query position of this band, allowed to be negative to ease band
    /// computation across nodes.
    qlo_signed: isize,

    /// Start index into `Bands::band_pred` for this band's predecessor edges.
    pred_ix_start: usize,

    /// Number of predecessor edges for this band.
    pred_num: usize,
}

impl Band {
    fn new(node_rank: usize, qlo_signed: isize, qhi: usize) -> Band {
        Band {
            node_rank,
            qlo_signed,
            qlo: qlo_signed.max(0) as usize,
            qhi,
            pred_ix_start: 0,
            pred_num: 0,
        }
    }

    fn new_with_pred(
        node_rank: usize,
        qlo_signed: isize,
        qhi: usize,
        pred_ix_start: usize,
    ) -> Self {
        Band {
            node_rank,
            qlo_signed,
            qlo: qlo_signed.max(0) as usize,
            qhi,
            pred_ix_start,
            pred_num: 1,
        }
    }

    #[inline]
    pub(crate) fn width(&self) -> usize {
        self.qhi - self.qlo
    }
}

impl PartialOrd for Band {
    fn partial_cmp(&self, other: &Band) -> Option<Ordering> {
        (self.qlo, self.qhi).partial_cmp(&(other.qlo, other.qhi))
    }
}

impl Eq for Band {}

impl PartialEq for Band {
    fn eq(&self, other: &Band) -> bool {
        (self.qlo, self.qhi) == (other.qlo, other.qhi)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BandEdge {
    /// Index of the predecessor band. Set at band-construction time; used by
    /// `BandMatrix::reset_for_bands` to fill `data_offset`, and retained as a
    /// debug breadcrumb.
    pred: usize,

    /// Offset into `BandMatrix::data` where the predecessor band's scores
    /// start. Filled in by `BandMatrix::reset_for_bands` once `n_states` is
    /// known; left as 0 during band construction.
    data_offset: usize,

    /// Width (`qhi - qlo`) of the predecessor band.
    width: usize,

    /// First query position of the predecessor band.
    qlo: usize,
}

impl BandEdge {
    fn new(pred: usize, width: usize, qlo: usize) -> Self {
        BandEdge {
            pred,
            data_offset: 0,
            width,
            qlo,
        }
    }
}

/// A set of bands through the alignment matrix.
///
/// This struct represents all bands of per-node query positions to compute
/// as a directed acylic graph, keeping track per band from which predecessor
/// bands it originated.
#[derive(Default)]
pub(crate) struct Bands {
    /// All computed bands
    bands: Vec<Band>,

    /// List of "edges", i.e., all predecessors across all bands in a single vector.
    ///
    /// Not to be confused with the edges in the POA graph itself.
    ///
    /// Use `Band.pred_ix_start` and `Band.pred_num` to find the range of edges
    /// belonging to a specific band.
    band_pred: Vec<BandEdge>,

    /// Per-node range of bands, represented as [start, end) index in `bands`.
    ///
    /// This vector enables mapping a POA node rank to its list of bands
    node_bands: Vec<(usize, usize)>,

    /// Working buffer for `reset_for_global_alignment`'s per-node band
    /// construction. Reused across nodes and across calls to avoid allocation.
    new_bands_tmp: Vec<Band>,

    /// Working buffer for collecting and sorting predecessor band indices
    /// during `reset_for_global_alignment`. Reused across nodes.
    pred_bands_tmp: Vec<usize>,
}

impl Bands {
    /// Recompute the band layout for a global alignment, reusing existing
    /// allocations. Equivalent to constructing a fresh `Bands`, but `Vec`
    /// capacity is preserved across calls.
    fn reset_for_global_alignment<G: AlignableGraph>(&mut self, graph: &G, query: &[u8], k: isize) {
        let bands = &mut self.bands;
        let band_pred = &mut self.band_pred;
        let node_bands = &mut self.node_bands;
        let new_bands = &mut self.new_bands_tmp;
        let pred_bands = &mut self.pred_bands_tmp;

        bands.clear();
        band_pred.clear();
        node_bands.clear();

        let query_max = query.len();

        // Initial band on start node (half-open [qlo, qhi)).
        let qhi = (k.unsigned_abs() + 1).max(2).min(query_max + 1);
        bands.push(Band::new(0, -k, qhi));
        node_bands.push((0usize, 1usize));

        // Process all other nodes in topological order
        for rank in 1..graph.node_count() {
            let node = graph.rank_to_node(rank);
            new_bands.clear();

            // Compute bands for the current node based on bands on predecessor nodes.
            //
            // Each predecessor's bands are already sorted by query position, but we
            // need a single global sort to detect overlaps across predecessors. For
            // the typical small fan-in (1-3 predecessors with ≤1 band each), a
            // collect-and-sort beats a heap-based k-way merge.
            //
            // Fan-in == 1 (the common case for chain-like POA regions) is handled
            // without sorting: the single predecessor's band range is already
            // qlo-ordered by construction.
            pred_bands.clear();
            let mut n_preds = 0usize;
            for pred in graph.predecessors(node) {
                n_preds += 1;
                let pred_rank = graph.node_rank(pred);
                let (pred_band_start, pred_band_end) = node_bands[pred_rank];
                pred_bands.extend(pred_band_start..pred_band_end);
            }
            if n_preds > 1 {
                pred_bands.sort_unstable_by(|&a, &b| {
                    (bands[a].qlo, bands[a].qhi).cmp(&(bands[b].qlo, bands[b].qhi))
                });
            }

            let mut curr_new_band: Option<Band> = None;
            for &pred_ix in pred_bands.iter() {
                let pred_band = &bands[pred_ix];
                // Move predecessor band one to the right (half-open [qlo, qhi)).
                let target_qlo = (pred_band.qlo_signed + 1).min(query_max as isize);
                let target_qhi = (pred_band.qhi + 1).min(query_max + 1);

                // Skip zero-width bands (can happen when a predecessor's band
                // already reaches the end of the query sequence).
                let effective_qlo = target_qlo.max(0) as usize;
                if target_qhi <= effective_qlo {
                    continue;
                }

                if let Some(ref mut target) = curr_new_band {
                    if pred_band.qlo < target.qhi {
                        // Overlap, merge bands
                        let prev_qhi = target.qhi;
                        target.qhi = target.qhi.max(target_qhi);
                        tracing::debug!(
                            node_rank = rank,
                            target_qlo = target.qlo,
                            prev_qhi,
                            merged_qhi = target.qhi,
                            "merged overlapping bands",
                        );

                        // Add edge
                        let edge = BandEdge::new(pred_ix, pred_band.width(), pred_band.qlo);
                        band_pred.push(edge);
                        target.pred_num += 1;
                    } else {
                        // No overlap, create new band
                        let new_band =
                            Band::new_with_pred(rank, target_qlo, target_qhi, band_pred.len());
                        let edge = BandEdge::new(pred_ix, pred_band.width(), pred_band.qlo);
                        band_pred.push(edge);

                        // Store the band currently in `curr_new_band` in `new_bands`, and replace
                        // the old value with the newly created `new_band`
                        let to_store = std::mem::replace(target, new_band);
                        new_bands.push(to_store);
                    }
                } else {
                    curr_new_band = Some(Band::new_with_pred(
                        rank,
                        target_qlo,
                        target_qhi,
                        band_pred.len(),
                    ));
                    let edge = BandEdge::new(pred_ix, pred_band.width(), pred_band.qlo);
                    band_pred.push(edge);
                }
            }

            if let Some(target) = curr_new_band {
                new_bands.push(target);
            }

            bands.append(&mut *new_bands);

            let prev_ix_start = node_bands.last().unwrap().1;
            node_bands.push((prev_ix_start, bands.len()))
        }

        tracing::debug!(
            k,
            n_bands = bands.len(),
            n_edges = band_pred.len(),
            n_nodes = node_bands.len(),
            "computed bands for all nodes",
        );
    }

    fn len(&self) -> usize {
        self.bands.len()
    }

    /// Return the slice of `BandEdge` describing the predecessor bands of
    /// the band at index `band_ix`.
    #[inline]
    fn pred_for_band(&self, band_ix: usize) -> &[BandEdge] {
        let band = &self.bands[band_ix];
        &self.band_pred[band.pred_ix_start..band.pred_ix_start + band.pred_num]
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = &Band> {
        self.bands.iter()
    }

    /// Returns `(absolute_band_ix, &Band)` pairs for a given node rank.
    #[inline]
    fn bands_for_node_indexed(&self, node_rank: usize) -> impl Iterator<Item = (usize, &Band)> {
        let (lo, hi) = self.node_bands[node_rank];
        (lo..hi).map(move |bi| (bi, &self.bands[bi]))
    }
}

impl Index<usize> for Bands {
    type Output = Band;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        self.bands.index(index)
    }
}

// ─── BandMatrix (runtime n_states) ────────────────────────────────────────

/// Runtime-flexible banded DP matrix (n_states set at construction time).
/// Used by `align_banded` to support kernels with any number of states.
#[derive(Default)]
pub(crate) struct BandMatrix {
    pub(crate) bands: Bands,
    pub(crate) band_data_offset: Vec<usize>,
    pub(crate) data: Vec<u32>,
    pub(crate) n_states: usize,
}

impl BandMatrix {
    /// Recompute the prefix-sum offset table and resize the score buffer to
    /// match the bands currently stored in `self.bands`. Reuses existing Vec
    /// allocations across calls to avoid heap churn.
    pub(crate) fn reset_for_bands(&mut self, n_states: usize) {
        self.band_data_offset.clear();
        let mut offset = 0usize;
        for band in self.bands.iter() {
            self.band_data_offset.push(offset);
            offset += n_states * band.width();
        }
        self.data.clear();
        self.data.resize(offset, INF);
        self.n_states = n_states;

        // Now that the offset table exists, populate `data_offset` on each
        // predecessor edge so `BandPredIterator` can build slices into `data`
        // without consulting `band_data_offset` and `bands` at iteration time.
        for edge in &mut self.bands.band_pred {
            edge.data_offset = self.band_data_offset[edge.pred];
        }
    }

    pub(crate) fn compute_iter(&mut self) -> BandComputeIterator<'_> {
        BandComputeIterator {
            matrix: self,
            current: 0,
        }
    }

    /// Return the column (all states at one query position) for `node_rank` at
    /// query position `q`, or `None` if no band of this node covers `q`.
    ///
    /// Used on the backtrace hot path. The common case is `hi - lo == 1` (one
    /// band per node), handled inline without the per-state min loop or the
    /// `[INF; MAX_STATES]` zero-init scaffolding. Multi-band nodes (only at
    /// band merges) fall through to the cold loop.
    #[inline]
    pub(crate) fn node_col_at(&self, node_rank: usize, q: usize) -> Option<StateCol> {
        let n = self.n_states;
        let (lo, hi) = self.bands.node_bands[node_rank];

        if hi - lo == 1 {
            let bi = lo;
            let b = &self.bands[bi];
            if q < b.qlo || q >= b.qhi {
                return None;
            }
            let off = self.band_data_offset[bi];
            let w = b.width();
            let qi = q - b.qlo;
            let mut col = [INF; crate::align::engine::dp::MAX_STATES];
            for s in 0..n {
                col[s] = self.data[off + s * w + qi];
            }
            return Some(col);
        }

        self.node_col_at_multi(node_rank, q)
    }

    #[cold]
    #[inline(never)]
    fn node_col_at_multi(&self, node_rank: usize, q: usize) -> Option<StateCol> {
        let n = self.n_states;
        let (lo, hi) = self.bands.node_bands[node_rank];
        let mut col = [INF; crate::align::engine::dp::MAX_STATES];
        let mut any = false;
        for bi in lo..hi {
            let b = &self.bands[bi];
            if q < b.qlo || q >= b.qhi {
                continue;
            }
            let off = self.band_data_offset[bi];
            let w = b.width();
            let qi = q - b.qlo;
            for s in 0..n {
                let v = self.data[off + s * w + qi];
                if v < col[s] {
                    col[s] = v;
                }
            }
            any = true;
        }
        if any { Some(col) } else { None }
    }
}

/// Mutable view of one band's data in a BandMatrix.
///
/// Yielded by `BandComputeIterator::next_band` for the band currently being
/// computed. Holds a mutable slice into `BandMatrix::data` for this band's
/// scores plus the metadata needed to lazily construct read-only slices for
/// the band's predecessors.
pub(crate) struct BandDataSlice<'a> {
    /// The band being computed.
    pub band: &'a Band,

    /// Mutable scores for this band. Length = `n_states * band.width()`.
    pub data: &'a mut [u32],

    /// Predecessor edges feeding this band.
    pub pred_bands: &'a [BandEdge],

    /// Pointer to the start of `BandMatrix::data`, used by `pred_iter` to
    /// construct read-only slices into earlier bands' regions of the buffer.
    data_ptr: *const u32,

    pub n_states: usize,
}

impl<'a> BandDataSlice<'a> {
    pub fn pred_iter(&self) -> BandPredIterator<'a> {
        BandPredIterator {
            edges: self.pred_bands,
            data_ptr: self.data_ptr,
            n_states: self.n_states,
            current: 0,
        }
    }
}

pub(crate) struct BandPredIterator<'a> {
    edges: &'a [BandEdge],
    data_ptr: *const u32,
    n_states: usize,
    current: usize,
}

impl<'a> Iterator for BandPredIterator<'a> {
    type Item = (usize, &'a [u32]);

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.edges.get(self.current)?;
        self.current += 1;

        // SAFETY: Bands are built in topological order, so every predecessor
        // edge points at a band strictly earlier than the current band.
        // `BandMatrix::band_data_offset` is a strict prefix-sum over
        // `n_states * width` chunks, so regions for different bands do not
        // overlap. The predecessor slice `[edge.data_offset, edge.data_offset +
        // n_states * edge.width)` is therefore entirely disjoint from the
        // current band's mutable slice held by the parent `BandDataSlice`.
        unsafe {
            let slice: &[u32] = std::slice::from_raw_parts(
                self.data_ptr.add(edge.data_offset),
                self.n_states * edge.width,
            );
            Some((edge.qlo, slice))
        }
    }
}

/// The compute iterator enables computing alignment scores for each band
/// using the DP recursion.
///
/// This is a *lending* iterator (not `std::iter::Iterator`): each call to
/// `next_band` returns a `BandDataSlice` whose lifetime is tied to the
/// borrow of `self`. The previous slice must therefore be dropped before
/// the next call. Use a `while let Some(slice) = it.next_band() { ... }`
/// loop in the caller.
pub(crate) struct BandComputeIterator<'a> {
    matrix: &'a mut BandMatrix,
    current: usize,
}

impl<'a> BandComputeIterator<'a> {
    pub(crate) fn next_band(&mut self) -> Option<BandDataSlice<'_>> {
        let b = self.current;
        if b >= self.matrix.bands.len() {
            return None;
        }
        self.current += 1;

        let n = self.matrix.n_states;
        let off = self.matrix.band_data_offset[b];
        let pred_edges = self.matrix.bands.pred_for_band(b);
        let band = &self.matrix.bands[b];
        let w = band.width();

        let ptr = self.matrix.data.as_mut_ptr();

        // SAFETY: `ptr.add(off)` points at the start of band `b`'s exclusive
        // chunk in `BandMatrix::data`, sized `n * w`. The chunks for
        // different bands are disjoint by the prefix-sum offset layout, so
        // returning a `&mut [u32]` over this range does not alias any other
        // band's mutable region. Predecessor slices materialized via
        // `BandPredIterator` are read-only and target strictly earlier
        // bands' chunks, which are also disjoint from this one.
        let data = unsafe { std::slice::from_raw_parts_mut(ptr.add(off), n * w) };

        Some(BandDataSlice {
            band,
            data,
            pred_bands: pred_edges,
            data_ptr: ptr as *const u32,
            n_states: n,
        })
    }
}

// ─── Forward pass ────────────────────────────────────────────────────────────

/// Attempt global alignment within bandwidth `k`. Returns `None` if the band
/// was too narrow to reach the end of the alignment (triggers band doubling).
fn align_banded<K: DPKernel, Ix: IndexType>(
    costs: &K::Costs,
    graph: &POAGraph<Ix>,
    query: &[u8],
    k: usize,
    logger: Option<&DebugLogger>,
    scratch: &mut Scratch,
) -> Option<AlignOutput<POAGraph<Ix>>>
where
    K::Costs: AlignmentCostModel,
{
    let m = query.len();
    let mp1 = m + 1;

    tracing::debug!(
        k,
        query_len = m,
        n_states = K::STATES,
        "align_banded attempt"
    );

    // Split scratch into disjoint field borrows so the DP loop can hold
    // `&mut matrix` and `&mut diag_buf` simultaneously.
    let Scratch {
        matrix,
        init_data,
        node_to_rank,
        diag_buf,
    } = scratch;

    // 1. Build band structure and (re-)size the flat DP matrix in place.
    matrix
        .bands
        .reset_for_global_alignment(graph, query, k as isize);

    // Write band debug output (outside the hot DP loop).
    if let Some(log) = logger {
        for band in &matrix.bands.bands {
            log.write_band(k, band.node_rank, band.qlo, band.qhi);
        }
    }

    matrix.reset_for_bands(K::STATES);

    // 2. Initialize band 0 (start sentinel) via kernel.
    {
        let off = matrix.band_data_offset[0];
        let w = matrix.bands[0].width();
        let qlo = matrix.bands[0].qlo;
        K::init_start(&mut matrix.data[off..off + K::STATES * w], w, qlo, costs);
    }

    // 3. Build init_data for the start sentinel used by the generic backtrace.
    //    Layout: `init_data[s * mp1 + q]`, length = K::STATES * mp1.
    let init_data_len = K::STATES * mp1;
    init_data.clear();
    init_data.resize(init_data_len, INF);
    K::init_start(init_data, mp1, 0, costs);

    // 4. Build node_to_rank: node.index() → topological rank.
    let max_idx = graph
        .all_nodes_iter()
        .map(|nd| nd.index())
        .max()
        .unwrap_or(0);
    node_to_rank.clear();
    node_to_rank.resize(max_idx + 1, usize::MAX);
    for topo_idx in 0..graph.node_count() {
        let v = graph.rank_to_node(topo_idx);
        node_to_rank[v.index()] = topo_idx;
    }

    // 5. Forward DP over all bands in topological order.
    //
    // Band 0 is the start sentinel and is skipped via `.skip(1)` (it was
    // initialized above via `init_start`). The end sentinel's bands are the
    // last entries in `bands` (it has the highest topological rank), so we
    // stop the loop as soon as we encounter them — the final alignment score
    // is read separately from the end node's predecessors below.
    let end = graph.end_node();
    let end_rank = graph.node_rank(end);

    // Scratch buffer for the diagonal predecessor values used by M[q].
    // Pre-allocated at maximum band width to avoid per-band heap churn.
    let max_band_width = matrix
        .bands
        .bands
        .iter()
        .map(|b| b.width())
        .max()
        .unwrap_or(1);
    if diag_buf.len() < max_band_width {
        diag_buf.resize(max_band_width, INF);
    }

    let mut compute_iter = matrix.compute_iter();
    let _ = compute_iter.next_band(); // skip start sentinel (band 0)
    while let Some(band_slices) = compute_iter.next_band() {
        let node_rank = band_slices.band.node_rank;
        if node_rank == end_rank {
            break;
        }
        let node = graph.rank_to_node(node_rank);

        let symbol = graph.get_node_symbol(node);
        let qlo = band_slices.band.qlo;
        let qhi = band_slices.band.qhi;
        let w = qhi - qlo;

        // ── Fused per-predecessor pass: D-state + diagonal buffer ────────────
        let diag_buf = &mut diag_buf[..w];
        diag_buf.fill(INF);

        for (pred_qlo, pred_data) in band_slices.pred_iter() {
            let pred_w = pred_data.len() / K::STATES;
            K::accumulate_pred(
                band_slices.data,
                diag_buf,
                w,
                qlo,
                pred_data,
                pred_w,
                pred_qlo,
                costs,
            );
        }

        // ── Finalize: apply substitution and sequential M/I states ────────────
        K::finalize(band_slices.data, w, qlo, symbol, query, diag_buf, costs);

        // Write per-cell debug output (streamed via BufWriter).
        if let Some(log) = logger {
            let data: &[u32] = band_slices.data;
            let mut cell_buf = [0i32; 8];
            for qi in 0..w {
                for s in 0..K::STATES {
                    cell_buf[s] = data[s * w + qi] as i32;
                }
                log.write_cell(k, node_rank, qlo + qi, &cell_buf[..K::STATES]);
            }
        }

        tracing::trace!(
            node_rank,
            symbol,
            qlo,
            qhi,
            width = w,
            data = ?band_slices.data,
            "dp cell slice computed",
        );
    }

    // 6. Check alignment success: find the best score at (end's predecessor, m).
    let start = graph.start_node();
    let matrix_ref = &matrix;
    let (best_score, best_node, best_state) = graph
        .predecessors(end)
        .filter(|&pred| pred != start)
        .flat_map(|pred| {
            let pred_rank = graph.node_rank(pred);
            matrix_ref
                .bands
                .bands_for_node_indexed(pred_rank)
                .filter(move |(_, belem)| belem.qlo <= m && m < belem.qhi)
                .map(move |(bi, belem)| {
                    let off = matrix_ref.band_data_offset[bi];
                    let w = belem.width();
                    let qi = m - belem.qlo;
                    let data = &matrix_ref.data[off..off + K::STATES * w];
                    let col = crate::align::engine::dp::extract_col::<K>(data, w, qi);
                    let (score, state) = K::best_state(&col[..K::STATES]);
                    (score, pred, state)
                })
        })
        .fold((INF, start, 0u8), |(bs, bn, bst), (local, pred, state)| {
            if local < bs {
                (local, pred, state)
            } else {
                (bs, bn, bst)
            }
        });

    if best_score >= INF {
        tracing::debug!(k, "band too narrow to reach end; triggering doubling");
        return None; // band too narrow → trigger doubling
    }

    tracing::debug!(
        best_score,
        best_state,
        "alignment end reached; starting backtrace"
    );

    // 7. Backtrace using the generic kernel-aware backtrace from dp.rs.
    //    Column-view closure: reads n_states cells at column q from the bands
    //    that cover `q` for the given node. No per-step allocation.
    let get_node_col = |node: <POAGraph<Ix> as AlignableGraph>::Node,
                        q: usize|
     -> Option<crate::align::engine::dp::StateCol> {
        let rank = node_to_rank[node.index()];
        matrix_ref.node_col_at(rank, q)
    };

    let alignment = backtrace_generic::<K, POAGraph<Ix>, _>(
        costs,
        graph,
        &get_node_col,
        init_data,
        best_node,
        best_state,
        m,
    );

    let cells_computed = matrix.data.len();
    let n_real = graph.node_count().saturating_sub(2);
    let full_matrix_cells = mp1.saturating_mul(n_real + 1).saturating_mul(K::STATES);
    let fraction_of_full_matrix = if full_matrix_cells == 0 {
        0.0
    } else {
        cells_computed as f64 / full_matrix_cells as f64
    };

    Some(AlignOutput {
        score: best_score,
        alignment,
        stats: AlignmentStats {
            max_bandwidth: k,
            cells_computed,
            fraction_of_full_matrix,
        },
    })
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{
        cost_models::affine::Affine, engine::dp::CanonicalDP, traits::AlignmentEngine,
    };
    use crate::graph::{alignment::AddAlignment, poa::POAGraph};

    fn scores() -> Affine {
        // match=0, mismatch=1, gap_open=2, gap_extend=1
        Affine::new(1, 2, 1)
    }

    fn linear_graph(seq: &[u8]) -> POAGraph<u32> {
        let mut g = POAGraph::new();
        g.add_alignment("s0", seq, None, &vec![1; seq.len()])
            .unwrap();
        g
    }

    fn bubble_graph(seq1: &[u8], seq2: &[u8]) -> POAGraph<u32> {
        let mut g = POAGraph::new();
        g.add_alignment("s0", seq1, None, &vec![1; seq1.len()])
            .unwrap();
        g.add_alignment("s1", seq2, None, &vec![1; seq2.len()])
            .unwrap();
        g
    }

    fn run_banded(graph: &POAGraph<u32>, query: &[u8]) -> AlignOutput<POAGraph<u32>> {
        BandDoublingEngineScalar::<Affine, u32>::new(scores())
            .align(graph, query)
            .unwrap()
    }

    fn run_canonical(graph: &POAGraph<u32>, query: &[u8]) -> AlignOutput<POAGraph<u32>> {
        CanonicalDP::<Affine, POAGraph<u32>>::new(scores())
            .align(graph, query)
            .unwrap()
    }

    #[test]
    fn test_banded_exact_match() {
        let g = linear_graph(b"ACGT");
        let r = run_banded(&g, b"ACGT");
        assert_eq!(r.score, 0, "exact match should score 0");
    }

    #[test]
    fn test_banded_stats_populated() {
        let g = linear_graph(b"ACGTACGT");
        let r = run_banded(&g, b"ACGTACGT");
        let s = r.stats;
        let m = 8usize;
        let n_real = 8usize;
        assert!(s.cells_computed > 0, "cells_computed must be positive");
        assert!(
            s.max_bandwidth > m.abs_diff(n_real),
            "max_bandwidth={} must cover |m-n|+1",
            s.max_bandwidth
        );
        assert!(
            s.fraction_of_full_matrix > 0.0 && s.fraction_of_full_matrix <= 1.0,
            "fraction out of range: {}",
            s.fraction_of_full_matrix
        );
    }

    #[test]
    fn test_canonical_stats_are_full_matrix() {
        let g = linear_graph(b"ACGTACGT");
        let r = run_canonical(&g, b"ACGTACGT");
        let s = r.stats;
        assert!((s.fraction_of_full_matrix - 1.0).abs() < 1e-9);
        assert!(s.cells_computed > 0);
    }

    #[test]
    fn test_banded_one_mismatch() {
        let g = linear_graph(b"ACGT");
        let r = run_banded(&g, b"ACXT");
        assert_eq!(r.score, 1, "one mismatch should cost 1");
    }

    #[test]
    fn test_banded_one_deletion() {
        let g = linear_graph(b"ACGT");
        let r = run_banded(&g, b"ACT");
        assert_eq!(r.score, 3, "one deletion: open(2)+extend(1)=3");
    }

    #[test]
    fn test_banded_one_insertion() {
        let g = linear_graph(b"ACT");
        let r = run_banded(&g, b"ACGT");
        assert_eq!(r.score, 3, "one insertion: open(2)+extend(1)=3");
    }

    #[test]
    fn test_banded_long_deletion() {
        let g = linear_graph(b"ACCCGT");
        let r = run_banded(&g, b"AGT");
        assert_eq!(r.score, 5, "gap of 3: open(2)+3*extend(1)=5");
    }

    #[test]
    fn test_banded_bubble_longer_path() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        let r = run_banded(&g, b"ACGAT");
        assert_eq!(r.score, 0, "query matches longer path exactly");
    }

    #[test]
    fn test_banded_bubble_shorter_path() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        let r = run_banded(&g, b"ACGT");
        assert_eq!(r.score, 0, "query matches shorter path exactly");
    }

    #[test]
    fn test_banded_matches_canonical_linear() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGT", b"ACGT"),
            (b"ACGT", b"ACXT"),
            (b"ACGT", b"ACT"),
            (b"TTTT", b"TTTTTTTT"),
            (b"AAAA", b"AAAA"),
            (b"ACGTACGT", b"ACGTXCGT"),
        ];
        for (seq, query) in cases {
            let g = linear_graph(seq);
            let banded = run_banded(&g, query);
            let canonical = run_canonical(&g, query);
            assert_eq!(
                banded.score,
                canonical.score,
                "banded score mismatch for seq={} query={}",
                std::str::from_utf8(seq).unwrap(),
                std::str::from_utf8(query).unwrap()
            );
        }
    }

    #[test]
    fn test_banded_matches_canonical_bubble() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        let queries: &[&[u8]] = &[b"ACGAT", b"ACGT", b"ACT", b"ACGXT"];
        for query in queries {
            let banded = run_banded(&g, query);
            let canonical = run_canonical(&g, query);
            assert_eq!(
                banded.score,
                canonical.score,
                "bubble banded score mismatch for query={}",
                std::str::from_utf8(query).unwrap()
            );
        }
    }

    #[test]
    fn align_banded_with_kernel_matches_canonical_bubble() {
        let cases: &[(&[u8], &[u8], &[u8])] = &[
            (b"ACGAT", b"ACGT", b"ACGAT"),
            (b"ACGAT", b"ACGT", b"ACGT"),
            (b"ACGAT", b"ACGT", b"ACXGT"),
        ];
        for (s1, s2, q) in cases {
            let g = bubble_graph(s1, s2);
            let banded = run_banded(&g, q);
            let canonical = run_canonical(&g, q);
            assert_eq!(
                banded.score,
                canonical.score,
                "kernel mismatch: s1={} s2={} q={}",
                std::str::from_utf8(s1).unwrap(),
                std::str::from_utf8(s2).unwrap(),
                std::str::from_utf8(q).unwrap()
            );
        }
    }

    #[test]
    fn ukkonen_detects_internal_indels() {
        // Scores where mismatches are far more costly than gap pairs, so the
        // optimum has off-diagonal indels even though |m − n_real| == 0.
        // Without the Ukkonen termination check, align_banded at k=1 returns
        // a non-INF but suboptimal score (all mismatches on the diagonal)
        // and the old loop would accept it. With the check, required_k grows
        // with the score gap and forces doubling until the optimum is
        // provably inside the band.
        let sc = Affine::new(100, 1, 1);
        let mut g: POAGraph<u32> = POAGraph::new();
        g.add_alignment("s0", b"AAAACCCC", None, &[1; 8]).unwrap();
        let query = b"CCCCAAAA";

        let engine: BandDoublingEngineScalar<Affine, u32> =
            BandDoublingEngineScalar::new(sc).with_initial_k(1);
        let banded = engine.align(&g, query).unwrap();

        let canonical: CanonicalDP<Affine, POAGraph<u32>> = CanonicalDP::new(sc);
        let expected = canonical.align(&g, query).unwrap();

        assert_eq!(
            banded.score, expected.score,
            "banded must match canonical even when initial_k is too small"
        );
    }

    #[test]
    fn ukkonen_doubling_count_small_for_easy_input() {
        // Identical query and graph: score 0, required_k = 0/min_ge + 1 = 1,
        // so the doubling loop must return on the first iteration at k=1.
        let g = linear_graph(b"ACGTACGT");
        let query = b"ACGTACGT";
        let engine: BandDoublingEngineScalar<Affine, u32> =
            BandDoublingEngineScalar::new(scores()).with_initial_k(1);
        let result = engine.align(&g, query).unwrap();
        assert_eq!(result.score, 0);
        assert!(
            result.stats.max_bandwidth <= 2,
            "easy input should not double beyond k=1 (got {})",
            result.stats.max_bandwidth
        );
    }

    #[test]
    fn dyn_band_matrix_allocates_n_states_per_position() {
        let g = linear_graph(b"ACGT");
        let query = b"ACGT";
        let mut mat3 = BandMatrix::default();
        mat3.bands.reset_for_global_alignment(&g, query, 4isize);
        let total_3: usize = mat3.bands.iter().map(|b| 3 * b.width()).sum();
        mat3.reset_for_bands(3);
        assert_eq!(mat3.data.len(), total_3);

        let mut mat1 = BandMatrix::default();
        mat1.bands.reset_for_global_alignment(&g, query, 4isize);
        let total_1: usize = mat1.bands.iter().map(|b| b.width()).sum();
        mat1.reset_for_bands(1);
        assert_eq!(mat1.data.len(), total_1);
    }
}
