//! Canonical O(N×m) Gotoh affine-gap DP for POA graphs.
//!
//! Used as a ground-truth oracle to verify banded/SIMD aligner correctness.

use std::convert::Infallible;
use std::marker::PhantomData;

use crate::align::{
    cost_models::AlignmentCostModel,
    engine::{AlignResult, AlignedPair, AlignmentStats},
    kernels::{BacktraceOp, DPKernel},
    traits::{AlignableGraph, AlignmentEngine},
};
use crate::graph::traits::GraphNodeId;

/// Maximum DP state count across all supported kernels.
///
/// Used to size stack-allocated column snapshots (`[u32; MAX_STATES]`) passed
/// to the column-based backtrace API. Larger than all current `DPKernel::STATES`
/// values (1 linear, 3 affine, 5 two-piece); unused slots stay at [`INF`] and
/// are never read by kernel impls (which only index `col[0..STATES]`).
pub(crate) const MAX_STATES: usize = 5;

/// Snapshot of all `STATES` DP values at one `(node, q)` position.
pub(crate) type StateCol = [u32; MAX_STATES];

/// Extract the column at query position `q` from a node's flat band data.
///
/// `data` has layout `data[s * w + qi]`; the returned array holds
/// `col[s] = data[s * w + qi]` for `s in 0..STATES`, and leaves the remaining
/// slots at [`INF`]. Used by callers that already hold a dense band buffer
/// (forward-pass best-score extraction), not on the backtrace hot path.
#[inline]
pub(crate) fn extract_col<K: DPKernel>(data: &[u32], w: usize, qi: usize) -> StateCol {
    let mut col = [INF; MAX_STATES];
    for s in 0..K::STATES {
        col[s] = data[s * w + qi];
    }
    col
}

/// Extract the column at query position `q` from `init_data` (length
/// `STATES * mp1`).
#[inline]
pub(crate) fn extract_init_col<K: DPKernel>(init_data: &[u32], mp1: usize, q: usize) -> StateCol {
    let mut col = [INF; MAX_STATES];
    for s in 0..K::STATES {
        col[s] = init_data[s * mp1 + q];
    }
    col
}

/// Sentinel value for "unreachable / infinite cost".
///
/// Chosen so that `INF + any_realistic_cost` stays below `u32::MAX`, which lets
/// `sat` use a plain wrapping add + min instead of the more expensive
/// `saturating_add` (which requires overflow detection on every u32 lane and
/// maps to multiple instructions instead of `vpaddd` + `vpminud`).
///
/// Upper bound on a real score: `sequence_len × max_cost_per_step`.  For
/// sequences up to ~10 million bp and costs up to ~1 000, that is ≈10^10 —
/// far below `1<<28 ≈ 2.7×10^8`.  If you align sequences of length > 100 M or
/// use very large cost values, raise this constant accordingly.
pub(crate) const INF: u32 = 1 << 28;

/// Saturating-style cost addition, capped at [`INF`].
///
/// Uses `wrapping_add` rather than `saturating_add` so that the compiler can
/// lower the vectorised form to a single `vpaddd` + `vpminud` pair instead of
/// the multi-instruction overflow-detection sequence that `u32::saturating_add`
/// requires.  Correctness relies on the invariant that real scores never exceed
/// `u32::MAX - max_single_step_cost`; see [`INF`] for the bound.
#[inline(always)]
pub(crate) fn sat(a: u32, b: u32) -> u32 {
    a.wrapping_add(b).min(INF)
}

#[inline(always)]
pub(crate) fn min3(a: u32, b: u32, c: u32) -> u32 {
    a.min(b).min(c)
}

// ─── Public engine struct ────────────────────────────────────────────────────

pub struct CanonicalDP<C, G> {
    costs: C,
    _graph: PhantomData<G>,
}

impl<C: AlignmentCostModel, G: AlignableGraph> CanonicalDP<C, G> {
    pub fn new(costs: C) -> Self {
        Self {
            costs,
            _graph: PhantomData,
        }
    }
}

impl<C, G> AlignmentEngine<&[u8]> for CanonicalDP<C, G>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    type Graph = G;
    type Success = AlignResult<G>;
    type Error = Infallible;

    fn align(&self, graph: &G, query: &[u8]) -> Result<AlignResult<G>, Infallible> {
        Ok(align_dp::<C::Kernel, G>(&self.costs, graph, query))
    }
}

// ─── Generic forward pass ────────────────────────────────────────────────────

fn align_dp<K: DPKernel, G: AlignableGraph>(
    costs: &K::Costs,
    graph: &G,
    query: &[u8],
) -> AlignResult<G> {
    let m = query.len();
    let n = graph.node_count();
    let n_states = K::STATES;
    let mp1 = m + 1;
    let start = graph.start_node();
    let end = graph.end_node();

    let _seq_span = tracing::info_span!(
        "align_sequence",
        engine = "canonical_dp",
        query_len = m,
        graph_nodes = n,
        n_states,
    )
    .entered();

    // Build node_to_rank: maps node.index() → topological rank (0-indexed among
    // non-sentinel nodes). n_real = number of non-sentinel nodes.
    let max_idx = graph
        .all_nodes_iter()
        .map(|nd| nd.index())
        .max()
        .unwrap_or(0);
    let mut node_to_rank = vec![usize::MAX; max_idx + 1];
    let mut n_real = 0usize;
    for topo_idx in 0..n {
        let v = graph.rank_to_node(topo_idx);
        if v == start || v == end {
            continue;
        }
        node_to_rank[v.index()] = n_real;
        n_real += 1;
    }

    // Allocate flat DP storage: n_real * n_states * mp1.
    // Layout: dp_data[rank * n_states * mp1 + s * mp1 + q]
    let dp_cells = n_real * n_states * mp1;
    tracing::debug!(
        n_real,
        n_states,
        mp1,
        cells = dp_cells,
        "allocating full DP matrix",
    );
    let mut dp_data = vec![INF; dp_cells];

    // Build init_data for the start sentinel: n_states * mp1 elements.
    // We use a band of width mp1 starting at qlo=0.
    let mut init_data = vec![INF; n_states * mp1];
    K::init_start(&mut init_data, mp1, 0, costs);

    // Helper: compute base offset for a real node by rank.
    let slot_base = |rank: usize| rank * n_states * mp1;

    // Forward pass.
    let mut real_rank = 0usize;
    for topo_idx in 0..n {
        let v = graph.rank_to_node(topo_idx);
        if v == start || v == end {
            continue;
        }

        let slot = real_rank;
        real_rank += 1;
        let symbol = graph.get_node_symbol(v);

        // Collect predecessor data into temporary buffers to satisfy the borrow
        // checker (current slot and predecessors both live in dp_data; they do not
        // overlap in practice but Rust cannot verify that statically).
        let pred_bufs: Vec<Vec<u32>> = graph
            .predecessors(v)
            .filter(|&pred| pred != end)
            .map(|pred| {
                if pred == start {
                    init_data.to_vec()
                } else {
                    let r = node_to_rank[pred.index()];
                    let base = slot_base(r);
                    dp_data[base..base + n_states * mp1].to_vec()
                }
            })
            .collect();

        // Accumulate gap states from all predecessors.
        for src in &pred_bufs {
            K::accumulate_gap(
                &mut dp_data[slot_base(slot)..slot_base(slot) + n_states * mp1],
                mp1,
                0,
                src,
                mp1,
                0,
                costs,
            );
        }

        // Accumulate diagonal predecessor scores into diag_buf.
        let mut diag_buf = vec![INF; mp1];
        for src in &pred_bufs {
            K::accumulate_diag(&mut diag_buf, mp1, 0, src, mp1, 0);
        }

        // Finalize: substitution + sequential M/I recurrence.
        K::finalize(
            &mut dp_data[slot_base(slot)..slot_base(slot) + n_states * mp1],
            mp1,
            0,
            symbol,
            query,
            &mut diag_buf,
            costs,
        );

        tracing::trace!(
            node_rank = slot,
            symbol,
            n_preds = pred_bufs.len(),
            data = ?&dp_data[slot_base(slot)..slot_base(slot) + n_states * mp1],
            "dp cell slice computed",
        );
    }

    // Find best score at predecessors of end node at query position m.
    let mut best_score = INF;
    let mut best_node = start;
    let mut best_state_id = 0u8;

    for pred in graph.predecessors(end) {
        if pred == start {
            continue;
        }
        let rank = node_to_rank[pred.index()];
        let base = slot_base(rank);
        let data = &dp_data[base..base + n_states * mp1];
        let col = extract_col::<K>(data, mp1, m);
        let (score, state) = K::best_state(&col[..n_states]);
        if score < best_score {
            best_score = score;
            best_node = pred;
            best_state_id = state;
        }
    }

    tracing::debug!(
        best_score,
        best_state = best_state_id,
        "alignment end reached; starting backtrace",
    );

    // Backtrace using the generic backtrace — column-view closure reads only
    // the N_STATES cells at column `q`, no per-step allocation.
    let dp_data_ref = &dp_data;
    let node_to_rank_ref = &node_to_rank;
    let get_node_col = |node: G::Node, q: usize| -> Option<StateCol> {
        if q >= mp1 {
            return None;
        }
        let rank = node_to_rank_ref[node.index()];
        let base = rank * n_states * mp1;
        let mut col = [INF; MAX_STATES];
        for s in 0..n_states {
            col[s] = dp_data_ref[base + s * mp1 + q];
        }
        Some(col)
    };

    let alignment = backtrace_generic::<K, G>(
        costs,
        graph,
        &get_node_col,
        &init_data,
        best_node,
        best_state_id,
        m,
    );

    let cells_computed = dp_cells;
    let full_matrix_cells = mp1.saturating_mul(n_real).saturating_mul(n_states);
    let fraction_of_full_matrix = if full_matrix_cells == 0 {
        0.0
    } else {
        cells_computed as f64 / full_matrix_cells as f64
    };

    AlignResult {
        score: best_score,
        alignment,
        stats: AlignmentStats {
            max_bandwidth: m.saturating_add(n_real),
            cells_computed,
            fraction_of_full_matrix,
        },
    }
}

// ─── Generic backtrace ───────────────────────────────────────────────────────

/// Generic backtrace using `K::backtrace_op` for per-cell state machine decisions.
///
/// `get_node_col(node, q)` returns the column of length `STATES` at query
/// position `q` in the given node's band data, or `None` if `q` lies outside
/// every band of `node` (equivalent to "this cell was never computed").
///
/// `init_data` has `K::STATES * (m+1)` elements and represents the start
/// sentinel column (built via `K::init_start`).
pub(crate) fn backtrace_generic<K: DPKernel, G: AlignableGraph>(
    costs: &K::Costs,
    graph: &G,
    get_node_col: &dyn Fn(G::Node, usize) -> Option<StateCol>,
    init_data: &[u32],
    best_node: G::Node,
    best_state: u8,
    m: usize,
) -> Vec<AlignedPair<G::Node>> {
    let n_states = K::STATES;
    // Backtrace walks at most one step per graph node per query position;
    // reserving `m + 32` avoids most reallocations without over-allocating.
    let mut pairs = Vec::with_capacity(m + 32);
    let mut v = best_node;
    let mut state = best_state;
    let mut q = m;

    loop {
        // Column at q for the backtrace decision.
        let col = get_node_col(v, q).unwrap_or([INF; MAX_STATES]);
        let col_prev_storage;
        let col_prev = if q > 0 {
            col_prev_storage = get_node_col(v, q - 1).unwrap_or([INF; MAX_STATES]);
            Some(&col_prev_storage[..n_states])
        } else {
            None
        };

        let op = K::backtrace_op(&col[..n_states], col_prev, state, costs);

        match op {
            BacktraceOp::Insert { next_state } => {
                debug_assert!(q > 0, "q=0 in Insert state");
                pairs.push(AlignedPair::new(None, Some(q - 1)));
                q -= 1;
                state = next_state;
            }
            BacktraceOp::Diagonal { next_state } => {
                if q == 0 {
                    // At q=0 in M/Diagonal state, redirect to deletion path
                    state = next_state;
                    continue;
                }
                pairs.push(AlignedPair::new(Some(v), Some(q - 1)));
                q -= 1;
                match find_best_pred_generic::<K, G>(graph, get_node_col, init_data, v, q, costs) {
                    None => break,
                    Some((pv, ps)) => {
                        v = pv;
                        state = ps;
                    }
                }
            }
            BacktraceOp::Delete { .. } => {
                pairs.push(AlignedPair::new(Some(v), None));
                match find_del_pred_generic::<K, G>(graph, get_node_col, init_data, v, q, costs) {
                    None => break,
                    Some((pv, ps)) => {
                        v = pv;
                        state = ps;
                    }
                }
            }
            BacktraceOp::Done => break,
        }
    }

    pairs.reverse();
    pairs
}

/// Find the predecessor of `v` with the best (minimum) overall score at query
/// position `q`, across all states. Returns `None` if the best predecessor is
/// the start sentinel (alignment is complete), otherwise `Some((pred_node, best_state_id))`.
fn find_best_pred_generic<K: DPKernel, G: AlignableGraph>(
    graph: &G,
    get_node_col: &dyn Fn(G::Node, usize) -> Option<StateCol>,
    init_data: &[u32],
    v: G::Node,
    q: usize,
    _costs: &K::Costs,
) -> Option<(G::Node, u8)> {
    let start = graph.start_node();
    let end = graph.end_node();
    let mp1 = init_data.len() / K::STATES;
    let n_states = K::STATES;
    let mut best_score = INF;
    let mut best: Option<(G::Node, u8)> = None;
    let mut best_is_start = false;

    for pred in graph.predecessors(v) {
        if pred == end {
            continue;
        }
        if pred == start {
            if q < mp1 {
                let col = extract_init_col::<K>(init_data, mp1, q);
                let (score, _state) = K::best_state(&col[..n_states]);
                if score < best_score {
                    best_score = score;
                    best_is_start = true;
                    best = None;
                }
            }
        } else if let Some(col) = get_node_col(pred, q) {
            let (score, state) = K::best_state(&col[..n_states]);
            if score < best_score {
                best_score = score;
                best_is_start = false;
                best = Some((pred, state));
            }
        }
    }

    if best_is_start {
        None
    } else {
        best
    }
}

/// Find the predecessor of `v` that minimises the deletion transition cost at
/// query position `q`, using the exact same recurrence as the forward pass:
///
/// `D[v][q] = min(sat(M_pred[q].min(I_pred[q]), go+ge), sat(D_pred[q], ge))`
///
/// Returns `None` if the best predecessor is the start sentinel (alignment is
/// complete), otherwise `Some((pred_node, pred_state_id))` where `pred_state_id`
/// is the state in the predecessor that produced the minimum cost.
fn find_del_pred_generic<K: DPKernel, G: AlignableGraph>(
    graph: &G,
    get_node_col: &dyn Fn(G::Node, usize) -> Option<StateCol>,
    init_data: &[u32],
    v: G::Node,
    q: usize,
    costs: &K::Costs,
) -> Option<(G::Node, u8)> {
    let start = graph.start_node();
    let end = graph.end_node();
    let mp1 = init_data.len() / K::STATES;
    let n_states = K::STATES;
    let mut best_cost = INF;
    let mut best: Option<(G::Node, u8)> = None;

    for pred in graph.predecessors(v) {
        if pred == end {
            continue;
        }
        let col = if pred == start {
            if q >= mp1 {
                continue;
            }
            extract_init_col::<K>(init_data, mp1, q)
        } else {
            match get_node_col(pred, q) {
                Some(c) => c,
                None => continue,
            }
        };
        let (cost, pred_state) = K::del_transition_cost(&col[..n_states], costs);
        if cost < best_cost {
            best_cost = cost;
            if pred == start {
                best = None;
            } else {
                best = Some((pred, pred_state));
            }
        }
    }
    best
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{cost_models::affine::Affine, traits::AlignmentEngine};
    use crate::graph::{alignment::AddAlignment, poa::POAGraph};

    // Tests for LinearKernel/TwoPieceAffineKernel via CanonicalDP will be added in Tasks 6/7

    fn costs() -> Affine {
        // match=0, mismatch=1, gap_open=2, gap_extend=1
        Affine::new(0, 1, 2, 1)
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

    fn run(graph: &POAGraph<u32>, query: &[u8]) -> AlignResult<POAGraph<u32>> {
        let engine: CanonicalDP<Affine, POAGraph<u32>> = CanonicalDP::new(costs());
        engine.align(graph, query).unwrap()
    }

    #[test]
    fn test_exact_match() {
        let g = linear_graph(b"ACGT");
        let r = run(&g, b"ACGT");
        assert_eq!(r.score, 0, "exact match should score 0");
    }

    #[test]
    fn test_one_mismatch() {
        let g = linear_graph(b"ACGT");
        let r = run(&g, b"ACXT");
        assert_eq!(r.score, 1, "one mismatch should score 1");
    }

    #[test]
    fn test_one_deletion() {
        // ACGT vs ACT: G deleted → gap_open+gap_extend=3
        let g = linear_graph(b"ACGT");
        let r = run(&g, b"ACT");
        assert_eq!(r.score, 3, "one deletion: open(2)+extend(1)=3");
    }

    #[test]
    fn test_one_insertion() {
        // ACT vs ACGT: extra G in query → gap_open+gap_extend=3
        let g = linear_graph(b"ACT");
        let r = run(&g, b"ACGT");
        assert_eq!(r.score, 3, "one insertion: open(2)+extend(1)=3");
    }

    #[test]
    fn test_long_deletion() {
        // ACCCGT (6) vs AGT (3): A + del(CCC) + G + T → open(2)+3*extend(1)=5
        let g = linear_graph(b"ACCCGT");
        let r = run(&g, b"AGT");
        assert_eq!(r.score, 5, "gap of 3: open(2)+3*extend(1)=5");
    }

    #[test]
    fn test_bubble_match_longer_path() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        let r = run(&g, b"ACGAT");
        assert_eq!(r.score, 0, "query matches longer path exactly");
    }

    #[test]
    fn test_bubble_match_shorter_path() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        let r = run(&g, b"ACGT");
        assert_eq!(r.score, 0, "query matches shorter path exactly");
    }

    #[test]
    fn test_score_consistent_with_alignment() {
        let g = linear_graph(b"ACGT");
        let query = b"AXGT";
        let r = run(&g, query);
        let recomputed = score_from_alignment(&g, query, &r.alignment, &costs());
        assert_eq!(
            r.score, recomputed,
            "score from alignment pairs must match result score"
        );
    }

    fn score_from_alignment(
        graph: &POAGraph<u32>,
        query: &[u8],
        aln: &[AlignedPair<<POAGraph<u32> as AlignableGraph>::Node>],
        costs: &Affine,
    ) -> u32 {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        let mut score = 0u32;
        let mut in_gap = false;

        for pair in aln {
            match (pair.node(), pair.query_pos()) {
                (Some(node), Some(qpos)) => {
                    let gsym = graph.get_node_symbol(node);
                    let qsym = query[qpos];
                    score = sat(
                        score,
                        if gsym == qsym {
                            costs.equal() as u32
                        } else {
                            costs.mismatch() as u32
                        },
                    );
                    in_gap = false;
                }
                (Some(_), None) | (None, Some(_)) => {
                    score = sat(score, if in_gap { ge } else { go + ge });
                    in_gap = true;
                }
                (None, None) => {}
            }
        }
        score
    }
}
