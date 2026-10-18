//! The smooth gap-chaining seed heuristic, after A*PA (Groot Koerkamp and
//! Ivanov, Bioinformatics 2024), on POA graphs with affine and two-piece
//! costs. `docs/chain-heuristic-explained.md` explains it and
//! `docs/chain-heuristic-proofs.md` proves it admissible and consistent.
//!
//! In short: the query is cut into seeds of length k. Every node v has the
//! interval `[lo(v), hi(v)]` of path lengths to the graph end. A state x at
//! query offset i and node v gets a potential P(x), the charge for every
//! seed still ahead, and two numbers
//!
//! ```text
//! T(x) = ( e·(lo(v) + i) − P(x),  −e·(hi(v) + i) − P(x) )
//! ```
//!
//! with e the smallest gap extend. A seed match gets the T of the state
//! just after its first base. The heuristic is P(x) minus the savings of
//! the best chain of matches that the gap costs allow from x, which is a
//! dominance query on the contours.
//!
//! `seeds.rs` holds the query's seeds and the offset arithmetic around
//! them, `costs.rs` the seed charges, `trie.rs` the seeds as a trie,
//! `finder.rs` the sweep that finds every seed match, and `contours.rs` the
//! savings of chains of matches. `inspect.rs`, behind the `viz` feature,
//! splits `h` into its parts for tools that show the search.

pub(crate) mod contours;
pub(crate) mod costs;
pub(crate) mod finder;
#[cfg(any(test, feature = "viz"))]
mod inspect;
pub(crate) mod seeds;
pub(crate) mod trie;

use std::marker::PhantomData;
use std::sync::Arc;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

use contours::Contours;
pub use contours::Point;
use costs::SeedCosts;
use finder::{SeedFinder, SeedMatches};
use seeds::{Seeds, Walk};
use trie::SeedTrie;

#[cfg(feature = "viz")]
pub use contours::{MatchId, Staircase, SubPoint};
#[cfg(feature = "viz")]
pub use inspect::{ChainExplain, MatchInfo};

use super::{AstarHeuristic, HeuristicConfig, MinimumGapCost};
use crate::align::astar::aln_graph::{AlignmentGraphNode, resolve_bound};
use crate::align::cost_models::AlignmentCostModel;
use crate::align::offsets::OffsetType;
use crate::align::traits::AlignableGraph;
use crate::align::{AlignState, AlignmentMode};
use crate::graph::bubbles::index::BubbleIndex;
use crate::graph::traits::GraphNodeId;

/// Seed length, substitution budget and match pruning, as the CLI passes
/// them.
#[derive(Clone, Copy, Debug)]
pub struct SeedParams {
    pub k: usize,
    pub r: usize,
    /// Delete matches as the search passes them. Without it the heuristic
    /// is quadratic on long reads. Off only for diagnosis.
    pub prune_matches: bool,
}

/// Picks the smooth gap-chaining seed heuristic, taken as the larger of
/// itself and min-gap-cost.
#[derive(Clone, Copy, Debug)]
pub struct SmoothHeuristicConfig(pub SeedParams);

impl HeuristicConfig for SmoothHeuristicConfig {
    type Heuristic<N: GraphNodeId, O: OffsetType, C: AlignmentCostModel + Copy> = SmoothHeuristic<N, O, C>;

    fn instantiate<N, O, C>(&self, costs: C) -> SmoothHeuristic<N, O, C>
    where
        N: GraphNodeId,
        O: OffsetType,
        C: AlignmentCostModel + Copy,
    {
        SmoothHeuristic::new(costs, self.0)
    }
}

/// The parts of a state's potential.
#[derive(Clone, Copy, Debug, Default)]
pub struct Potential {
    /// Seeds after the state's offset, each charged in full.
    pub seeds_ahead: usize,
    /// Added for a half-aligned seed that no path can finish.
    pub partial: i64,
    /// Subtracted in an insertion state.
    pub ins_discount: i64,
}

impl Potential {
    pub fn value(&self, charge: i64) -> i64 {
        charge * self.seeds_ahead as i64 + self.partial - self.ins_discount
    }
}

/// What `prepare` builds for one alignment.
struct SeedIndex<N> {
    seeds: Seeds,
    costs: SeedCosts,
    bubbles: Arc<BubbleIndex<N>>,
    matches: SeedMatches<N>,
    contours: Contours,
}

impl<N: GraphNodeId> SeedIndex<N> {
    /// The potential at offset i in `state` if no path finishes the seed i
    /// is inside: an upper bound that needs no walk.
    #[inline]
    fn potential_bound(&self, i: usize, state: AlignState) -> Potential {
        let c = &self.costs;
        let seeds_ahead = self.seeds.ahead(i);
        match state {
            AlignState::Insertion | AlignState::Insertion2 => {
                Potential { seeds_ahead, partial: 0, ins_discount: c.ins_discount[state as usize][i % self.seeds.k] }
            }
            _ => {
                let partial = if self.seeds.inside(i) { c.partial[state as usize] } else { 0 };
                Potential { seeds_ahead, partial, ins_discount: 0 }
            }
        }
    }

    /// The potential of node v at offset i, from `p`, its `potential_bound`.
    fn potential<G>(&self, graph: &G, v: N, i: usize, mut p: Potential, walk: &mut Vec<Walk<N>>) -> Potential
    where
        G: AlignableGraph<Node = N>,
    {
        if p.partial > 0 && finder::rest_matches(graph, &self.seeds, v, i, walk) {
            p.partial = 0;
        }
        p
    }

    /// T of node v at offset i with potential p.
    fn point(&self, v: N, i: usize, p: i64) -> Point {
        let (lo, hi) = (self.bubbles.get_min_dist_to_end(v) as i64, self.bubbles.get_max_dist_to_end(v) as i64);
        let e = self.costs.e_min;
        Point { t1: e * (lo + i as i64) - p, t2: -e * (hi + i as i64) - p }
    }
}

/// Buffers kept between alignments.
struct Scratch<N> {
    trie: SeedTrie,
    finder: SeedFinder<N>,
    walk: Vec<Walk<N>>,
    /// The seed and match buffers while `index` is None.
    spare: (Seeds, SeedMatches<N>),
}

/// See the module docs.
pub struct SmoothHeuristic<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    params: SeedParams,
    costs: C,
    min_gap: MinimumGapCost<N, O, C>,
    /// None when the heuristic is off for this alignment. `h` is then
    /// min-gap-cost.
    index: Option<SeedIndex<N>>,
    scratch: Scratch<N>,
    _marker: PhantomData<O>,
}

impl<N, O, C> SmoothHeuristic<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel + Copy,
{
    pub fn new(costs: C, params: SeedParams) -> Self {
        assert!((1..=255).contains(&params.k), "seed length 1 to 255");
        Self {
            params,
            costs,
            min_gap: MinimumGapCost::unprepared(costs),
            index: None,
            scratch: Scratch { trie: SeedTrie::default(), finder: SeedFinder::default(), walk: Vec::new(), spare: Default::default() },
            _marker: PhantomData,
        }
    }

    /// Alignment-graph states for this cost model: 1, 3 or 5.
    // Read off the costs, since `AlignmentCostModel` does not
    // carry the state count. An affine model with a zero gap open counts as
    // linear here, which only lowers the charge.
    fn states(&self) -> usize {
        if self.costs.gap_extend2() > 0 {
            5
        } else if self.costs.gap_open() > 0 {
            3
        } else {
            1
        }
    }
}

impl<N, O, C> AstarHeuristic<N, O> for SmoothHeuristic<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel + Copy,
{
    fn prepare<G>(&mut self, graph: &G, query: &[u8], bubbles: &Arc<BubbleIndex<N>>, mode: AlignmentMode)
    where
        G: AlignableGraph<Node = N>,
    {
        // wasm32-unknown-unknown has no clock: `Instant::now` panics there.
        #[cfg(not(target_arch = "wasm32"))]
        let started = Instant::now();
        self.min_gap.prepare(graph, query, bubbles, mode);
        let SeedParams { k, r, .. } = self.params;
        let free_end = match mode {
            AlignmentMode::Global => 0,
            AlignmentMode::EndsFree { qry_free_end, .. } => resolve_bound(qry_free_end, query.len()),
        };
        let n = query.len().saturating_sub(free_end) / k;
        let costs = SeedCosts::new(&self.costs, self.states(), k, r, n);
        // Reuse the previous alignment's buffers.
        let (mut seeds, mut matches) = match self.index.take() {
            Some(ix) => (ix.seeds, ix.matches),
            None => std::mem::take(&mut self.scratch.spare),
        };
        if n == 0 || costs.charge <= 0 {
            tracing::debug!(target: "poasta::astar::seeds", n_seeds = n, "seed heuristic off for this query");
            self.scratch.spare = (seeds, matches);
            return;
        }

        seeds.reset(query, k, r, n);
        self.scratch.trie.rebuild(&seeds);
        self.scratch.finder.find(graph, &seeds, &self.scratch.trie, &mut matches);
        let mut ix = SeedIndex { seeds, costs, bubbles: bubbles.clone(), matches, contours: Contours::default() };

        // A match of seed j sits at the state just after its first base,
        // with the potential of the seeds after it.
        let c = ix.costs.charge;
        let mut points = Vec::with_capacity(ix.matches.len());
        for j in 0..n {
            let (i, p) = (ix.seeds.match_offset(j), c * (n - j - 1) as i64);
            points.extend(ix.matches.of_seed(j).iter().map(|m| ix.point(m.start, i, p)));
        }
        // Fresh contours per alignment. Reuse their allocations if
        // profiles of many short alignments show them.
        ix.contours = Contours::new(points, c, ix.costs.lattice);
        #[cfg(not(target_arch = "wasm32"))]
        let build_us = started.elapsed().as_micros() as u64;
        #[cfg(target_arch = "wasm32")]
        let build_us = 0u64;
        tracing::debug!(
            target: "poasta::astar::seeds",
            n_seeds = n,
            n_matches = ix.matches.len(),
            miss = ix.costs.miss,
            charge = c,
            build_us,
            "seed heuristic built",
        );
        self.index = Some(ix);
    }

    fn h<G>(&mut self, graph: &G, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize
    where
        G: AlignableGraph<Node = N>,
    {
        let h_min_gap = self.min_gap.bound(aln_node, aln_state);
        let Some(ix) = &mut self.index else {
            return h_min_gap;
        };
        let (v, i) = (aln_node.node(), aln_node.offset().as_usize());
        let charge = ix.costs.charge;
        // The chain part is at most the potential. Where even its upper
        // bound cannot beat min-gap-cost, skip the walk and the query.
        let bound = ix.potential_bound(i, aln_state);
        if bound.value(charge) <= h_min_gap as i64 {
            return h_min_gap;
        }
        let p = ix.potential(graph, v, i, bound, &mut self.scratch.walk).value(charge);
        if p <= h_min_gap as i64 {
            return h_min_gap;
        }
        let credit = ix.contours.credit(ix.point(v, i, p));
        (p - credit).max(h_min_gap as i64) as usize
    }

    fn on_expand<G>(&mut self, graph: &G, aln_node: &AlignmentGraphNode<N, O>, h_limit: usize) -> bool
    where
        G: AlignableGraph<Node = N>,
    {
        if !self.params.prune_matches {
            return false;
        }
        let Some(ix) = &self.index else {
            return false;
        };
        let (v, i) = (aln_node.node(), aln_node.offset().as_usize());
        let Some(j) = ix.seeds.seed_matched_at(i) else {
            return false;
        };
        let ids = ix.matches.starting_at(j, v);
        if ids.is_empty() {
            return false;
        }
        if h_limit != usize::MAX && self.h(graph, aln_node, AlignState::Match) > h_limit {
            return false;
        }
        let contours = &mut self.index.as_mut().expect("checked above").contours;
        let mut any = false;
        for m in ids {
            if contours.is_live(m) {
                contours.prune(m);
                tracing::trace!(target: "poasta::astar::match_pruned", m, node = v.index(), offset = i);
                any = true;
            }
        }
        any
    }

    fn is_dynamic(&self) -> bool {
        self.index.is_some() && self.params.prune_matches
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::astar::{AstarEngine, Offset};
    use crate::align::cost_models::affine::Affine;
    use crate::align::traits::AlignmentEngine;
    use crate::graph::alignment::AddAlignment;
    use crate::graph::poa::POAGraph;

    /// On real sequences the random test's graphs are too small. Build a
    /// dnaG graph and check every score against min-gap-cost, globally and
    /// with the middle of each sequence in semi-global mode.
    #[test]
    fn dna_g_scores_equal_min_gap_cost() {
        use crate::graph::io::seq::open_sequences;
        use std::ops::Bound;
        let seqs: Vec<Vec<u8>> = open_sequences(std::path::Path::new("tests/dnaG/align_seq.fna"))
            .unwrap()
            .take(12)
            .map(|r| r.unwrap().1)
            .collect();
        let semi = AlignmentMode::EndsFree {
            qry_free_begin: Bound::Excluded(0),
            qry_free_end: Bound::Excluded(0),
            graph_free_begin: Bound::Unbounded,
            graph_free_end: Bound::Unbounded,
        };
        let costs = Affine::new(4, 6, 2);
        let params = SeedParams { k: 8, r: 0, prune_matches: true };
        let mut aligner = crate::align::PoastaAligner::new(AstarEngine::<POAGraph<u32>, Affine>::new(costs), POAGraph::new());
        aligner.align_one_named("s0", &seqs[0]).unwrap();
        for (i, s) in seqs.iter().enumerate().skip(1) {
            let g = aligner.graph();
            for (mode, q) in [(AlignmentMode::Global, &s[..]), (semi, &s[s.len() / 4..3 * s.len() / 4])] {
                let want = AstarEngine::<POAGraph<u32>, Affine>::new(costs).with_mode(mode).align(g, q).unwrap().score;
                let got = AstarEngine::<POAGraph<u32>, Affine>::new(costs)
                    .with_mode(mode)
                    .with_heuristic(SmoothHeuristicConfig(params))
                    .align(g, q)
                    .unwrap()
                    .score;
                assert_eq!(got, want, "sequence {i}, {mode:?}");
            }
            aligner.align_one_named(&format!("s{i}"), s).unwrap();
        }
    }

    /// Greedy extension can pass states whose f has already grown past the
    /// popped priority, some of them at a suboptimal score. Pruning there
    /// deleted a match the better path needed: 12 instead of 11. The engine
    /// prunes only where f still equals the popped priority.
    #[test]
    fn greedy_extension_prunes_only_at_the_popped_priority() {
        let seqs: [&[u8]; 5] = [
            b"CAAACACAACACCCCCCACCCC",
            b"CAAACCAACACACCCCCCAACCCC",
            b"CAAACACAACACACCCCCCAACCCC",
            b"CAAACAAAACACACCCCCCAAACCCC",
            b"CAACACAACCACACCCCCCAACCCCC",
        ];
        let engine = AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(4, 6, 2)).with_heuristic(crate::align::astar::DijkstraHeuristic);
        let mut aligner = crate::align::PoastaAligner::new(engine, POAGraph::new());
        for (i, s) in seqs.iter().enumerate() {
            aligner.align_one_named(&format!("s{i}"), s).unwrap();
        }
        let g = aligner.into_graph();
        let costs = Affine::new(1, 2, 1);
        let q = b"AAAACCCCCACACC";
        let got = AstarEngine::<POAGraph<u32>, Affine>::new(costs)
            .with_heuristic(SmoothHeuristicConfig(SeedParams { k: 3, r: 0, prune_matches: true }))
            .without_pruning()
            .align(&g, q)
            .unwrap()
            .score;
        assert_eq!(got, 11);
    }

    #[test]
    fn short_query_turns_the_heuristic_off() {
        let mut g = POAGraph::<u32>::new();
        g.add_alignment("s", b"ACGTACGT", None, &[1; 8]).unwrap();
        let mut h = SmoothHeuristic::<_, Offset, _>::new(Affine::new(4, 6, 2), SeedParams { k: 8, r: 0, prune_matches: true });
        let bubbles = Arc::new(BubbleIndex::new(&g));
        h.prepare(&g, b"ACG", &bubbles, AlignmentMode::Global);
        assert!(h.index.is_none());
        let with = AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(4, 6, 2))
            .with_heuristic(SmoothHeuristicConfig(SeedParams { k: 8, r: 0, prune_matches: true }))
            .align(&g, b"ACG")
            .unwrap();
        let without = AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(4, 6, 2)).align(&g, b"ACG").unwrap();
        assert_eq!(with.score, without.score);
    }
}
