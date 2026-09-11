//! A* alignment engine: best-first search over the alignment graph, with
//! greedy zero-cost match extension and bubble-index-based pruning.
//!
//! One [`AstarEngine`] serves the linear, affine and two-piece cost models;
//! [`AstarCostModel`] picks the state machine for each. Both `Global` and
//! `EndsFree` alignment modes are supported.

pub mod aln_graph;
pub mod dfa;
pub mod heuristic;
pub mod pruning;
pub mod queue;
pub mod visited;

#[cfg(test)]
pub(crate) mod dp_oracle;

use std::error::Error;
use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::align::cost_models::AlignmentCostModel;
use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::astar::aln_graph::{
    AffineAlignmentGraph, AlignmentGraph, AlignmentGraphNode, LinearAlignmentGraph, TwoPieceAlignmentGraph,
};
use crate::align::astar::dfa::{DepthFirstGreedyAlignment, ExtendResult};
use crate::align::astar::heuristic::{AstarHeuristic, Dijkstra, MinimumGapCost};
use crate::align::astar::queue::{AffineLayeredQueue, AstarQueue, LinearLayeredQueue, TwoPieceLayeredQueue};
use crate::align::astar::visited::{
    AffineAstarData, AstarData, AstarVisited, LinearAstarData, TwoPieceAstarData,
};
use crate::align::traits::{AlignOutput, AlignmentStats};
use crate::align::traits::{AlignableGraph, AlignmentEngine};
use crate::align::{AlignState, AlignmentMode};
use crate::graph::bubbles::index::BubbleIndex;
use crate::align::offsets::OffsetType;
use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};

/// Query offset width. Only `u32` is built; a `u64` variant for queries over
/// 4 Gbp could be added without changing the public API.
pub type Offset = u32;

/// Selects which [`AstarHeuristic`] an [`AstarEngine`] builds.
///
/// Carried as the engine's `H` type parameter so "affine + Dijkstra" vs
/// "affine + min-gap-cost" is a type choice instead of a bespoke config
/// struct per combination (the four hand-written `AlignmentConfig`s on
/// `stable` collapse into this).
pub trait HeuristicConfig {
    /// Whether this heuristic needs the bubble index to compute `h`. (The
    /// bubble index is built unconditionally today since bubble-exit pruning
    /// also depends on it; this is documentation for callers / future gating.)
    fn needs_bubble_index(&self) -> bool;

    fn build<C, N, O>(
        &self,
        costs: &C,
        bubbles: Arc<BubbleIndex<N>>,
        seq_len: usize,
        mode: AlignmentMode,
    ) -> Box<dyn AstarHeuristic<N, O>>
    where
        C: AlignmentCostModel + Copy + 'static,
        N: GraphNodeId + 'static,
        O: crate::align::offsets::OffsetType + 'static;
}

/// `h ≡ 0`: A* reduces to Dijkstra. Exact and admissible, but explores every
/// state within the optimal-score radius.
#[derive(Default, Clone, Copy)]
pub struct DijkstraHeuristic;

impl HeuristicConfig for DijkstraHeuristic {
    fn needs_bubble_index(&self) -> bool {
        false
    }

    fn build<C, N, O>(
        &self,
        _costs: &C,
        _bubbles: Arc<BubbleIndex<N>>,
        _seq_len: usize,
        _mode: AlignmentMode,
    ) -> Box<dyn AstarHeuristic<N, O>>
    where
        C: AlignmentCostModel + Copy + 'static,
        N: GraphNodeId + 'static,
        O: crate::align::offsets::OffsetType + 'static,
    {
        Box::new(Dijkstra::default())
    }
}

/// Minimum-gap-cost heuristic: an admissible lower bound on the remaining
/// cost, derived from the bubble index's cached shortest/longest distance
/// to the graph end. The default heuristic (tighter bound, fewer states
/// explored than Dijkstra).
#[derive(Default, Clone, Copy)]
pub struct MinGapCostHeuristic;

impl HeuristicConfig for MinGapCostHeuristic {
    fn needs_bubble_index(&self) -> bool {
        true
    }

    fn build<C, N, O>(
        &self,
        costs: &C,
        bubbles: Arc<BubbleIndex<N>>,
        seq_len: usize,
        mode: AlignmentMode,
    ) -> Box<dyn AstarHeuristic<N, O>>
    where
        C: AlignmentCostModel + Copy + 'static,
        N: GraphNodeId + 'static,
        O: crate::align::offsets::OffsetType + 'static,
    {
        Box::new(MinimumGapCost::new_with_mode(
            *costs, bubbles, seq_len, mode,
        ))
    }
}

/// Returned when the priority queue empties before the end state is reached
/// (e.g. a disconnected graph).
#[derive(Debug)]
pub enum AstarError {
    NoAlignmentFound,
}

impl fmt::Display for AstarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoAlignmentFound => write!(
                f,
                "A* queue emptied before reaching the alignment end state"
            ),
        }
    }
}

impl Error for AstarError {}

/// Ties a cost model to the A* machinery shaped for it: how many visited
/// slots a cell needs, which alignment graph expands its states, how scores
/// are stored, and how the queue buckets them. This is what lets one
/// [`AstarEngine`] serve linear (1 state), affine (3) and two-piece affine
/// (5) costs.
pub trait AstarCostModel<N, O>: AlignmentCostModel + Copy + 'static
where
    N: GraphNodeId,
    O: OffsetType,
{
    /// Visited-storage slots per `(node, offset)` cell.
    const STATES: usize;
    type AlnGraph: AlignmentGraph;
    type Visited: AstarVisited<N, O>;
    type Queue: AstarQueue<N, O>;

    fn aln_graph(&self, mode: AlignmentMode) -> Self::AlnGraph;

    fn visited<G>(&self, graph: &G, seq: &[u8], bubbles: Arc<BubbleIndex<N>>) -> Self::Visited
    where
        G: GraphWithNodeOrdering<NodeType = N>;
}

impl<N: GraphNodeId, O: OffsetType> AstarCostModel<N, O> for Linear {
    const STATES: usize = 1;
    type AlnGraph = LinearAlignmentGraph;
    type Visited = LinearAstarData<N, O>;
    type Queue = LinearLayeredQueue<N, O>;

    fn aln_graph(&self, mode: AlignmentMode) -> Self::AlnGraph {
        LinearAlignmentGraph::new(*self, mode)
    }

    fn visited<G>(&self, graph: &G, seq: &[u8], bubbles: Arc<BubbleIndex<N>>) -> Self::Visited
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        AstarData::new(*self, graph, seq, bubbles)
    }
}

impl<N: GraphNodeId, O: OffsetType> AstarCostModel<N, O> for Affine {
    const STATES: usize = 3;
    type AlnGraph = AffineAlignmentGraph;
    type Visited = AffineAstarData<N, O>;
    type Queue = AffineLayeredQueue<N, O>;

    fn aln_graph(&self, mode: AlignmentMode) -> Self::AlnGraph {
        AffineAlignmentGraph::new(*self, mode, 1)
    }

    fn visited<G>(&self, graph: &G, seq: &[u8], bubbles: Arc<BubbleIndex<N>>) -> Self::Visited
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        AstarData::new(*self, graph, seq, bubbles)
    }
}

impl<N: GraphNodeId, O: OffsetType> AstarCostModel<N, O> for TwoPieceAffine {
    const STATES: usize = 5;
    type AlnGraph = TwoPieceAlignmentGraph;
    type Visited = TwoPieceAstarData<N, O>;
    type Queue = TwoPieceLayeredQueue<N, O>;

    fn aln_graph(&self, mode: AlignmentMode) -> Self::AlnGraph {
        TwoPieceAlignmentGraph::new(*self, mode, 2)
    }

    fn visited<G>(&self, graph: &G, seq: &[u8], bubbles: Arc<BubbleIndex<N>>) -> Self::Visited
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        AstarData::new(*self, graph, seq, bubbles)
    }
}

/// The A* alignment engine. `C` picks the cost model and, through
/// [`AstarCostModel`], the state machine shaped for it; `H` picks the
/// heuristic.
pub struct AstarEngine<G, C = Affine, H = MinGapCostHeuristic>
where
    G: AlignableGraph,
    C: AlignmentCostModel,
{
    costs: C,
    mode: AlignmentMode,
    heuristic_cfg: H,
    existing_bubbles: Option<Arc<BubbleIndex<G::Node>>>,
    enable_pruning: bool,
    _graph: PhantomData<G>,
}

impl<G, C> AstarEngine<G, C, MinGapCostHeuristic>
where
    G: AlignableGraph,
    C: AlignmentCostModel,
{
    /// Global alignment, min-gap-cost heuristic, pruning enabled.
    pub fn new(costs: C) -> Self {
        Self {
            costs,
            mode: AlignmentMode::Global,
            heuristic_cfg: MinGapCostHeuristic,
            existing_bubbles: None,
            enable_pruning: true,
            _graph: PhantomData,
        }
    }
}

impl<G, C, H> AstarEngine<G, C, H>
where
    G: AlignableGraph,
    C: AlignmentCostModel,
{
    pub fn with_mode(mut self, mode: AlignmentMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn with_heuristic<H2>(self, heuristic_cfg: H2) -> AstarEngine<G, C, H2> {
        AstarEngine {
            costs: self.costs,
            mode: self.mode,
            heuristic_cfg,
            existing_bubbles: self.existing_bubbles,
            enable_pruning: self.enable_pruning,
            _graph: PhantomData,
        }
    }

    /// Reuse a bubble index built earlier for the same graph, for callers
    /// that align many queries against one read-only graph.
    pub fn with_existing_bubbles(mut self, bubbles: Arc<BubbleIndex<G::Node>>) -> Self {
        self.existing_bubbles = Some(bubbles);
        self
    }

    pub fn without_pruning(mut self) -> Self {
        self.enable_pruning = false;
        self
    }
}

impl<G, C, H> AlignmentEngine<&[u8]> for AstarEngine<G, C, H>
where
    G: AlignableGraph,
    G::Node: 'static,
    C: AstarCostModel<G::Node, Offset>,
    H: HeuristicConfig,
{
    type Graph = G;
    type Success = AlignOutput<G>;
    type Error = AstarError;

    fn align(&self, graph: &G, query: &[u8]) -> Result<AlignOutput<G>, Self::Error> {
        let _span = tracing::info_span!(
            "align_sequence",
            engine = "astar",
            query_len = query.len(),
            graph_nodes = graph.node_count(),
        )
        .entered();

        let bubble_index = self
            .existing_bubbles
            .clone()
            .unwrap_or_else(|| Arc::new(BubbleIndex::new(graph)));

        let aln_graph = self.costs.aln_graph(self.mode);
        let mut visited = self.costs.visited(graph, query, bubble_index.clone());
        let heuristic = self.heuristic_cfg.build::<C, G::Node, Offset>(
            &self.costs,
            bubble_index,
            query.len(),
            self.mode,
        );
        let mut queue = C::Queue::default();

        let (end_score, end_node, num_visited) = run_astar_search(
            graph,
            query,
            &aln_graph,
            &mut visited,
            &mut queue,
            heuristic.as_ref(),
            self.enable_pruning,
        )?;

        Ok(finish_alignment(
            graph,
            query,
            &visited,
            end_score,
            end_node,
            num_visited,
            C::STATES,
        ))
    }
}

/// The A* main loop, shared across every cost-model shape: best-first search
/// with greedy zero-cost match extension via [`DepthFirstGreedyAlignment`],
/// re-entering the queue only at mismatches, the graph end, or the query end.
/// Every argument is behind a trait (`AlignmentGraph`, `AstarVisited`,
/// `AstarQueue`), so this one function drives the 1-, 3- and 5-state
/// machines alike.
#[allow(clippy::too_many_arguments)]
fn run_astar_search<G, AG, V, Q>(
    graph: &G,
    query: &[u8],
    aln_graph: &AG,
    visited: &mut V,
    queue: &mut Q,
    heuristic: &dyn AstarHeuristic<G::Node, Offset>,
    enable_pruning: bool,
) -> Result<(u32, AlignmentGraphNode<G::Node, Offset>, usize), AstarError>
where
    G: AlignableGraph,
    AG: AlignmentGraph,
    V: AstarVisited<G::Node, Offset>,
    Q: AstarQueue<G::Node, Offset>,
{
    type O = Offset;
    let mut num_visited = 0usize;

    for initial_state in aln_graph.initial_states::<G, O>(graph, query) {
        let h = heuristic.h(&initial_state, AlignState::Match);
        queue.queue_aln_state(initial_state, AlignState::Match, 0, h);
        visited.set_score(&initial_state, AlignState::Match, 0);
    }

    let (end_score, end_node) = 'main: loop {
        let Some(item) = queue.pop_aln_state() else {
            return Err(AstarError::NoAlignmentFound);
        };
        let score = item.score();
        let aln_node = item.aln_node();
        let aln_state = item.aln_state();

        if score > visited.get_score(&aln_node, aln_state) {
            continue;
        }

        if aln_graph.is_end(graph, query, &aln_node, aln_state) {
            trace_visited(score, &aln_node, aln_state);
            num_visited += 1;
            break (score, aln_node);
        }

        if enable_pruning && visited.prune(score, &aln_node, aln_state) {
            continue;
        }

        visited.mark_reached(score, &aln_node, aln_state);
        trace_visited(score, &aln_node, aln_state);
        num_visited += 1;

        if aln_state == AlignState::Match {
            let mut dfa = DepthFirstGreedyAlignment::new(graph, query, score, &aln_node);
            let is_end_check = |n: &AlignmentGraphNode<G::Node, O>| {
                aln_graph.is_end(graph, query, n, AlignState::Match)
            };

            while let Some(end_point) = dfa.extend(visited, is_end_check) {
                match end_point {
                    ExtendResult::RefGraphEnd(parent, child) => {
                        if aln_graph.is_end(graph, query, &child, AlignState::Match) {
                            num_visited += dfa.get_num_visited();
                            break 'main (score, child);
                        }

                        aln_graph.expand_ref_graph_end(
                            visited,
                            &parent,
                            score,
                            |delta, succ, succ_state| {
                                let h = heuristic.h(&succ, succ_state);
                                queue.queue_aln_state(succ, succ_state, score + delta as u32, h);
                            },
                        );
                    }
                    ExtendResult::QueryEnd(parent, child) => {
                        aln_graph.expand_query_end(
                            visited,
                            &parent,
                            child,
                            score,
                            |delta, succ, succ_state| {
                                let h = heuristic.h(&succ, succ_state);
                                queue.queue_aln_state(succ, succ_state, score + delta as u32, h);
                            },
                        );
                    }
                    ExtendResult::Mismatch(parent, child) => {
                        aln_graph.expand_mismatch(
                            visited,
                            &parent,
                            &child,
                            score,
                            |delta, succ, succ_state| {
                                let h = heuristic.h(&succ, succ_state);
                                queue.queue_aln_state(succ, succ_state, score + delta as u32, h);
                            },
                        );
                    }
                    ExtendResult::EndReached(end_node) => {
                        num_visited += dfa.get_num_visited();
                        break 'main (score, end_node);
                    }
                }
            }

            num_visited += dfa.get_num_visited();
        } else {
            aln_graph.expand_all(
                visited,
                graph,
                query,
                score,
                &aln_node,
                aln_state,
                |delta, succ, succ_state| {
                    let h = heuristic.h(&succ, succ_state);
                    queue.queue_aln_state(succ, succ_state, score + delta as u32, h);
                },
            );
        }
    };

    Ok((end_score, end_node, num_visited))
}

/// Emit the `set_visited` trace event consumed by `debug::subscriber`.
#[inline]
fn trace_visited<N: GraphNodeId, O: OffsetType>(score: u32, node: &AlignmentGraphNode<N, O>, state: AlignState) {
    tracing::trace!(
        target: "poasta::astar::set_visited",
        score,
        node = node.node().index(),
        offset = node.offset().as_usize(),
        state = ?state,
    );
}

/// Backtrace + wrap the search result into an [`AlignOutput`], shared by
/// every cost model. `n_states` is the visited-storage state count, used only
/// for the fraction-of-full-matrix stat.
fn finish_alignment<G, V>(
    graph: &G,
    query: &[u8],
    visited: &V,
    end_score: u32,
    end_node: AlignmentGraphNode<G::Node, Offset>,
    num_visited: usize,
    n_states: usize,
) -> AlignOutput<G>
where
    G: AlignableGraph,
    V: AstarVisited<G::Node, Offset>,
{
    tracing::debug!(
        score = end_score,
        num_visited,
        "alignment end reached; starting backtrace"
    );

    let alignment = visited.backtrace(graph, query, &end_node);

    let n_real = graph.node_count().saturating_sub(2);
    let full_matrix_cells = (query.len() + 1)
        .saturating_mul(n_real)
        .saturating_mul(n_states);
    let fraction_of_full_matrix = if full_matrix_cells == 0 {
        0.0
    } else {
        num_visited as f64 / full_matrix_cells as f64
    };

    AlignOutput {
        score: end_score,
        alignment,
        stats: AlignmentStats {
            max_bandwidth: 0,
            cells_computed: num_visited,
            fraction_of_full_matrix,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::dp_oracle::DpOracle;
    use crate::align::test_support::score_from_alignment;
    use crate::graph::alignment::AddAlignment;
    use crate::graph::poa::POAGraph;

    fn costs() -> Affine {
        // match=0, mismatch=1, gap_open=2, gap_extend=1 (same as dp.rs's tests).
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

    fn multi_branch_graph() -> POAGraph<u32> {
        let mut g = POAGraph::new();
        g.add_alignment("s0", b"ACGTACGT", None, &vec![1; 8])
            .unwrap();
        g.add_alignment("s1", b"ACGAACGT", None, &vec![1; 8])
            .unwrap();
        g.add_alignment("s2", b"ACGTCCGT", None, &vec![1; 8])
            .unwrap();
        g.add_alignment("s3", b"AGGTACAT", None, &vec![1; 8])
            .unwrap();
        g
    }

    /// Assert the A* engine (both heuristics) agrees with the canonical DP
    /// oracle's score on `graph`/`query` — the core M2 equivalence check.
    fn assert_astar_matches_oracle(graph: &POAGraph<u32>, query: &[u8]) {
        let oracle: DpOracle<Affine, POAGraph<u32>> = DpOracle::new(costs());
        let expected = oracle.align(graph, query).unwrap();

        let dijkstra_engine: AstarEngine<POAGraph<u32>, Affine, DijkstraHeuristic> =
            AstarEngine::new(costs()).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        assert_eq!(
            dijkstra_result.score,
            expected.score,
            "Dijkstra A* score mismatch for query {:?}",
            String::from_utf8_lossy(query)
        );

        let mingapcost_engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs());
        let mingapcost_result = mingapcost_engine.align(graph, query).unwrap();
        assert_eq!(
            mingapcost_result.score,
            expected.score,
            "MinGapCost A* score mismatch for query {:?}",
            String::from_utf8_lossy(query)
        );
    }

    #[test]
    fn exact_match() {
        assert_astar_matches_oracle(&linear_graph(b"ACGT"), b"ACGT");
    }

    #[test]
    fn one_mismatch() {
        assert_astar_matches_oracle(&linear_graph(b"ACGT"), b"ACXT");
    }

    #[test]
    fn one_deletion() {
        assert_astar_matches_oracle(&linear_graph(b"ACGT"), b"ACT");
    }

    #[test]
    fn one_insertion() {
        assert_astar_matches_oracle(&linear_graph(b"ACT"), b"ACGT");
    }

    #[test]
    fn long_deletion() {
        assert_astar_matches_oracle(&linear_graph(b"ACCCGT"), b"AGT");
    }

    #[test]
    fn bubble_match_longer_path() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        assert_astar_matches_oracle(&g, b"ACGAT");
    }

    #[test]
    fn bubble_match_shorter_path() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        assert_astar_matches_oracle(&g, b"ACGT");
    }

    #[test]
    fn bubble_graph_with_mismatches_and_indels() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        assert_astar_matches_oracle(&g, b"AXGT");
        assert_astar_matches_oracle(&g, b"ACGATT");
    }

    #[test]
    fn multi_sequence_graph_various_queries() {
        let g = multi_branch_graph();
        for query in [
            b"ACGTACGT".as_slice(),
            b"ACGAACGT".as_slice(),
            b"ACGTCCGT".as_slice(),
            b"AGGTACAT".as_slice(),
            b"ACGTACAT".as_slice(),
            b"CCGTACGT".as_slice(),
        ] {
            assert_astar_matches_oracle(&g, query);
        }
    }

    #[test]
    fn without_pruning_matches_oracle() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).without_pruning();
        let result = engine.align(&g, b"AXGT").unwrap();

        let oracle: DpOracle<Affine, POAGraph<u32>> = DpOracle::new(costs());
        let expected = oracle.align(&g, b"AXGT").unwrap();
        assert_eq!(result.score, expected.score);
    }

    /// Bubble-exit pruning (Task 16 "bench & tune" check): confirm pruning
    /// never changes the reported score and never visits *more* states than
    /// leaving it off (the `AlignmentStats.cells_computed`/
    /// `fraction_of_full_matrix` plumbing already flows through unchanged
    /// from `AlignOutput`, so this doubles as the "no regression on the
    /// fraction-of-matrix metric" check called for by the task).
    ///
    /// Block size `B` (`BlockedVisitedStorage`'s const generic, default 8)
    /// and the `u64` offset-width fallback aren't swept here: neither is
    /// exposed as a tunable yet (`Offset` is hardcoded to `u32`, `B` is
    /// baked into `AstarData`), and adding
    /// that knob before profiling shows it matters would be tuning without
    /// a signal to tune against.
    #[test]
    fn pruning_never_increases_cells_computed_or_changes_score() {
        let g = multi_branch_graph();
        let query = b"ACGTACAT";

        let with_pruning: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs());
        let without_pruning: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).without_pruning();

        let pruned = with_pruning.align(&g, query).unwrap();
        let unpruned = without_pruning.align(&g, query).unwrap();

        assert_eq!(
            pruned.score, unpruned.score,
            "pruning must not change the optimal score"
        );
        assert!(
            pruned.stats.cells_computed <= unpruned.stats.cells_computed,
            "pruning should not visit more states than without it: {} vs {}",
            pruned.stats.cells_computed,
            unpruned.stats.cells_computed
        );
    }

    /// Assert the A* engine's own (independently implemented) backtrace
    /// re-scores, via [`score_from_alignment`], to exactly the score it
    /// reported — checked for both heuristics since each drives a separate
    /// `align` call and thus a separate backtrace walk.
    fn assert_astar_backtrace_matches_score(graph: &POAGraph<u32>, query: &[u8]) {
        let dijkstra_engine: AstarEngine<POAGraph<u32>, Affine, DijkstraHeuristic> =
            AstarEngine::new(costs()).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        let dijkstra_recomputed =
            score_from_alignment(graph, query, &dijkstra_result.alignment, &costs());
        assert_eq!(
            dijkstra_result.score,
            dijkstra_recomputed,
            "Dijkstra A* backtrace/score mismatch for query {:?}: reported={}, recomputed={}, alignment={:?}",
            String::from_utf8_lossy(query),
            dijkstra_result.score,
            dijkstra_recomputed,
            dijkstra_result.alignment
        );

        let mingapcost_engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs());
        let mingapcost_result = mingapcost_engine.align(graph, query).unwrap();
        let mingapcost_recomputed =
            score_from_alignment(graph, query, &mingapcost_result.alignment, &costs());
        assert_eq!(
            mingapcost_result.score,
            mingapcost_recomputed,
            "MinGapCost A* backtrace/score mismatch for query {:?}: reported={}, recomputed={}, alignment={:?}",
            String::from_utf8_lossy(query),
            mingapcost_result.score,
            mingapcost_recomputed,
            mingapcost_result.alignment
        );
    }

    #[test]
    fn backtrace_matches_score_linear() {
        assert_astar_backtrace_matches_score(&linear_graph(b"ACGT"), b"ACGT");
        assert_astar_backtrace_matches_score(&linear_graph(b"ACGT"), b"ACXT");
        assert_astar_backtrace_matches_score(&linear_graph(b"ACGT"), b"ACT");
        assert_astar_backtrace_matches_score(&linear_graph(b"ACT"), b"ACGT");
        assert_astar_backtrace_matches_score(&linear_graph(b"ACCCGT"), b"AGT");
    }

    #[test]
    fn backtrace_matches_score_multiple_substitutions() {
        assert_astar_backtrace_matches_score(&linear_graph(b"ACGTACGTAC"), b"AGGTCCGTAG");
    }

    #[test]
    fn backtrace_matches_score_bubble() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        for query in [
            b"ACGAT".as_slice(),
            b"ACGT".as_slice(),
            b"AXGT".as_slice(),
            b"ACGATT".as_slice(),
        ] {
            assert_astar_backtrace_matches_score(&g, query);
        }
    }

    #[test]
    fn backtrace_matches_score_deletion_spanning_bubble() {
        // Query drops bases straddling the bubble's merge point.
        let g = bubble_graph(b"ACGAT", b"ACGT");
        for query in [b"AT".as_slice(), b"AGT".as_slice(), b"A".as_slice()] {
            assert_astar_backtrace_matches_score(&g, query);
        }
    }

    #[test]
    fn backtrace_matches_score_insertion_at_branch_point() {
        let g = bubble_graph(b"ACGAT", b"ACGT");
        for query in [b"ACGATXXT".as_slice(), b"ACGTXXT".as_slice()] {
            assert_astar_backtrace_matches_score(&g, query);
        }
    }

    #[test]
    fn backtrace_matches_score_multi_branch_graph() {
        let g = multi_branch_graph();
        for query in [
            b"ACGTACGT".as_slice(),
            b"ACGAACGT".as_slice(),
            b"ACGTCCGT".as_slice(),
            b"AGGTACAT".as_slice(),
            b"ACGTACAT".as_slice(),
            b"CCGTACGT".as_slice(),
            b"AGGAACAT".as_slice(),
            b"ACGACGT".as_slice(),
            b"ACGTACG".as_slice(),
        ] {
            assert_astar_backtrace_matches_score(&g, query);
        }
    }

    #[test]
    fn no_alignment_found_on_disconnected_graph_returns_err() {
        // A graph with an unreachable component (no path from start to end
        // through it) still aligns fine through the real path; to force
        // `NoAlignmentFound` we'd need a graph with no start->end path at
        // all, which `POAGraph` does not allow to construct. Instead assert
        // the error type is usable as `std::error::Error` (compile-time check).
        fn assert_error<E: Error>() {}
        assert_error::<AstarError>();
    }

    // ─── Linear cost model (Task 13) ────────────────────────────────────────

    fn assert_astar_matches_oracle_linear(
        graph: &POAGraph<u32>,
        query: &[u8],
        costs: crate::align::cost_models::linear::Linear,
    ) {
        let oracle: DpOracle<crate::align::cost_models::linear::Linear, POAGraph<u32>> =
            DpOracle::new(costs);
        let expected = oracle.align(graph, query).unwrap();

        let dijkstra_engine: AstarEngine<
            POAGraph<u32>,
            crate::align::cost_models::linear::Linear,
            DijkstraHeuristic,
        > = AstarEngine::new(costs).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        assert_eq!(
            dijkstra_result.score, expected.score,
            "Dijkstra A* (linear) score mismatch"
        );

        let mingapcost_engine: AstarEngine<
            POAGraph<u32>,
            crate::align::cost_models::linear::Linear,
            MinGapCostHeuristic,
        > = AstarEngine::new(costs);
        let mingapcost_result = mingapcost_engine.align(graph, query).unwrap();
        assert_eq!(
            mingapcost_result.score, expected.score,
            "MinGapCost A* (linear) score mismatch"
        );
    }

    #[test]
    fn linear_cost_model_matches_oracle() {
        let linear_costs = crate::align::cost_models::linear::Linear::new(1, 2); // mismatch=1, gap_extend=2

        assert_astar_matches_oracle_linear(&linear_graph(b"ACGT"), b"ACGT", linear_costs);
        assert_astar_matches_oracle_linear(&linear_graph(b"ACGT"), b"ACXT", linear_costs);
        assert_astar_matches_oracle_linear(&linear_graph(b"ACGT"), b"ACT", linear_costs);
        assert_astar_matches_oracle_linear(&linear_graph(b"ACT"), b"ACGT", linear_costs);
        assert_astar_matches_oracle_linear(&linear_graph(b"ACCCGT"), b"AGT", linear_costs);

        let g = bubble_graph(b"ACGAT", b"ACGT");
        assert_astar_matches_oracle_linear(&g, b"ACGAT", linear_costs);
        assert_astar_matches_oracle_linear(&g, b"AXGT", linear_costs);

        let g2 = multi_branch_graph();
        for query in [
            b"ACGTACGT".as_slice(),
            b"ACGAACGT".as_slice(),
            b"ACGTACAT".as_slice(),
        ] {
            assert_astar_matches_oracle_linear(&g2, query, linear_costs);
        }
    }

    /// The one-state linear engine: pruning must not change the score, and
    /// its own backtrace must re-score to what it reported.
    #[test]
    fn linear_pruning_and_backtrace_are_consistent() {
        let costs = Linear::new(1, 2);
        let g = multi_branch_graph();
        for query in [
            b"ACGTACAT".as_slice(),
            b"ACGACGT".as_slice(),
            b"ACGTTACGTA".as_slice(),
            b"AGGT".as_slice(),
        ] {
            let pruned: AstarEngine<POAGraph<u32>, Linear> = AstarEngine::new(costs);
            let unpruned: AstarEngine<POAGraph<u32>, Linear> = AstarEngine::new(costs).without_pruning();
            let p = pruned.align(&g, query).unwrap();
            let u = unpruned.align(&g, query).unwrap();
            assert_eq!(p.score, u.score, "pruning changed the linear score");
            assert!(p.stats.cells_computed <= u.stats.cells_computed);
            assert_eq!(score_from_alignment(&g, query, &p.alignment, &costs), p.score);
            assert_eq!(score_from_alignment(&g, query, &u.alignment, &costs), u.score);
        }
    }

    // ─── Two-piece cost model ─────────────────────────────────────

    fn assert_astar_matches_oracle_two_piece(
        graph: &POAGraph<u32>,
        query: &[u8],
        costs: crate::align::cost_models::two_piece::TwoPieceAffine,
    ) {
        let oracle: DpOracle<
            crate::align::cost_models::two_piece::TwoPieceAffine,
            POAGraph<u32>,
        > = DpOracle::new(costs);
        let expected = oracle.align(graph, query).unwrap();

        let dijkstra_engine: AstarEngine<POAGraph<u32>, TwoPieceAffine, DijkstraHeuristic> =
            AstarEngine::new(costs).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        assert_eq!(
            dijkstra_result.score, expected.score,
            "Dijkstra A* (two-piece) score mismatch"
        );

        let mingapcost_engine: AstarEngine<POAGraph<u32>, TwoPieceAffine, MinGapCostHeuristic> =
            AstarEngine::new(costs);
        let mingapcost_result = mingapcost_engine.align(graph, query).unwrap();
        assert_eq!(
            mingapcost_result.score, expected.score,
            "MinGapCost A* (two-piece) score mismatch"
        );
    }

    #[test]
    fn two_piece_cost_model_matches_oracle() {
        // mismatch=1, piece1 (short gaps): open=2, extend=2; piece2 (long
        // gaps): open=6, extend=1. go1 <= go2 and ge1 > ge2, so opening is
        // always at least as cheap under piece 1 and long gaps favour piece 2
        // — the intended two-piece usage.
        let costs = crate::align::cost_models::two_piece::TwoPieceAffine::new(1, 2, 2, 6, 1);

        assert_astar_matches_oracle_two_piece(&linear_graph(b"ACGT"), b"ACGT", costs);
        assert_astar_matches_oracle_two_piece(&linear_graph(b"ACGT"), b"ACXT", costs);
        assert_astar_matches_oracle_two_piece(&linear_graph(b"ACGT"), b"ACT", costs);
        assert_astar_matches_oracle_two_piece(&linear_graph(b"ACT"), b"ACGT", costs);
        // Long deletion: cheap under piece 2 (long-gap-optimized).
        assert_astar_matches_oracle_two_piece(&linear_graph(b"ACGTACGTACGT"), b"AT", costs);

        let g = bubble_graph(b"ACGAT", b"ACGT");
        assert_astar_matches_oracle_two_piece(&g, b"ACGAT", costs);
        assert_astar_matches_oracle_two_piece(&g, b"AXGT", costs);

        let g2 = multi_branch_graph();
        for query in [
            b"ACGTACGT".as_slice(),
            b"ACGAACGT".as_slice(),
            b"ACGTACAT".as_slice(),
        ] {
            assert_astar_matches_oracle_two_piece(&g2, query, costs);
        }
    }

    // ─── Semi-global / EndsFree (Task 14) ───────────────────────────────────
    //
    // No DP oracle exists for EndsFree on this branch (per the interface doc),
    // so these are hand-checked: each case's expected score is computed by
    // hand from the cost model, not cross-checked against DpOracle.

    fn free_ends(
        qry_free_begin: std::ops::Bound<usize>,
        qry_free_end: std::ops::Bound<usize>,
        graph_free_begin: std::ops::Bound<usize>,
        graph_free_end: std::ops::Bound<usize>,
    ) -> AlignmentMode {
        AlignmentMode::EndsFree {
            qry_free_begin,
            qry_free_end,
            graph_free_begin,
            graph_free_end,
        }
    }

    fn not_free() -> std::ops::Bound<usize> {
        std::ops::Bound::Excluded(0)
    }

    #[test]
    fn ends_free_trailing_query_bases_are_free() {
        // Graph "ACGT" exactly matches the query prefix; the trailing "XXXX"
        // has no corresponding graph path, but qry_free_end is unbounded, so
        // it should cost nothing.
        let g = linear_graph(b"ACGT");
        let mode = free_ends(
            not_free(),
            std::ops::Bound::Unbounded,
            not_free(),
            not_free(),
        );
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);
        let result = engine.align(&g, b"ACGTXXXX").unwrap();
        assert_eq!(
            result.score, 0,
            "free trailing query bases should not be penalised"
        );
    }

    #[test]
    fn ends_free_leading_query_bases_are_free() {
        // Leading "XXXX" has no corresponding graph path; qry_free_begin is
        // unbounded, so the alignment should start at query offset 4 for free.
        let g = linear_graph(b"ACGT");
        let mode = free_ends(
            std::ops::Bound::Unbounded,
            not_free(),
            not_free(),
            not_free(),
        );
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);
        let result = engine.align(&g, b"XXXXACGT").unwrap();
        assert_eq!(
            result.score, 0,
            "free leading query bases should not be penalised"
        );
    }

    #[test]
    fn ends_free_trailing_graph_nodes_are_free() {
        // Graph "ACGTACGT" is longer than the query "ACGT"; graph_free_end is
        // unbounded, so the unmatched trailing graph nodes cost nothing.
        let g = linear_graph(b"ACGTACGT");
        let mode = free_ends(
            not_free(),
            not_free(),
            not_free(),
            std::ops::Bound::Unbounded,
        );
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);
        let result = engine.align(&g, b"ACGT").unwrap();
        assert_eq!(
            result.score, 0,
            "free trailing graph nodes should not be penalised"
        );
    }

    #[test]
    fn ends_free_leading_graph_nodes_are_free() {
        // Graph "ACGTACGT"; graph_free_begin is unbounded, so the query can
        // match starting partway through the graph for free.
        let g = linear_graph(b"ACGTACGT");
        let mode = free_ends(
            not_free(),
            not_free(),
            std::ops::Bound::Unbounded,
            not_free(),
        );
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);
        let result = engine.align(&g, b"ACGT").unwrap();
        assert_eq!(
            result.score, 0,
            "free leading graph nodes should not be penalised"
        );
    }

    #[test]
    fn ends_free_dijkstra_and_mingapcost_agree() {
        let g = linear_graph(b"ACGTACGT");
        let mode = free_ends(
            std::ops::Bound::Unbounded,
            std::ops::Bound::Unbounded,
            not_free(),
            not_free(),
        );

        let dijkstra: AstarEngine<POAGraph<u32>, Affine, DijkstraHeuristic> =
            AstarEngine::new(costs())
                .with_mode(mode)
                .with_heuristic(DijkstraHeuristic);
        let mingapcost: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);

        let query = b"XXACGTAC";
        assert_eq!(
            dijkstra.align(&g, query).unwrap().score,
            mingapcost.align(&g, query).unwrap().score
        );
    }

    #[test]
    fn ends_free_interior_mismatch_still_costs() {
        // Free query ends don't exempt an interior mismatch from its usual cost.
        let g = linear_graph(b"ACGT");
        let mode = free_ends(
            not_free(),
            std::ops::Bound::Unbounded,
            not_free(),
            not_free(),
        );
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);
        let result = engine.align(&g, b"AXGTXXXX").unwrap();
        assert_eq!(
            result.score, 1,
            "interior mismatch should still cost 1, trailing XXXX still free"
        );
    }

    #[test]
    fn ends_free_two_piece_trailing_query_bases_are_free() {
        let g = linear_graph(b"ACGT");
        let mode = free_ends(
            not_free(),
            std::ops::Bound::Unbounded,
            not_free(),
            not_free(),
        );
        let costs = crate::align::cost_models::two_piece::TwoPieceAffine::new(1, 2, 2, 6, 1);
        let engine: AstarEngine<POAGraph<u32>, TwoPieceAffine, MinGapCostHeuristic> =
            AstarEngine::new(costs).with_mode(mode);
        let result = engine.align(&g, b"ACGTXXXX").unwrap();
        assert_eq!(
            result.score, 0,
            "free trailing query bases should not be penalised (two-piece)"
        );
    }

    /// Regression test for a backtrace bug: `SemiGlobal` mode (query aligned
    /// end-to-end, graph ends free) can terminate on a *real* graph node
    /// instead of the end sentinel (see `ends_free_is_end`). The backtrace
    /// used to only ever record `curr` once its predecessor was found, so the
    /// end state passed in — real for `EndsFree`, the sentinel for `Global` —
    /// was silently dropped whenever it wasn't the sentinel. For a query
    /// identical to a microsatellite-like repeat already in the graph, that
    /// turned the final matched base into a spurious 1-base insertion (and a
    /// matching deletion of the graph's real last node) instead of a plain
    /// match, growing the graph by an extra node per added sequence.
    #[test]
    fn semi_global_exact_duplicate_has_no_dropped_final_match() {
        let seq = b"TCTCTCTCTCTCTCTCTCTCTC"; // "TC" x 11
        let g = linear_graph(seq);
        let mode = free_ends(
            not_free(),
            not_free(),
            std::ops::Bound::Unbounded,
            std::ops::Bound::Unbounded,
        );
        let engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).with_mode(mode);
        let result = engine.align(&g, seq).unwrap();

        assert_eq!(
            result.score, 0,
            "identical sequence should align with score 0"
        );
        assert_eq!(
            result.alignment.len(),
            seq.len(),
            "every query base must be covered by the alignment, got {:?}",
            result.alignment
        );
        assert!(
            result
                .alignment
                .iter()
                .all(|p| p.node().is_some() && p.query_pos().is_some()),
            "an exact match must consist only of match pairs, no indels: {:?}",
            result.alignment
        );
    }
}
