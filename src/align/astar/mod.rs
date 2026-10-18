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
mod differential_tests;
#[cfg(test)]
pub(crate) mod dp_oracle;

use std::error::Error;
use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::align::astar::aln_graph::{
    AffineAlignmentGraph, AlignmentGraph, AlignmentGraphNode, EndsFreeRegion, LinearAlignmentGraph,
    TwoPieceAlignmentGraph,
};
use crate::align::astar::dfa::{DepthFirstGreedyAlignment, ExtendResult};
use crate::align::astar::heuristic::AstarHeuristic;
pub use crate::align::astar::heuristic::{
    DijkstraHeuristic, HeuristicConfig, MinGapCostHeuristic, SeedParams, SmoothHeuristicConfig,
};
use crate::align::astar::queue::{
    AffineLayeredQueue, AstarQueue, LinearLayeredQueue, TwoPieceLayeredQueue,
};
use crate::align::astar::visited::{
    AffineAstarData, AstarData, AstarVisited, LinearAstarData, TwoPieceAstarData,
};
use crate::align::cost_models::AlignmentCostModel;
use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::offsets::OffsetType;
use crate::align::traits::{AlignOutput, AlignmentStats};
use crate::align::traits::{AlignableGraph, AlignmentEngine};
use crate::align::{AlignState, AlignmentMode};
use crate::graph::bubbles::index::BubbleIndex;
use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};

/// Default query offset width, the one the CLI uses. Library callers pick
/// another with `AstarEngine`'s `O` parameter: `u16` for short reads, `u64`
/// for queries over 4 Gbp.
pub type Offset = u32;

/// Returned when the priority queue empties before the end state is reached
/// (e.g. a disconnected graph).
#[derive(Debug)]
pub enum AstarError {
    NoAlignmentFound,
    /// Returned when the query is longer than `O::MAX`, the largest offset
    /// the engine's offset type holds.
    QueryTooLong {
        len: usize,
        max: usize,
    },
}

impl fmt::Display for AstarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoAlignmentFound => write!(
                f,
                "A* queue emptied before reaching the alignment end state"
            ),
            Self::QueryTooLong { len, max } => write!(
                f,
                "query of length {len} exceeds the offset type's maximum of {max}"
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
/// heuristic; `O` the query offset type, which caps the query length.
pub struct AstarEngine<G, C = Affine, H = MinGapCostHeuristic, O = Offset>
where
    G: AlignableGraph,
    O: OffsetType,
    C: AstarCostModel<G::Node, O>,
    H: HeuristicConfig,
{
    costs: C,
    mode: AlignmentMode,
    heuristic_cfg: H,
    /// Built from `heuristic_cfg` on the first `align` and kept, like
    /// `visited`, so a heuristic keeps its buffers between alignments.
    heuristic: Option<H::Heuristic<G::Node, O, C>>,
    existing_bubbles: Option<Arc<BubbleIndex<G::Node>>>,
    enable_pruning: bool,
    /// Count re-expanded states in `AlignmentStats::states_reopened`.
    track_reopen: bool,
    /// Visited storage and queue kept across `align` calls. Per alignment
    /// they would otherwise allocate one block window per touched rank group
    /// and a queue layer per priority, and free it all at the end. That free
    /// storm was 5 to 8 % of the run on inputs with many short alignments.
    /// Reset walks only the rows the previous alignment touched, so a short
    /// query after a long one does not pay for the long one's rows. So
    /// `align` takes `&mut self`.
    visited: Option<C::Visited>,
    queue: C::Queue,
    _graph: PhantomData<(G, O)>,
}

impl<G, C, O> AstarEngine<G, C, MinGapCostHeuristic, O>
where
    G: AlignableGraph,
    O: OffsetType,
    C: AstarCostModel<G::Node, O>,
{
    /// Global alignment, min-gap-cost heuristic, pruning enabled.
    pub fn new(costs: C) -> Self {
        Self {
            costs,
            mode: AlignmentMode::Global,
            heuristic_cfg: MinGapCostHeuristic,
            heuristic: None,
            existing_bubbles: None,
            enable_pruning: true,
            track_reopen: false,
            visited: None,
            queue: C::Queue::default(),
            _graph: PhantomData,
        }
    }
}

impl<G, C, H, O> AstarEngine<G, C, H, O>
where
    G: AlignableGraph,
    O: OffsetType,
    C: AstarCostModel<G::Node, O>,
    H: HeuristicConfig,
{
    pub fn with_mode(mut self, mode: AlignmentMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn with_heuristic<H2: HeuristicConfig>(
        self,
        heuristic_cfg: H2,
    ) -> AstarEngine<G, C, H2, O> {
        AstarEngine {
            costs: self.costs,
            mode: self.mode,
            heuristic_cfg,
            heuristic: None,
            existing_bubbles: self.existing_bubbles,
            enable_pruning: self.enable_pruning,
            track_reopen: self.track_reopen,
            visited: self.visited,
            queue: self.queue,
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

    /// Count re-expanded states in `AlignmentStats::states_reopened`. Off by
    /// default, because it keeps a hash set of every expanded state.
    pub fn with_reopen_tracking(mut self) -> Self {
        self.track_reopen = true;
        self
    }
}

impl<G, C, H, O> AlignmentEngine<&[u8]> for AstarEngine<G, C, H, O>
where
    G: AlignableGraph,
    G::Node: 'static,
    O: OffsetType,
    C: AstarCostModel<G::Node, O>,
    H: HeuristicConfig,
{
    type Graph = G;
    type Success = AlignOutput<G>;
    type Error = AstarError;

    fn align(&mut self, graph: &G, query: &[u8]) -> Result<AlignOutput<G>, Self::Error> {
        let _span = tracing::info_span!(
            "align_sequence",
            engine = "astar",
            query_len = query.len(),
            graph_nodes = graph.node_count(),
        )
        .entered();

        if query.len() > O::MAX {
            return Err(AstarError::QueryTooLong {
                len: query.len(),
                max: O::MAX,
            });
        }

        let bubble_index = self
            .existing_bubbles
            .clone()
            .unwrap_or_else(|| Arc::new(BubbleIndex::new(graph)));

        let aln_graph = self.costs.aln_graph(self.mode);
        let costs = self.costs;
        let heuristic_cfg = &self.heuristic_cfg;
        let heuristic = self
            .heuristic
            .get_or_insert_with(|| heuristic_cfg.instantiate(costs));
        heuristic.prepare(graph, query, &bubble_index, self.mode);

        // `costs` is copied out so the closures do not borrow `self` while
        // `self.visited` and `self.heuristic` are borrowed. The first `align`
        // resets a fresh storage, one pass over the nodes, once per engine.
        let visited = self
            .visited
            .get_or_insert_with(|| costs.visited(graph, query, bubble_index.clone()));
        visited.reset(graph, query, bubble_index);
        self.queue.clear();

        let result = run_astar_search(
            graph,
            query,
            &aln_graph,
            visited,
            &mut self.queue,
            heuristic,
            self.enable_pruning,
            self.track_reopen,
        )?;

        Ok(finish_alignment(graph, query, visited, result, C::STATES))
    }
}

/// The A* main loop, shared across every cost-model shape: best-first search
/// with greedy zero-cost match extension via [`DepthFirstGreedyAlignment`],
/// re-entering the queue only at mismatches, the graph end, or the query end.
/// Every argument is behind a trait (`AlignmentGraph`, `AstarVisited`,
/// `AstarQueue`), so this one function drives the 1-, 3- and 5-state
/// machines alike.
#[allow(clippy::too_many_arguments)]
fn run_astar_search<G, O, AG, V, Q, Hs>(
    graph: &G,
    query: &[u8],
    aln_graph: &AG,
    visited: &mut V,
    queue: &mut Q,
    heuristic: &mut Hs,
    enable_pruning: bool,
    track_reopen: bool,
) -> Result<SearchResult<G::Node, O>, AstarError>
where
    G: AlignableGraph,
    O: OffsetType,
    AG: AlignmentGraph,
    V: AstarVisited<G::Node, O>,
    Q: AstarQueue<G::Node, O>,
    Hs: AstarHeuristic<G::Node, O>,
{
    let mut num_visited = 0usize;
    let mut num_expanded = 0usize;
    let mut num_pruned = 0usize;
    // A hash set of every expanded state, only when tracking is
    // on. For tests and diagnosis, far too slow for production runs.
    let mut expanded: Option<rustc_hash::FxHashSet<(usize, O, u8)>> =
        track_reopen.then(Default::default);
    let mut num_reopened = 0usize;

    // Start states get their score in the visited storage when popped, not
    // here. In semi-global mode every node is a start state at offset 0, and
    // writing them all up front stretches each block row's offset window
    // from 0 to wherever the search goes later: gigabytes on a 50 kb read,
    // for states a good heuristic never expands.
    let region = EndsFreeRegion::new(graph, query, aln_graph.mode());
    for initial_state in aln_graph.initial_states::<G, O>(graph, &region) {
        let h = heuristic.h(graph, &initial_state, AlignState::Match);
        queue.queue_aln_state(initial_state, AlignState::Match, 0, h);
    }

    // One stack allocation per alignment. The stack borrows the graph through
    // its successor iterators, so keeping it in the engine would put a graph
    // lifetime on the engine. Not worth it for one malloc.
    let mut dfa = DepthFirstGreedyAlignment::new(graph, query, enable_pruning);
    let dynamic = heuristic.is_dynamic();
    if track_reopen || dynamic {
        dfa.matched_log = Some(Vec::new());
    }
    let (end_score, end_node) = 'main: loop {
        let priority = queue.min_priority();
        let Some(item) = queue.pop_aln_state() else {
            return Err(AstarError::NoAlignmentFound);
        };
        let score = item.score();
        let aln_node = item.aln_node();
        let aln_state = item.aln_state();

        let stored = visited.get_score(&aln_node, aln_state);
        if score > stored {
            continue;
        }
        if score < stored {
            visited.set_score(&aln_node, aln_state, score);
        }

        // Match pruning raises h behind the search front. A state queued
        // before that goes back in at its new priority.
        if dynamic {
            let h = heuristic.h(graph, &aln_node, aln_state);
            if score as usize + h > priority {
                tracing::trace!(
                    target: "poasta::astar::requeued",
                    score,
                    node = aln_node.node().index(),
                    offset = aln_node.offset().as_usize(),
                    state = ?aln_state,
                    priority = score as usize + h,
                );
                queue.queue_aln_state(aln_node, aln_state, score, h);
                continue;
            }
        }

        if aln_graph.is_end(graph, query, &region, &aln_node, aln_state) {
            trace_visited(score, &aln_node, aln_state);
            num_visited += 1;
            break (score, aln_node);
        }

        if enable_pruning && visited.prune(score, &aln_node, aln_state) {
            trace_pruned(score, &aln_node, aln_state);
            num_pruned += 1;
            continue;
        }

        if let Some(set) = &mut expanded
            && !set.insert((aln_node.node().index(), aln_node.offset(), aln_state as u8))
        {
            num_reopened += 1;
        }

        visited.mark_reached(score, &aln_node, aln_state);
        trace_visited(score, &aln_node, aln_state);
        num_visited += 1;
        num_expanded += 1;

        if aln_state == AlignState::Match {
            heuristic.on_expand(graph, &aln_node, usize::MAX);
            dfa.start(score, &aln_node);
            let is_end_check = |n: &AlignmentGraphNode<G::Node, O>| {
                aln_graph.is_end(graph, query, &region, n, AlignState::Match)
            };

            while let Some(end_point) = dfa.extend(visited, is_end_check) {
                match end_point {
                    ExtendResult::RefGraphEnd(parent, child) => {
                        if aln_graph.is_end(graph, query, &region, &child, AlignState::Match) {
                            num_visited += dfa.get_num_visited();
                            break 'main (score, child);
                        }

                        aln_graph.expand_ref_graph_end(
                            visited,
                            &parent,
                            score,
                            |delta, succ, succ_state| {
                                let h = heuristic.h(graph, &succ, succ_state);
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
                                let h = heuristic.h(graph, &succ, succ_state);
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
                                let h = heuristic.h(graph, &succ, succ_state);
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
            if let Some(log) = &mut dfa.matched_log {
                for n in log.drain(..) {
                    // A passed state whose f still equals the popped priority
                    // could have been the next pop, so it is at its optimal
                    // score too and may prune (proofs note, Theorem 3).
                    if dynamic {
                        heuristic.on_expand(graph, &n, priority - score as usize);
                    }
                    if let Some(set) = &mut expanded
                        && !set.insert((n.node().index(), n.offset(), AlignState::Match as u8))
                    {
                        num_reopened += 1;
                    }
                }
            }
        } else {
            aln_graph.expand_all(
                visited,
                graph,
                query,
                score,
                &aln_node,
                aln_state,
                |delta, succ, succ_state| {
                    let h = heuristic.h(graph, &succ, succ_state);
                    queue.queue_aln_state(succ, succ_state, score + delta as u32, h);
                },
            );
        }
    };

    Ok(SearchResult {
        end_score,
        end_node,
        num_visited,
        num_expanded,
        num_pruned: num_pruned + dfa.num_pruned,
        num_reopened,
    })
}

/// What `run_astar_search` hands to `finish_alignment`.
struct SearchResult<N: GraphNodeId, O> {
    end_score: u32,
    end_node: AlignmentGraphNode<N, O>,
    num_visited: usize,
    num_expanded: usize,
    num_pruned: usize,
    num_reopened: usize,
}

/// Emit the `set_visited` trace event consumed by `debug::subscriber`.
#[inline]
fn trace_visited<N: GraphNodeId, O: OffsetType>(
    score: u32,
    node: &AlignmentGraphNode<N, O>,
    state: AlignState,
) {
    tracing::trace!(
        target: "poasta::astar::set_visited",
        score,
        node = node.node().index(),
        offset = node.offset().as_usize(),
        state = ?state,
    );
}

/// Emit the `pruned` trace event: a state dropped by bubble-exit pruning.
#[inline]
pub(crate) fn trace_pruned<N: GraphNodeId, O: OffsetType>(
    score: u32,
    node: &AlignmentGraphNode<N, O>,
    state: AlignState,
) {
    tracing::trace!(
        target: "poasta::astar::pruned",
        score,
        node = node.node().index(),
        offset = node.offset().as_usize(),
        state = ?state,
    );
}

/// Backtrace + wrap the search result into an [`AlignOutput`], shared by
/// every cost model. `n_states` is the visited-storage state count, used only
/// for the fraction-of-full-matrix stat.
fn finish_alignment<G, O, V>(
    graph: &G,
    query: &[u8],
    visited: &V,
    result: SearchResult<G::Node, O>,
    n_states: usize,
) -> AlignOutput<G>
where
    G: AlignableGraph,
    O: OffsetType,
    V: AstarVisited<G::Node, O>,
{
    let SearchResult {
        end_score,
        end_node,
        num_visited,
        num_expanded,
        num_pruned,
        num_reopened,
    } = result;
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
            cells_computed: num_visited,
            states_expanded: num_expanded,
            fraction_of_full_matrix,
            states_pruned: num_pruned,
            states_reopened: num_reopened,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::dp_oracle::DpOracle;
    use super::*;
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

    /// Assert the A* engine (both heuristics) agrees with the DP oracle's
    /// score on `graph`/`query`.
    fn assert_astar_matches_oracle(graph: &POAGraph<u32>, query: &[u8]) {
        let mut oracle: DpOracle<Affine, POAGraph<u32>> = DpOracle::new(costs());
        let expected = oracle.align(graph, query).unwrap();

        let mut dijkstra_engine: AstarEngine<POAGraph<u32>, Affine, DijkstraHeuristic> =
            AstarEngine::new(costs()).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        assert_eq!(
            dijkstra_result.score,
            expected.score,
            "Dijkstra A* score mismatch for query {:?}",
            String::from_utf8_lossy(query)
        );

        let mut mingapcost_engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs()).without_pruning();
        let result = engine.align(&g, b"AXGT").unwrap();

        let mut oracle: DpOracle<Affine, POAGraph<u32>> = DpOracle::new(costs());
        let expected = oracle.align(&g, b"AXGT").unwrap();
        assert_eq!(result.score, expected.score);
    }

    /// Bubble-exit pruning must leave the reported score unchanged and must
    /// not visit more states than leaving it off. The `cells_computed` and
    /// `fraction_of_full_matrix` statistics come through unchanged from
    /// `AlignOutput`, so this also guards the fraction-of-matrix metric.
    ///
    /// Block size `B` (`BlockStorage`'s const generic, default 8) isn't
    /// swept here: it is baked into `AstarData`, and adding that knob before
    /// profiling shows it matters would be tuning without a signal to tune
    /// against. Offset widths are covered by `offset_widths_agree`.
    #[test]
    fn pruning_never_increases_cells_computed_or_changes_score() {
        let g = multi_branch_graph();
        let query = b"ACGTACAT";

        let mut with_pruning: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
            AstarEngine::new(costs());
        let mut without_pruning: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        assert_eq!(unpruned.stats.states_pruned, 0);
    }

    /// `without_pruning` must switch off the check in greedy match
    /// extension too, not only at queue pops. A graph of real sequences has
    /// enough bubbles for the default engine to prune.
    #[test]
    fn without_pruning_prunes_nothing() {
        use crate::align::PoastaAligner;
        use crate::graph::io::seq::open_sequences;
        let seqs: Vec<Vec<u8>> = open_sequences(std::path::Path::new("tests/dnaG/align_seq.fna"))
            .unwrap()
            .take(6)
            .map(|r| r.unwrap().1)
            .collect();
        let mut aligner = PoastaAligner::new(
            AstarEngine::<POAGraph<u32>, Affine>::new(costs()),
            POAGraph::new(),
        );
        for (i, s) in seqs[..5].iter().enumerate() {
            aligner.align_one_named(&format!("s{i}"), s).unwrap();
        }
        let g = aligner.into_graph();

        let pruned = AstarEngine::<POAGraph<u32>, Affine>::new(costs())
            .align(&g, &seqs[5])
            .unwrap();
        let unpruned = AstarEngine::<POAGraph<u32>, Affine>::new(costs())
            .without_pruning()
            .align(&g, &seqs[5])
            .unwrap();
        assert_eq!(pruned.score, unpruned.score);
        assert!(pruned.stats.states_pruned > 0);
        assert_eq!(unpruned.stats.states_pruned, 0);
    }

    /// Assert the A* engine's own (independently implemented) backtrace
    /// re-scores, via [`score_from_alignment`], to exactly the score it
    /// reported. Checked for both heuristics, since each is a separate
    /// `align` call and so a separate backtrace walk.
    fn assert_astar_backtrace_matches_score(graph: &POAGraph<u32>, query: &[u8]) {
        let mut dijkstra_engine: AstarEngine<POAGraph<u32>, Affine, DijkstraHeuristic> =
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

        let mut mingapcost_engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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

    /// One engine reused across queries and graphs must give the same score
    /// and do the same work as a fresh engine per query. A reset that leaves
    /// a stale cell or reached bit shows up as a different score or count.
    #[test]
    fn engine_reuse_gives_same_scores() {
        let multi = multi_branch_graph();
        let small = linear_graph(b"ACGTTGCA");
        let mut reused: AstarEngine<POAGraph<u32>, Affine> = AstarEngine::new(costs());

        let cases = [
            (&multi, b"ACGTACGT".as_slice()),
            (&multi, b"ACGAACGT".as_slice()),
            (&multi, b"ACGTCCGT".as_slice()),
            (&multi, b"AGGTACAT".as_slice()),
            (&multi, b"ACGTACAT".as_slice()),
            (&multi, b"CCGTACGT".as_slice()),
            // A long query, then a short one, then a different graph.
            (&multi, b"ACGTTACGTAACGTTTACGTACGGT".as_slice()),
            (&multi, b"AC".as_slice()),
            (&small, b"ACGATGCA".as_slice()),
            (&multi, b"ACGTACAT".as_slice()),
        ];
        for (graph, query) in cases {
            let fresh = AstarEngine::<POAGraph<u32>, Affine>::new(costs())
                .align(graph, query)
                .unwrap();
            let again = reused.align(graph, query).unwrap();
            assert_eq!(
                (again.score, again.stats.cells_computed),
                (fresh.score, fresh.stats.cells_computed),
                "reused engine differs for query {:?}",
                String::from_utf8_lossy(query)
            );
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

    // ─── Linear cost model ───────────────────────────────────────────────────────────────────────────────────────

    fn assert_astar_matches_oracle_linear(
        graph: &POAGraph<u32>,
        query: &[u8],
        costs: crate::align::cost_models::linear::Linear,
    ) {
        let mut oracle: DpOracle<crate::align::cost_models::linear::Linear, POAGraph<u32>> =
            DpOracle::new(costs);
        let expected = oracle.align(graph, query).unwrap();

        let mut dijkstra_engine: AstarEngine<
            POAGraph<u32>,
            crate::align::cost_models::linear::Linear,
            DijkstraHeuristic,
        > = AstarEngine::new(costs).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        assert_eq!(
            dijkstra_result.score, expected.score,
            "Dijkstra A* (linear) score mismatch"
        );

        let mut mingapcost_engine: AstarEngine<
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

    /// For the one-state linear engine, pruning must not change the score, and
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
            let mut pruned: AstarEngine<POAGraph<u32>, Linear> = AstarEngine::new(costs);
            let mut unpruned: AstarEngine<POAGraph<u32>, Linear> =
                AstarEngine::new(costs).without_pruning();
            let p = pruned.align(&g, query).unwrap();
            let u = unpruned.align(&g, query).unwrap();
            assert_eq!(p.score, u.score, "pruning changed the linear score");
            assert!(p.stats.cells_computed <= u.stats.cells_computed);
            assert_eq!(
                score_from_alignment(&g, query, &p.alignment, &costs),
                p.score
            );
            assert_eq!(
                score_from_alignment(&g, query, &u.alignment, &costs),
                u.score
            );
        }
    }

    // ─── Two-piece cost model ─────────────────────────────────────

    fn assert_astar_matches_oracle_two_piece(
        graph: &POAGraph<u32>,
        query: &[u8],
        costs: crate::align::cost_models::two_piece::TwoPieceAffine,
    ) {
        let mut oracle: DpOracle<
            crate::align::cost_models::two_piece::TwoPieceAffine,
            POAGraph<u32>,
        > = DpOracle::new(costs);
        let expected = oracle.align(graph, query).unwrap();

        let mut dijkstra_engine: AstarEngine<POAGraph<u32>, TwoPieceAffine, DijkstraHeuristic> =
            AstarEngine::new(costs).with_heuristic(DijkstraHeuristic);
        let dijkstra_result = dijkstra_engine.align(graph, query).unwrap();
        assert_eq!(
            dijkstra_result.score, expected.score,
            "Dijkstra A* (two-piece) score mismatch"
        );

        let mut mingapcost_engine: AstarEngine<POAGraph<u32>, TwoPieceAffine, MinGapCostHeuristic> =
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
        // always at least as cheap under piece 1 and long gaps favour piece 2,
        // the intended two-piece usage.
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

    // ─── Semi-global / EndsFree ──────────────────────────────────────────────
    //
    // The DP oracle is global only,
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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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

        let mut dijkstra: AstarEngine<POAGraph<u32>, Affine, DijkstraHeuristic> =
            AstarEngine::new(costs())
                .with_mode(mode)
                .with_heuristic(DijkstraHeuristic);
        let mut mingapcost: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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
        let mut engine: AstarEngine<POAGraph<u32>, TwoPieceAffine, MinGapCostHeuristic> =
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
    /// end state passed in (real for `EndsFree`, the sentinel for `Global`)
    /// was silently dropped whenever it wasn't the sentinel. For a query
    /// identical to a microsatellite-like repeat already in the graph, that
    /// turned the final matched base into a spurious 1-base insertion (and a
    /// matching deletion of the graph's real last node) instead of a plain
    /// match, growing the graph by an extra node per added sequence.
    /// Bubble pruning priced the gap from a reached exit with the cheaper
    /// gap piece, also for a candidate inside a gap of the other piece. That
    /// pruned the first state of this optimal alignment, which opens with
    /// five piece-two deletions.
    #[test]
    fn two_piece_pruning_prices_gaps_with_the_candidate_piece() {
        let costs = TwoPieceAffine::new(1, 2, 2, 6, 1);
        let g = linear_graph(b"AAAAACAAAACACCCAACAACACCCAACAAAACACAAA");
        let q = b"CAAAACCCCAACCAACACCC";
        let want = DpOracle::<TwoPieceAffine, POAGraph<u32>>::new(costs)
            .align(&g, q)
            .unwrap()
            .score;
        let got = AstarEngine::<POAGraph<u32>, TwoPieceAffine>::new(costs)
            .with_heuristic(DijkstraHeuristic)
            .align(&g, q)
            .unwrap()
            .score;
        assert_eq!(want, 33);
        assert_eq!(got, want);
    }

    /// A deletion candidate reaches its bubble exit by deleting, below the
    /// diagonal window, and then keeps deleting at its own rate. A path from
    /// a reached exit state would have to open the deletion again. The check
    /// ignored that and pruned the first state of this optimal alignment
    /// (min-gap-cost with pruning scored 95).
    #[test]
    fn deletion_candidates_are_only_dominated_after_a_reopen() {
        let costs = TwoPieceAffine::new(4, 6, 2, 24, 1);
        let g = linear_graph(b"AACAGTCCCATGACGTTCAGACGTGAAGCGCCTAGGGTAAATATACTGATTA");
        let q = b"ACCGTAGGCTAGA";
        let want = DpOracle::<TwoPieceAffine, POAGraph<u32>>::new(costs)
            .align(&g, q)
            .unwrap()
            .score;
        let got = AstarEngine::<POAGraph<u32>, TwoPieceAffine>::new(costs)
            .align(&g, q)
            .unwrap()
            .score;
        assert_eq!(want, 94);
        assert_eq!(got, want);
    }

    /// Insertions before the first graph node sit on the start sentinel at
    /// offsets above 0. The gap backtrace used to stop as soon as it reached
    /// the start node, so those query bases were left out of the alignment,
    /// and out of the MSA.
    #[test]
    fn backtrace_keeps_leading_insertions() {
        fn graph(seqs: &[&[u8]]) -> POAGraph<u32> {
            let engine = AstarEngine::<POAGraph<u32>, Affine>::new(costs());
            let mut aligner = crate::align::PoastaAligner::new(engine, POAGraph::new());
            for (i, s) in seqs.iter().enumerate() {
                aligner.align_one_named(&format!("s{i}"), s).unwrap();
            }
            aligner.into_graph()
        }
        fn covered(out: &AlignOutput<POAGraph<u32>>) -> usize {
            out.alignment
                .iter()
                .filter(|p| p.query_pos().is_some())
                .count()
        }

        let g = graph(&[b"TGCTACC", b"TGCCTACC", b"TGCCTACC", b"TGCCTACC"]);
        let q = b"CATTGAGAGCT";
        let out = AstarEngine::<POAGraph<u32>, Affine>::new(costs())
            .align(&g, q)
            .unwrap();
        assert_eq!(covered(&out), q.len(), "{:?}", out.alignment);
        assert_eq!(
            score_from_alignment(&g, q, &out.alignment, &costs()),
            out.score
        );

        let g = graph(&[b"CCACACCAA", b"CACAAACCA", b"CACCACCAA"]);
        let q = b"ACACCCCCC";
        let mode = free_ends(
            not_free(),
            not_free(),
            std::ops::Bound::Unbounded,
            std::ops::Bound::Unbounded,
        );
        let out = AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(4, 6, 2))
            .with_mode(mode)
            .align(&g, q)
            .unwrap();
        assert_eq!(covered(&out), q.len(), "{:?}", out.alignment);
    }

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
        let mut engine: AstarEngine<POAGraph<u32>, Affine, MinGapCostHeuristic> =
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

    /// `u16`, `u32` and `u64` offsets must give the same score and the same
    /// alignment, in both modes, and the backtrace must re-score to it.
    #[test]
    fn offset_widths_agree() {
        use std::ops::Bound::Unbounded;
        fn run<O: OffsetType>(g: &POAGraph<u32>, q: &[u8], mode: AlignmentMode) -> (u32, String) {
            let out = AstarEngine::<POAGraph<u32>, Affine, MinGapCostHeuristic, O>::new(costs())
                .with_mode(mode)
                .align(g, q)
                .unwrap();
            assert_eq!(
                out.score,
                score_from_alignment(g, q, &out.alignment, &costs())
            );
            (out.score, format!("{:?}", out.alignment))
        }

        let graphs = [
            linear_graph(b"ACGTACGT"),
            bubble_graph(b"ACGAT", b"ACGT"),
            multi_branch_graph(),
        ];
        for g in &graphs {
            for q in [&b"ACGTACGT"[..], b"AXGTCCGA", b"ACGAACGTT", b"CGT"] {
                for mode in [
                    AlignmentMode::Global,
                    free_ends(not_free(), not_free(), Unbounded, Unbounded),
                ] {
                    let want = run::<u32>(g, q, mode);
                    assert_eq!(run::<u16>(g, q, mode), want);
                    assert_eq!(run::<u64>(g, q, mode), want);
                }
            }
        }
    }

    #[test]
    fn query_longer_than_offset_type_is_rejected() {
        let g = linear_graph(b"ACGT");
        let q = vec![b'A'; u16::MAX as usize + 1];
        let err = AstarEngine::<POAGraph<u32>, Affine, MinGapCostHeuristic, u16>::new(costs())
            .align(&g, &q);
        assert!(matches!(err, Err(AstarError::QueryTooLong { len, max: 65535 }) if len == q.len()));
    }
}
