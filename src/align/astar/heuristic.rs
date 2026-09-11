//! A* heuristics, ported from `stable`'s `aligner::heuristic`.

use std::marker::PhantomData;
use std::sync::Arc;

use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::astar::aln_graph::{resolve_bound, AlignmentGraphNode};
use crate::align::offsets::OffsetType;
use crate::align::{AlignState, AlignmentMode};
use crate::graph::bubbles::index::BubbleIndex;
use crate::graph::traits::GraphNodeId;

/// Gap-affine-shaped heuristic, specialised to [`Affine`] (kept as the
/// concrete name most callers reach for).
pub type MinimumGapCostAffine<N, O> = MinimumGapCost<N, O, Affine>;

pub trait AstarHeuristic<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize;
}

/// Heuristic that always returns 0, reducing A* to plain Dijkstra: exact and
/// admissible, but explores every state within the optimal-score radius.
pub struct Dijkstra<N, O>(PhantomData<(N, O)>);

// Manual impl: `#[derive(Default)]` would incorrectly require `N: Default, O:
// Default` even though `PhantomData` doesn't need it.
impl<N, O> Default for Dijkstra<N, O> {
    fn default() -> Self {
        Self(PhantomData)
    }
}

impl<N, O> AstarHeuristic<N, O> for Dijkstra<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn h(&self, _aln_node: &AlignmentGraphNode<N, O>, _aln_state: AlignState) -> usize {
        0
    }
}

/// Admissible lower bound on the remaining alignment cost: the minimum gap
/// cost required to cover the (min, max) distance-to-end mismatch between
/// the remaining graph path and the remaining query length, using the
/// bubble index's cached shortest/longest path lengths to the graph end.
///
/// Generic over any cost model `C`: `C::gap_cost` already knows how to
/// compute an admissible per-state gap cost for its own shape (single-piece
/// affine/linear, or the cheaper-of-two-pieces minimum for two-piece), so
/// this heuristic doesn't need a model-specific variant.
///
/// Mode-aware: under `AlignmentMode::EndsFree`, a gap that falls entirely
/// within the free leading/trailing region costs nothing, so `h` returns 0
/// for that dimension rather than charging `C::gap_cost` — otherwise the
/// heuristic would overestimate (charge for a gap the alignment mode doesn't
/// actually penalise) and stop being admissible.
pub struct MinimumGapCost<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    costs: C,
    bubble_index: Arc<BubbleIndex<N>>,
    seq_length: usize,
    mode: AlignmentMode,
    _marker: PhantomData<O>,
}

impl<N, O, C> MinimumGapCost<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    pub fn new(costs: C, bubble_index: Arc<BubbleIndex<N>>, seq_length: usize) -> Self {
        Self::new_with_mode(costs, bubble_index, seq_length, AlignmentMode::Global)
    }

    pub fn new_with_mode(costs: C, bubble_index: Arc<BubbleIndex<N>>, seq_length: usize, mode: AlignmentMode) -> Self {
        Self {
            costs,
            bubble_index,
            seq_length,
            mode,
            _marker: PhantomData,
        }
    }
}

impl<N, O, C> AstarHeuristic<N, O> for MinimumGapCost<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, mut aln_state: AlignState) -> usize {
        let min_dist_to_exit = self.bubble_index.get_min_dist_to_end(aln_node.node()).saturating_sub(1);
        let max_dist_to_exit = self.bubble_index.get_max_dist_to_end(aln_node.node()).saturating_sub(1);

        let target_offset_min = aln_node.offset().as_usize() + min_dist_to_exit;
        let target_offset_max = aln_node.offset().as_usize() + max_dist_to_exit;

        let (qry_free_end, graph_free_end) = match self.mode {
            AlignmentMode::Global => (std::ops::Bound::Excluded(0), std::ops::Bound::Excluded(0)),
            AlignmentMode::EndsFree { qry_free_end, graph_free_end, .. } => (qry_free_end, graph_free_end),
        };

        let min_gap_length;
        if target_offset_min > self.seq_length {
            // Requires deletions (graph runs further before the query catches up). Free if
            // the remaining graph path to the true end fits within the free-end budget.
            let graph_free_hops = resolve_bound(graph_free_end, usize::MAX);
            if min_dist_to_exit <= graph_free_hops {
                return 0;
            }

            min_gap_length = target_offset_min - self.seq_length;

            // Already being in *either* piece's deletion state means no re-open cost.
            if !matches!(aln_state, AlignState::Deletion | AlignState::Deletion2) {
                aln_state = AlignState::Match;
            }
        } else if target_offset_max < self.seq_length {
            // Requires insertions (query has more left after the graph path ends). Free if
            // the remaining query fits within the free-end budget.
            let remaining_query = self.seq_length - target_offset_max;
            let max_qry_free = resolve_bound(qry_free_end, self.seq_length);
            if remaining_query <= max_qry_free {
                return 0;
            }

            min_gap_length = remaining_query;

            // Already being in *either* piece's insertion state means no re-open cost.
            if !matches!(aln_state, AlignState::Insertion | AlignState::Insertion2) {
                aln_state = AlignState::Match;
            }
        } else {
            min_gap_length = 0;
        }

        self.costs.gap_cost(aln_state, min_gap_length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::mock::create_linear_graph;
    use crate::graph::traits::GraphBase;

    type O = u32;

    #[test]
    fn test_min_gap_cost() {
        // mismatch=4, gap_open=6, gap_extend=2 (matches stable's
        // `GapAffine::new(4, 2, 6)`, whose constructor order is
        // mismatch/extend/open rather than ours).
        let costs = Affine::new(4, 6, 2);

        // Chain of 7 nodes; the 2nd node (label 2) is exactly 5 edges from
        // the auto-detected end node, so its (min, max) dist-to-end is (5, 5) —
        // matching stable's hand-supplied `vec![(5, 5)]` bubble-index stub.
        let g = create_linear_graph(7);
        let node_under_test = g.all_nodes_iter().nth(1).unwrap();
        let bubble_index = Arc::new(BubbleIndex::new(&g));

        let heuristic: MinimumGapCostAffine<_, O> = MinimumGapCostAffine::new(costs, bubble_index, 10);

        let node1 = AlignmentGraphNode::new(node_under_test, 2u32);
        assert_eq!(heuristic.h(&node1, AlignState::Match), 14);
        assert_eq!(heuristic.h(&node1, AlignState::Deletion), 14);
        // If already in insertion state, we wouldn't need to incur the gap-open cost.
        assert_eq!(heuristic.h(&node1, AlignState::Insertion), 8);

        let node2 = AlignmentGraphNode::new(node_under_test, 7u32);
        assert_eq!(heuristic.h(&node2, AlignState::Match), 8);
        // If already in deletion state, we wouldn't need to incur the gap-open cost.
        assert_eq!(heuristic.h(&node2, AlignState::Deletion), 2);
        assert_eq!(heuristic.h(&node2, AlignState::Insertion), 8);

        let node3 = AlignmentGraphNode::new(node_under_test, 6u32);
        assert_eq!(heuristic.h(&node3, AlignState::Match), 0);
        assert_eq!(heuristic.h(&node3, AlignState::Deletion), 0);
        assert_eq!(heuristic.h(&node3, AlignState::Insertion), 0);
    }
}
