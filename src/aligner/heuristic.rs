use std::marker::PhantomData;
use std::rc::Rc;
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, GapAffine};
use crate::bubbles::index::BubbleIndex;
use crate::graphs::NodeIndexType;

pub trait AstarHeuristic<N, O> {
    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType;
}

/// A* heuristic that always returns 0, such that
/// A* reduces to standard Dijkstra's algorithm.
#[derive(Default)]
pub struct Dijkstra<N, O>(PhantomData<(N, O)>);


impl<N, O> AstarHeuristic<N, O> for Dijkstra<N, O> {
    fn h(&self, _: &AlignmentGraphNode<N, O>, _: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType
    {
        0
    }
}


pub struct MinimumGapCostAffine<N, O> {
    costs: GapAffine,
    bubble_index: Rc<BubbleIndex<N>>,
    seq_length: usize,
    dummy: PhantomData<O>,
}

impl<N, O> MinimumGapCostAffine<N, O> {
    pub fn new(costs: GapAffine, bubble_index: Rc<BubbleIndex<N>>, seq_length: usize) -> Self {
        Self {
            costs,
            bubble_index,
            seq_length,
            dummy: PhantomData,
        }
    }
}

impl<N, O> AstarHeuristic<N, O> for MinimumGapCostAffine<N, O> {

    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, mut aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType,
    {
        let min_dist_to_exit = self.bubble_index.get_min_dist_to_end(aln_node.node())
            .saturating_sub(1);
        let max_dist_to_exit = self.bubble_index.get_max_dist_to_end(aln_node.node())
            .saturating_sub(1);

        let target_offset_min = aln_node.offset().as_usize() + min_dist_to_exit;
        let target_offset_max = aln_node.offset().as_usize() + max_dist_to_exit;

        let min_gap_length;
        if target_offset_min > self.seq_length {
            min_gap_length = target_offset_min - self.seq_length;

            // requires deletions, so if in insertion or match state, we need to open a new gap
            if aln_state != AlignState::Deletion {
                aln_state = AlignState::Match;
            }
        } else if target_offset_max < self.seq_length {
            min_gap_length = self.seq_length - target_offset_max;

            // requires insertions, so if in deletion or match state, we need to open a new gap
            if aln_state != AlignState::Insertion {
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
    use super::AstarHeuristic;
    use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
    use crate::aligner::heuristic::MinimumGapCostAffine;
    use crate::aligner::scoring::GapAffine;

    #[test]
    fn test_min_gap_cost() {
        let costs = GapAffine::new(4, 2, 6);

        let heuristic = MinimumGapCostAffine::new(costs, vec![(5, 5)], 10);

        let node1 = AlignmentGraphNode::new(0u32, 2u32);
        assert_eq!(heuristic.h(&node1, AlignState::Match), 14);
        assert_eq!(heuristic.h(&node1, AlignState::Deletion), 14);
        // If already in insertion state, we wouldn't need to incur the gap-open cost
        assert_eq!(heuristic.h(&node1, AlignState::Insertion), 8);

        let node2 = AlignmentGraphNode::new(0u32, 7u32);
        assert_eq!(heuristic.h(&node2, AlignState::Match), 8);
        // If already in deletion state, we wouldn't need to incur the gap-open cost
        assert_eq!(heuristic.h(&node2, AlignState::Deletion), 2);
        assert_eq!(heuristic.h(&node2, AlignState::Insertion), 8);

        let node3 = AlignmentGraphNode::new(0u32, 6u32);
        assert_eq!(heuristic.h(&node3, AlignState::Match), 0);
        // If already in deletion state, we wouldn't need to incur the gap-open cost
        assert_eq!(heuristic.h(&node3, AlignState::Deletion), 0);
        assert_eq!(heuristic.h(&node3, AlignState::Insertion), 0);
    }
}
