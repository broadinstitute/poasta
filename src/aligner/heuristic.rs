use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, GapAffine};
use crate::graphs::NodeIndexType;

pub trait AstarHeuristic {
    fn h<N, O>(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType;
}

/// A* heuristic that always returns 0, such that
/// A* reduces to standard Dijkstra's algorithm.
#[derive(Default)]
pub struct Dijkstra;


impl AstarHeuristic for Dijkstra {
    fn h<N, O>(&self, _: &AlignmentGraphNode<N, O>, _: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType
    {
        0
    }
}


pub struct MinimumGapCostAffine {
    costs: GapAffine,
    dist_to_end: Vec<usize>,
    seq_length: usize,
}

impl MinimumGapCostAffine {
    pub fn new(costs: GapAffine, dist_to_end: Vec<usize>, seq_length: usize) -> Self {
        Self {
            costs,
            dist_to_end,
            seq_length,
        }
    }
}

impl AstarHeuristic for MinimumGapCostAffine {

    fn h<N, O>(&self, aln_node: &AlignmentGraphNode<N, O>, mut aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType,
    {
        let dist_to_exit =  self.dist_to_end[aln_node.node().index()].saturating_sub(1);
        if dist_to_exit == 0 && aln_node.offset().as_usize() == self.seq_length {
            return 0;
        }

        let target_offset = aln_node.offset().as_usize() + dist_to_exit;

        let gap_length;
        if target_offset < self.seq_length {
            gap_length = self.seq_length - target_offset;

            // requires insertions, so if in deletion or match state, we need to open a new gap
            if aln_state != AlignState::Insertion {
                aln_state = AlignState::Match;
            }
        } else {
            gap_length = target_offset - self.seq_length;

            // requires deletions, so if in insertion or match state, we need to open a new gap
            if aln_state != AlignState::Deletion {
                aln_state = AlignState::Match;
            }
        }

        self.costs.gap_cost(aln_state, gap_length)
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

        let heuristic = MinimumGapCostAffine::new(costs, vec![5], 10);

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
