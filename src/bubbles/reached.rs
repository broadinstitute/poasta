use std::collections::BTreeMap;
use std::ops::{RangeFrom, RangeTo};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::bubbles::index::{BubbleIndex, NodeBubbleMap};
use crate::graphs::{AlignableRefGraph, NodeIndexType};


pub struct ReachedBubbleExits<C, O>
where
    C: AlignmentCosts,
    O: OffsetType,
{
    costs: C,
    reached_exits: Vec<BTreeMap<O, Score>>,
}

impl<C, O> ReachedBubbleExits<C, O>
where
    C: AlignmentCosts,
    O: OffsetType,
{
    pub fn new<G: AlignableRefGraph>(costs: C, ref_graph: &G) -> Self {
        Self {
            costs,
            reached_exits: vec![BTreeMap::default();
                                ref_graph.node_count_with_start_and_end()],
        }
    }

    pub fn mark_reached<N: NodeIndexType>(&mut self, node: N, offset: O, score: Score) {
        self.reached_exits[node.index()]
            .insert(offset, score);
    }

    pub fn can_improve_alignment<N>(
        &self, bubble_index: &BubbleIndex<N>,
        aln_node: &AlignmentGraphNode<N, O>,
        current_score: Score
    ) -> bool
        where N: NodeIndexType,
    {
        if !bubble_index.node_is_part_of_bubble(aln_node.node()) {
            return true;
        }

        bubble_index.get_node_bubbles(aln_node.node())
            .iter()
            .all(|bubble|
                self.can_improve_bubble(bubble, aln_node, current_score))
    }

    pub fn can_improve_bubble<N>(
        &self, bubble: &NodeBubbleMap<N>,
        aln_node: &AlignmentGraphNode<N, O>,
        current_score: Score
    ) -> bool
        where N: NodeIndexType
    {
        if self.reached_exits[bubble.bubble_exit.index()].is_empty() {
            return true;
        }

        if aln_node.node() == bubble.bubble_exit {
            return true;
        }

        let target_offset = aln_node.offset() + O::new(bubble.dist_to_exit);

        // CASE 1: We could open an insertion from a reached bubble exit. Check if the most
        // optimal path from the current state (i.e., all matches, thus going diagonally and no
        // increase in alignment score) could improve over the score as reached when opening a
        // insertion from the bubble exit.
        if let Some((prev_offset, prev_score)) = self.reached_exits[bubble.bubble_exit.index()]
            .range(RangeTo { end: target_offset }).next_back()
        {
            let gap_length = target_offset - *prev_offset;
            let ins_score_from_exit = *prev_score
                + self.costs.gap_cost(AlignState::Match, gap_length.as_usize());

            if ins_score_from_exit <= current_score {
                // Previous offset that reached this exit can reach this new offset with a lower or equal score,
                return false;
            }
        }

        // CASE 2: We could open a deletion from a reached bubble exit. Check if the most
        // optimal path from the current state (i.e., all matches, thus going diagonally and no
        // increase in alignment score) could improve over the score as reached when opening a
        // deletion from the bubble exit.
        if let Some((next_offset, next_score)) = self.reached_exits[bubble.bubble_exit.index()]
            .range(RangeFrom { start: target_offset })
            .next()
        {
            let gap_length = *next_offset - target_offset;
            let del_score_from_exit = *next_score
                + self.costs.gap_cost(AlignState::Match, gap_length.as_usize());

            if del_score_from_exit <= current_score {
                return false;
            }
        }

        true
    }
}


#[cfg(test)]
mod tests {
    use petgraph::graph::NodeIndex;
    use crate::aligner::aln_graph::AlignmentGraphNode;
    use crate::aligner::scoring::{GapAffine, Score};
    use crate::bubbles::index::NodeBubbleMap;
    use crate::bubbles::reached::ReachedBubbleExits;
    use crate::graphs::mock::create_test_graph1;

    #[test]
    fn test_is_able_to_improve_bubble() {
        let g = create_test_graph1();

        let mut reached = ReachedBubbleExits::new(
            GapAffine::new(4, 2, 6),
            &g
        );

        reached.mark_reached(NodeIndex::<u32>::new(3), 5u32, Score::Score(4));

        let bubble_map = NodeBubbleMap::new(NodeIndex::<u32>::new(3), 3);

        // Score of 4 improves over opening a deletion from bubble exit
        assert!(reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 0u32), Score::Score(4),
        ));

        // Score of 14 and 16 do not improve over opening a deletion from bubble exit
        assert!(!reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 0u32), Score::Score(14),
        ));
        assert!(!reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 0u32), Score::Score(16),
        ));

        // Would reach bubble exit at same score and offset
        assert!(!reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 2u32), Score::Score(4),
        ));

        // Improves over opening an insertion from bubble exit
        assert!(reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 4u32), Score::Score(4),
        ));

        // Score of 14 and 16 do not improve over opening an insertion from bubble exit
        assert!(!reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 4u32), Score::Score(14),
        ));
        assert!(!reached.can_improve_bubble(
            &bubble_map, &AlignmentGraphNode::new(NodeIndex::<u32>::new(0), 4u32), Score::Score(14),
        ));

    }
}
