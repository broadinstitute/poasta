use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::{RangeFrom, RangeTo};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::bubbles::index::{BubbleIndex, NodeBubbleMap};
use crate::graphs::{AlignableRefGraph, NodeIndexType};


pub struct ReachedBubbleExits<C, O, S>
    where C: AlignmentCosts,
          O: OffsetType,
          S: CanImproveAlnCheck,
{
    costs: C,
    reached_exits: Vec<BTreeMap<O, Score>>,
    seq_len: usize,
    dummy: PhantomData<S>,
}

impl<C, O, S> ReachedBubbleExits<C, O, S>
    where C: AlignmentCosts,
          O: OffsetType,
          S: CanImproveAlnCheck,
{
    pub fn new<G: AlignableRefGraph>(costs: C, ref_graph: &G, seq_len: usize) -> Self {
        Self {
            costs,
            reached_exits: vec![BTreeMap::default();
                                ref_graph.node_count_with_start_and_end()],
            seq_len,
            dummy: PhantomData,
        }
    }

    pub fn mark_reached<N: NodeIndexType>(&mut self, node: N, offset: O, score: Score) {
        self.reached_exits[node.index()]
            .insert(offset, score);
    }

    pub fn can_improve_alignment<N>(
        &self, bubble_index: &BubbleIndex<N>,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
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
                self.can_improve_bubble(bubble_index, bubble, aln_node, aln_state, current_score))
    }

    pub fn can_improve_bubble<N>(
        &self, bubble_index: &BubbleIndex<N>,
        bubble: &NodeBubbleMap<N>,
        aln_node: &AlignmentGraphNode<N, O>,
        _: AlignState,
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

        let target_offset_min = aln_node.offset() + O::new(bubble.min_dist_to_exit);
        let target_offset_max = aln_node.offset() + O::new(bubble.max_dist_to_exit);
        let min_dist_to_end = bubble_index.get_min_dist_to_end(bubble.bubble_exit)
            .saturating_sub(1);

        // eprintln!("Checking for improving over bubble, [o^min, o^max] = [{target_offset_min:?}, {target_offset_max:?}]");

        let mut prev_reached = self.reached_exits[bubble.bubble_exit.index()]
            .range(RangeTo { end: target_offset_min }).next_back();

        // eprintln!("- prev reached before loop: {prev_reached:?}");

        for next_reached in self.reached_exits[bubble.bubble_exit.index()]
            .range(target_offset_min..=target_offset_max)
        {
            if S::can_improve_bubble(&self.costs, self.seq_len, min_dist_to_end, &current_score,
                target_offset_min, target_offset_max, prev_reached, Some(next_reached))
            {
                return true;
            }

            prev_reached = Some(next_reached);
        }

        // Let's check if a bubble has been reached at an offset > target_offset_max
        let next_reached = self.reached_exits[bubble.bubble_exit.index()]
            .range(RangeFrom { start: target_offset_max.saturating_add(&O::one()) }).next();

        if S::can_improve_bubble(&self.costs, self.seq_len, min_dist_to_end, &current_score,
            target_offset_min, target_offset_max, prev_reached, next_reached)
        {
            return true;
        }

        false
    }
}

pub trait CanImproveAlnCheck {
    fn can_improve_bubble<C, O>(
        costs: &C,
        seq_len: usize,
        min_dist_to_end: usize,
        current_score: &Score,
        target_offset_min: O,
        target_offset_max: O,
        prev_reached_bubble: Option<(&O, &Score)>,
        next_reached_bubble: Option<(&O, &Score)>
    ) -> bool
        where C: AlignmentCosts,
              O: OffsetType;
}

pub struct ReachedMatch;

impl CanImproveAlnCheck for ReachedMatch {
    fn can_improve_bubble<C, O>(
        costs: &C,
        seq_len: usize,
        min_dist_to_end: usize,
        current_score: &Score,
        target_offset_min: O,
        target_offset_max: O,
        prev_reached_bubble: Option<(&O, &Score)>,
        next_reached_bubble: Option<(&O, &Score)>
    ) -> bool
        where C: AlignmentCosts,
              O: OffsetType
    {
        match (prev_reached_bubble, next_reached_bubble) {
            (None, None) => true,
            (Some((prev_offset, prev_score)), Some((next_offset, next_score))) => {
                // eprintln!("- checking bubbles ({prev_offset:?}, score: {prev_score:?}), ({next_offset:?}, score: {next_score:?})");
                let gap_length = *next_offset - *prev_offset - O::one();
                let ins_score_from_prev = *prev_score
                    + costs.gap_cost(AlignState::Match, gap_length.as_usize());

                let can_improve_ins = *current_score < ins_score_from_prev;

                let del_score_from_next = *next_score
                    + costs.gap_cost(AlignState::Match, gap_length.as_usize());

                let can_improve_del = *current_score < del_score_from_next
                    && gap_length.as_usize() <= min_dist_to_end;

                // eprintln!("- ins: {can_improve_ins}, del: {can_improve_del}");
                can_improve_ins && can_improve_del
            },
            (None, Some((next_offset, next_score))) => {
                // eprintln!("- checking bubbles (None), ({next_offset:?}, score: {next_score:?})");
                let gap_length = *next_offset - target_offset_min;
                let del_score_from_next = *next_score +
                    costs.gap_cost(AlignState::Match, gap_length.as_usize());

                if gap_length.as_usize() > min_dist_to_end {
                    true
                } else {
                    *current_score < del_score_from_next
                }
            },
            (Some((prev_offset, prev_score)), None) => {
                // eprintln!("- checking bubbles ({prev_offset:?}, score: {prev_score:?}), (None)");
                if target_offset_max.as_usize() > seq_len {
                    true
                } else {
                    let gap_length = target_offset_max - *prev_offset;
                    let ins_score_from_prev = *prev_score
                        + costs.gap_cost(AlignState::Match, gap_length.as_usize());

                    // eprintln!("- ins score: {ins_score_from_prev:?}, curr: {current_score:?}");

                    *current_score < ins_score_from_prev
                }
            },
        }
    }
}


pub struct ReachedInsertion;

impl CanImproveAlnCheck for ReachedInsertion {
    fn can_improve_bubble<C, O>(
        costs: &C,
        seq_len: usize,
        _: usize,
        current_score: &Score,
        _: O,
        target_offset_max: O,
        prev_reached_bubble: Option<(&O, &Score)>,
        next_reached_bubble: Option<(&O, &Score)>
    ) -> bool
        where C: AlignmentCosts,
              O: OffsetType
    {
        match (prev_reached_bubble, next_reached_bubble) {
            (None, None) => true,
            (Some((prev_offset, prev_score)), Some((next_offset, _))) => {
                // eprintln!("- checking bubbles ({prev_offset:?}, score: {prev_score:?}), ({next_offset:?})");
                let gap_length = *next_offset - *prev_offset - O::one();
                let ins_score_from_prev = *prev_score
                    + costs.gap_cost(AlignState::Insertion, gap_length.as_usize());

                // eprintln!("- ins score from prev: {ins_score_from_prev:?}, current: {current_score:?}");
                *prev_score < *current_score && *current_score < ins_score_from_prev
            },
            (None, Some(_)) => true,
            (Some((prev_offset, prev_score)), None) => {
                // eprintln!("- checking bubbles ({prev_offset:?}, score: {prev_score:?}), (None)");
                if target_offset_max.as_usize() > seq_len {
                    true
                } else {
                    let gap_length = target_offset_max - *prev_offset;
                    let ins_score_from_prev = *prev_score
                        + costs.gap_cost(AlignState::Insertion, gap_length.as_usize());
                    // eprintln!("- ins score from prev: {ins_score_from_prev:?}, current: {current_score:?}");

                    *prev_score < *current_score && *current_score < ins_score_from_prev
                }
            },
        }
    }
}


pub struct ReachedDeletion;

impl CanImproveAlnCheck for ReachedDeletion {
    fn can_improve_bubble<C, O>(
        costs: &C,
        _: usize,
        min_dist_to_end: usize,
        current_score: &Score,
        target_offset_min: O,
        _: O,
        prev_reached_bubble: Option<(&O, &Score)>,
        next_reached_bubble: Option<(&O, &Score)>
    ) -> bool
        where C: AlignmentCosts,
              O: OffsetType
    {
        match (prev_reached_bubble, next_reached_bubble) {
            (None, None) => true,
            (Some((prev_offset, _)), Some((next_offset, next_score))) => {
                // eprintln!("- checking bubbles ({prev_offset:?}), ({next_offset:?}, score: {next_score:?})");
                let gap_length = *next_offset - *prev_offset - O::one();
                let del_score_from_next = *next_score
                    + costs.gap_cost(AlignState::Deletion, gap_length.as_usize());

                // eprintln!("- del score from next: {del_score_from_next:?}, current: {current_score:?}");
                *current_score < del_score_from_next
                    && gap_length.as_usize() <= min_dist_to_end
            },
            (None, Some((next_offset, next_score))) => {
                // eprintln!("- checking bubbles (None), ({next_offset:?}, score: {next_score:?})");
                let gap_length = *next_offset - target_offset_min;
                let del_score_from_next = *next_score +
                    costs.gap_cost(AlignState::Deletion, gap_length.as_usize());

                // eprintln!("- del score from next: {del_score_from_next:?}, current: {current_score:?}");
                if gap_length.as_usize() > min_dist_to_end {
                    true
                } else {
                    *current_score < del_score_from_next
                }
            },
            (Some(_), None) => {
                true
            },
        }
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
