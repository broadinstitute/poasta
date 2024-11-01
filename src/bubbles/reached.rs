use std::cmp::{max, min};
use std::collections::BTreeSet;
use std::marker::PhantomData;
use std::ops::{RangeFrom, RangeTo};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, GetAlignmentCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::astar::AstarVisited;
use crate::bubbles::index::{BubbleIndex, NodeBubbleMap};
use crate::graphs::NodeIndexType;


pub struct ReachedBubbleExitsMatch<'a, V, N, O>
    where V: AstarVisited<N, O> + GetAlignmentCosts,
          N: NodeIndexType,
          O: OffsetType,
{
    visited: &'a V,
    reached_offsets: &'a BTreeSet<O>,
    seq_len: usize,
    dummy: PhantomData<N>,
}

impl<'a, V, N, O> ReachedBubbleExitsMatch<'a, V, N, O>
    where V: AstarVisited<N, O> + GetAlignmentCosts,
          N: NodeIndexType,
          O: OffsetType,
{
    pub fn new(visited: &'a V, reached_offsets: &'a BTreeSet<O>, seq_len: usize) -> Self {
        Self {
            visited,
            reached_offsets,
            seq_len,
            dummy: PhantomData,
        }
    }

    pub fn can_improve_bubble(
        &self, bubble_index: &BubbleIndex<N>,
        bubble: &NodeBubbleMap<N>,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
        current_score: &Score
    ) -> bool {
        if self.reached_offsets.is_empty() {
            return true;
        }

        if aln_node.node() == bubble.bubble_exit {
            return true;
        }

        // j^min and j^max
        let target_offset_min = aln_node.offset() + O::new(bubble.min_dist_to_exit);
        let target_offset_max = aln_node.offset() + O::new(bubble.max_dist_to_exit);

        let min_dist_to_end = bubble_index.get_min_dist_to_end(bubble.bubble_exit)
            .saturating_sub(1);

        if target_offset_max.as_usize() > self.seq_len {
            return true;
        }

        // eprintln!("Checking for improving over bubble, [o, o'^min, o'^max] = [{:?}, {target_offset_min:?}, {target_offset_max:?}]", aln_node.offset());
        // eprintln!("Bubble exit: {:?}", bubble.bubble_exit);

        let mut prev_reached = self.reached_offsets
            .range(RangeTo { end: target_offset_min }).next_back();

        let mut last_offset = None;

        for next_reached in self.reached_offsets
            .range(target_offset_min..=target_offset_max)
        {
            let offset1 = prev_reached.map_or(
                target_offset_min,
                |v| max(target_offset_min, *v + O::one())
            );

            // Gap-affine and 2-piece alignment costs only.
            // If already in a deletion or insertion state, we need to take into account that
            // the current state under test would not have to incur the gap open cost again,
            // while an implicitly opened gap from a reached Match state would.
            if aln_state == AlignState::Deletion {
                let next_reached_cost = self.visited.get_score(&AlignmentGraphNode::new(
                    bubble.bubble_exit,
                    *next_reached
                ), AlignState::Match);

                if next_reached_cost + self.visited.get_costs().gap_open() > *current_score {
                    return true
                }
            } else if aln_state == AlignState::Deletion2 {
                let next_reached_cost = self.visited.get_score(&AlignmentGraphNode::new(
                    bubble.bubble_exit,
                    *next_reached
                ), AlignState::Match);

                if next_reached_cost + self.visited.get_costs().gap_open2() > *current_score {
                    return true
                }
            }

            if let Some(prev) = prev_reached {
                if aln_state == AlignState::Insertion {
                    let prev_reached_cost = self.visited.get_score(&AlignmentGraphNode::new(
                        bubble.bubble_exit,
                        *prev
                    ), AlignState::Match);

                    if prev_reached_cost + self.visited.get_costs().gap_open() > *current_score {
                        return true
                    }
                } else if aln_state == AlignState::Insertion2 {
                    let prev_reached_cost = self.visited.get_score(&AlignmentGraphNode::new(
                        bubble.bubble_exit,
                        *prev
                    ), AlignState::Match);

                    if prev_reached_cost + self.visited.get_costs().gap_open2() > *current_score {
                        return true
                    }
                }
            }

            // eprintln!("Checking o1... {offset1:?} (at score {current_score:?}");
            if self
                .can_improve_at_offset(bubble.bubble_exit, offset1, current_score, prev_reached, Some(next_reached), min_dist_to_end)
            {
                return true;
            }

            let offset2 = min(target_offset_max, max(target_offset_min, *next_reached - O::one()));
            if offset2 != offset1 {
                // eprintln!("Checking o2... {offset2:?} (at score {current_score:?}");

                if self
                    .can_improve_at_offset(bubble.bubble_exit, offset2, current_score, prev_reached, Some(next_reached), min_dist_to_end)
                {
                    return true;
                }
            }

            prev_reached = Some(next_reached);
            last_offset = Some(offset2);
        }

        let next_reached = self.reached_offsets
            .range(RangeFrom { start: target_offset_max.saturating_add(&O::one()) })
            .next();

        if last_offset.is_none()
            && self.can_improve_at_offset(bubble.bubble_exit, target_offset_min, current_score, prev_reached, next_reached, min_dist_to_end)
        {
            return true;
        }


        if last_offset.map_or(true, |v| v < target_offset_max)
            && self.can_improve_at_offset(bubble.bubble_exit, target_offset_max, current_score, prev_reached, next_reached, min_dist_to_end)
        {
            return true;
        }

        if let Some(prev) = prev_reached {
            if aln_state == AlignState::Insertion {
                let prev_reached_cost = self.visited.get_score(&AlignmentGraphNode::new(
                    bubble.bubble_exit,
                    *prev
                ), AlignState::Match);

                if prev_reached_cost + self.visited.get_costs().gap_open() > *current_score {
                    return true
                }
            } else if aln_state == AlignState::Insertion2 {
                let prev_reached_cost = self.visited.get_score(&AlignmentGraphNode::new(
                    bubble.bubble_exit,
                    *prev
                ), AlignState::Match);

                if prev_reached_cost + self.visited.get_costs().gap_open2() > *current_score {
                    return true
                }
            }
        }

        // eprintln!("PRUNE");
        false
    }

    fn can_improve_at_offset(
        &self,
        bubble_node: N,
        offset_to_check: O,
        score: &Score,
        reached_bubble_left: Option<&O>,
        reached_bubble_right: Option<&O>,
        min_dist_to_end: usize,
    ) -> bool {
        let implicit_score = match (reached_bubble_left, reached_bubble_right) {
            (None, None) => None,
            (Some(offset_left), Some(offset_right)) => {
                let left_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, *offset_left), AlignState::Match);
                let right_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, *offset_right), AlignState::Match);

                let gap_from_left = offset_to_check - *offset_left;
                let gap_from_right = *offset_right - offset_to_check;

                let implicit_score_from_left = left_score
                    + self.visited.get_costs().gap_cost(AlignState::Match, gap_from_left.as_usize());

                let implicit_score_from_right = right_score
                    + self.visited.get_costs().gap_cost(AlignState::Match, gap_from_right.as_usize());

                // eprintln!("Bubbles visited ({offset_left:?}, {offset_right:?}), scores: ({implicit_score_from_left:?}, {implicit_score_from_right:?})");

                // Only take deletion score into consideration if the deletion length is smaller than
                // the path length to the POA graph end node
                if gap_from_right.as_usize() > min_dist_to_end {
                    Some(implicit_score_from_left)
                } else {
                    Some(min(implicit_score_from_left, implicit_score_from_right))
                }
            },
            (None, Some(offset_right)) => {
                let right_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, *offset_right), AlignState::Match);
                let gap_from_right = *offset_right - offset_to_check;
                let implicit_score_from_right = right_score
                    + self.visited.get_costs().gap_cost(AlignState::Match, gap_from_right.as_usize());

                // eprintln!("Bubbles visited (None, {offset_right:?}), scores: (None, {implicit_score_from_right:?})");

                // Only take deletion score into consideration if the deletion length is smaller than
                // the path length to the POA graph end node
                if gap_from_right.as_usize() > min_dist_to_end {
                    None
                } else {
                    Some(implicit_score_from_right)
                }
            },
            (Some(offset_left), None) => {
                let left_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, *offset_left), AlignState::Match);
                let gap_from_left = offset_to_check - *offset_left;
                let implicit_score_from_left = left_score
                    + self.visited.get_costs().gap_cost(AlignState::Match, gap_from_left.as_usize());

                // eprintln!("Bubbles visited ({offset_left:?}, None), scores: ({implicit_score_from_left:?}, None)");

                Some(implicit_score_from_left)
            }
        };

        // eprintln!("Implicit gap score: {implicit_score:?}, score from align state: {score:?}");
        implicit_score.map_or(true, |s| *score < s)
    }
}
