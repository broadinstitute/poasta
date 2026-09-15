//! Bubble-exit pruning check.
//!
//! Given the set of query offsets at which a bubble's exit node has already
//! been reached (in the `Match` state) at some score, decide whether a
//! candidate state could still improve on that, that is, whether A*
//! needs to keep exploring it, or whether it's provably dominated.
//!
//! The check reads the exit's reached bits. A masked word scan finds the
//! neighbouring reached offsets, and the `Match` score at each is one
//! block lookup through `get_match_score`.

use std::cmp::{max, min};

use crate::align::AlignState;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::astar::visited::ReachedBits;
use crate::align::offsets::OffsetType;
use crate::graph::bubbles::index::{BubbleIndex, NodeBubbleMap};
use crate::graph::traits::GraphNodeId;

pub struct ReachedBubbleExitsMatch<'a, N, O, C, F>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
    F: Fn(N, O) -> u32,
{
    get_match_score: F,
    costs: &'a C,
    /// Query offsets where this bubble's exit was reached in the `Match`
    /// state.
    reached: &'a ReachedBits,
    seq_len: usize,
    _marker: std::marker::PhantomData<(N, O)>,
}

impl<'a, N, O, C, F> ReachedBubbleExitsMatch<'a, N, O, C, F>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
    F: Fn(N, O) -> u32,
{
    pub fn new(get_match_score: F, costs: &'a C, reached: &'a ReachedBits, seq_len: usize) -> Self {
        Self {
            get_match_score,
            costs,
            reached,
            seq_len,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn can_improve_bubble(
        &self,
        bubble_index: &BubbleIndex<N>,
        bubble: &NodeBubbleMap<N>,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
        current_score: u32,
    ) -> bool {
        if self.reached.is_empty() {
            return true;
        }

        if aln_node.node() == bubble.bubble_exit {
            return true;
        }

        // j^min and j^max
        let target_offset_min = aln_node.offset() + O::new(bubble.min_dist_to_exit as usize);
        let target_offset_max = aln_node.offset() + O::new(bubble.max_dist_to_exit as usize);

        let min_dist_to_end = bubble_index
            .get_min_dist_to_end(bubble.bubble_exit)
            .saturating_sub(1);

        if target_offset_max.as_usize() > self.seq_len {
            return true;
        }

        // The window `[target_offset_min, target_offset_max]` used to be located
        // with two `partition_point`s over the sorted offsets. It is now a bit
        // scan in the exit's bitmap.
        let (window_min, window_max) = (target_offset_min.as_usize() as u32, target_offset_max.as_usize() as u32);
        let mut prev_reached: Option<O> = self.reached.last_reached_below(window_min).map(|o| O::new(o as usize));

        let mut last_offset = None;

        for next_reached in self.reached.reached_in(window_min, window_max).map(|o| O::new(o as usize)) {
            let offset1 =
                prev_reached.map_or(target_offset_min, |v| max(target_offset_min, v + O::one()));

            // Gap-affine and 2-piece alignment costs only.
            // If already in a deletion or insertion state, we need to take into account that
            // the current state under test would not have to incur the gap open cost again,
            // while an implicitly opened gap from a reached Match state would.
            if aln_state == AlignState::Deletion {
                let next_reached_cost = (self.get_match_score)(bubble.bubble_exit, next_reached);
                if next_reached_cost + self.costs.gap_open() as u32 > current_score {
                    return true;
                }
            } else if aln_state == AlignState::Deletion2 {
                let next_reached_cost = (self.get_match_score)(bubble.bubble_exit, next_reached);
                if next_reached_cost + self.costs.gap_open2() as u32 > current_score {
                    return true;
                }
            }

            if let Some(prev) = prev_reached {
                if aln_state == AlignState::Insertion {
                    let prev_reached_cost = (self.get_match_score)(bubble.bubble_exit, prev);
                    if prev_reached_cost + self.costs.gap_open() as u32 > current_score {
                        return true;
                    }
                } else if aln_state == AlignState::Insertion2 {
                    let prev_reached_cost = (self.get_match_score)(bubble.bubble_exit, prev);
                    if prev_reached_cost + self.costs.gap_open2() as u32 > current_score {
                        return true;
                    }
                }
            }

            if self.can_improve_at_offset(
                bubble.bubble_exit,
                offset1,
                current_score,
                prev_reached,
                Some(next_reached),
                min_dist_to_end,
            ) {
                return true;
            }

            let offset2 = min(
                target_offset_max,
                max(target_offset_min, next_reached - O::one()),
            );
            if offset2 != offset1
                && self.can_improve_at_offset(
                    bubble.bubble_exit,
                    offset2,
                    current_score,
                    prev_reached,
                    Some(next_reached),
                    min_dist_to_end,
                )
            {
                return true;
            }

            prev_reached = Some(next_reached);
            last_offset = Some(offset2);
        }

        let next_reached = self.reached.first_reached_above(window_max).map(|o| O::new(o as usize));

        if last_offset.is_none()
            && self.can_improve_at_offset(
                bubble.bubble_exit,
                target_offset_min,
                current_score,
                prev_reached,
                next_reached,
                min_dist_to_end,
            )
        {
            return true;
        }

        if last_offset.is_none_or(|v| v < target_offset_max)
            && self.can_improve_at_offset(
                bubble.bubble_exit,
                target_offset_max,
                current_score,
                prev_reached,
                next_reached,
                min_dist_to_end,
            )
        {
            return true;
        }

        if let Some(prev) = prev_reached {
            if aln_state == AlignState::Insertion {
                let prev_reached_cost = (self.get_match_score)(bubble.bubble_exit, prev);
                if prev_reached_cost + self.costs.gap_open() as u32 > current_score {
                    return true;
                }
            } else if aln_state == AlignState::Insertion2 {
                let prev_reached_cost = (self.get_match_score)(bubble.bubble_exit, prev);
                if prev_reached_cost + self.costs.gap_open2() as u32 > current_score {
                    return true;
                }
            }
        }

        false
    }

    fn can_improve_at_offset(
        &self,
        bubble_node: N,
        offset_to_check: O,
        score: u32,
        reached_bubble_left: Option<O>,
        reached_bubble_right: Option<O>,
        min_dist_to_end: usize,
    ) -> bool {
        let implicit_score = match (reached_bubble_left, reached_bubble_right) {
            (None, None) => None,
            (Some(offset_left), Some(offset_right)) => {
                let left_score = (self.get_match_score)(bubble_node, offset_left);
                let right_score = (self.get_match_score)(bubble_node, offset_right);

                let gap_from_left = offset_to_check - offset_left;
                let gap_from_right = offset_right - offset_to_check;

                let implicit_score_from_left = left_score
                    + self
                        .costs
                        .gap_cost(AlignState::Match, gap_from_left.as_usize())
                        as u32;
                let implicit_score_from_right = right_score
                    + self
                        .costs
                        .gap_cost(AlignState::Match, gap_from_right.as_usize())
                        as u32;

                // Only take deletion score into consideration if the deletion length is smaller
                // than the path length to the graph end node
                if gap_from_right.as_usize() > min_dist_to_end {
                    Some(implicit_score_from_left)
                } else {
                    Some(implicit_score_from_left.min(implicit_score_from_right))
                }
            }
            (None, Some(offset_right)) => {
                let right_score = (self.get_match_score)(bubble_node, offset_right);
                let gap_from_right = offset_right - offset_to_check;
                let implicit_score_from_right = right_score
                    + self
                        .costs
                        .gap_cost(AlignState::Match, gap_from_right.as_usize())
                        as u32;

                if gap_from_right.as_usize() > min_dist_to_end {
                    None
                } else {
                    Some(implicit_score_from_right)
                }
            }
            (Some(offset_left), None) => {
                let left_score = (self.get_match_score)(bubble_node, offset_left);
                let gap_from_left = offset_to_check - offset_left;
                let implicit_score_from_left = left_score
                    + self
                        .costs
                        .gap_cost(AlignState::Match, gap_from_left.as_usize())
                        as u32;

                Some(implicit_score_from_left)
            }
        };

        implicit_score.is_none_or(|s| score < s)
    }
}
