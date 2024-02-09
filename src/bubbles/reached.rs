use std::cmp::{max, min};
use std::marker::PhantomData;
use std::ops::{RangeFrom, RangeTo};
use range_set_blaze::{RangeSetBlaze, Integer};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, GetAlignmentCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::astar::AstarVisited;
use crate::bubbles::index::{BubbleIndex, NodeBubbleMap};
use crate::graphs::NodeIndexType;


pub struct ReachedBubbleExitsMatch<'a, V, N, O>
    where V: AstarVisited<N, O> + GetAlignmentCosts,
          N: NodeIndexType,
          O: OffsetType + Integer,
{
    visited: &'a V,
    reached_offsets: &'a RangeSetBlaze<O>,
    seq_len: usize,
    dummy: PhantomData<N>,
}

impl<'a, V, N, O> ReachedBubbleExitsMatch<'a, V, N, O>
    where V: AstarVisited<N, O> + GetAlignmentCosts,
          N: NodeIndexType,
          O: OffsetType + Integer,
{
    pub fn new(visited: &'a V, reached_offsets: &'a RangeSetBlaze<O>, seq_len: usize) -> Self {
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

        let target_offset_min = aln_node.offset() + O::new(bubble.min_dist_to_exit);
        let target_offset_max = aln_node.offset() + O::new(bubble.max_dist_to_exit);
        let min_dist_to_end = bubble_index.get_min_dist_to_end(bubble.bubble_exit)
            .saturating_sub(1);
        let max_dist_to_end = bubble_index.get_max_dist_to_end(bubble.bubble_exit)
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
                |v| max(target_offset_min, v + O::one())
            );

            // Gap-affine and 2-piece alignment costs only:
            // For the first reached bubble to the right of the current alignment state offset,
            // we additionally check if the current state is already in deletion state and if extending that deletion
            // can ultimately improve over a newly opened deletion from the bubble exit (that will incur the gap open cost).
            if aln_state == AlignState::Deletion || aln_state == AlignState::Deletion2
                && offset1 == target_offset_min
            {
                // Extend the gap as far as we can, such that we can still reach the query offset of the
                // next reached bubble with match edges.
                let num_required_match = next_reached - aln_node.offset();
                if num_required_match.as_usize() < max_dist_to_end {
                    let gap_ext = max_dist_to_end - (next_reached - aln_node.offset()).as_usize();
                    let longest_deletion = bubble_index.get_max_dist_to_end(bubble.bubble_exit)
                        .saturating_sub(1);

                    let gap_cost_ext = self.visited.get_costs().gap_cost(aln_state, gap_ext);
                    let gap_cost_open = self.visited.get_costs().gap_cost(AlignState::Match, longest_deletion);

                    // eprintln!("INDEL Extension test: gap ext len: {gap_ext:?} (cost: {gap_cost_ext:?})");
                    // eprintln!("INDEL Extension test: gap open len: {longest_deletion:?} (cost: {gap_cost_open:?})");

                    if gap_cost_ext < gap_cost_open {
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

            let offset2 = min(target_offset_max, max(target_offset_min, next_reached - O::one()));
            if offset2 != offset1 {
                // eprintln!("Checking o2... {offset2:?} (at score {current_score:?}");

                if self
                    .can_improve_at_offset(bubble.bubble_exit, offset2, current_score, prev_reached, Some(next_reached), min_dist_to_end)
                {
                    return true;
                }
            }

            if aln_state == AlignState::Insertion || aln_state == AlignState::Insertion2 {
                if let Some(offset_left) = prev_reached {
                    // For gap-affine and 2-piece alignment costs only:
                    // If the current alignment state is already in insertion state, check if this state
                    // could ultimately improve over an existing bubble by extending this insertion
                    // up until we reach `next_reached` with match edges
                    let ext_end_offset = next_reached.as_usize() - bubble.max_dist_to_exit;
                    if ext_end_offset > aln_node.offset().as_usize() {
                        let gap_ext = (ext_end_offset - aln_node.offset().as_usize()).saturating_sub(1);
                        let longest_insertion = (next_reached - offset_left).as_usize().saturating_sub(1);

                        let gap_cost_ext = self.visited.get_costs().gap_cost(aln_state, gap_ext);
                        let gap_cost_open = self.visited.get_costs().gap_cost(AlignState::Match, longest_insertion);

                        // eprintln!("INS Extension test: gap ext len: {gap_ext:?} (cost: {gap_cost_ext:?})");
                        // eprintln!("INS Extension test: gap open len: {longest_insertion:?} (cost: {gap_cost_open:?})");

                        if gap_cost_ext < gap_cost_open {
                            return true
                        }

                    }
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

        if let Some(offset_left) = prev_reached {
            if aln_state == AlignState::Insertion || aln_state == AlignState::Insertion2 {
                // For gap-affine and 2-piece alignment costs only:
                // If the current alignment state is already in insertion state, check if this state
                // could ultimately improve over an existing bubble by extending this insertion as far as we can,
                // such that we still can reach the bubble exit with match edges.
                let ext_end_offset = self.seq_len - bubble.max_dist_to_exit;
                if ext_end_offset > aln_node.offset().as_usize() {
                    let gap_ext = ext_end_offset - aln_node.offset().as_usize();
                    let longest_insertion = self.seq_len - offset_left.as_usize();

                    let gap_cost_ext = self.visited.get_costs().gap_cost(aln_state, gap_ext);
                    let gap_cost_open = self.visited.get_costs().gap_cost(AlignState::Match, longest_insertion);

                    // eprintln!("INS Extension extreme test: gap ext len: {gap_ext:?} (cost: {gap_cost_ext:?})");
                    // eprintln!("INS Extension extreme test: gap open len: {longest_insertion:?} (cost: {gap_cost_open:?})");

                    if gap_cost_ext < gap_cost_open {
                        return true
                    }
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
        reached_bubble_left: Option<O>,
        reached_bubble_right: Option<O>,
        min_dist_to_end: usize,
    ) -> bool {
        let implicit_score = match (reached_bubble_left, reached_bubble_right) {
            (None, None) => None,
            (Some(offset_left), Some(offset_right)) => {
                let left_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, offset_left), AlignState::Match);
                let right_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, offset_right), AlignState::Match);

                let gap_from_left = offset_to_check - offset_left;
                let gap_from_right = offset_right - offset_to_check;

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
                let right_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, offset_right), AlignState::Match);
                let gap_from_right = offset_right - offset_to_check;
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
                let left_score = self.visited.get_score(&AlignmentGraphNode::new(bubble_node, offset_left), AlignState::Match);
                let gap_from_left = offset_to_check - offset_left;
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
