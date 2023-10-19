use std::collections::BTreeMap;
use std::ops::{RangeFrom, RangeTo};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentCosts;
use crate::aligner::state::{AlignState, Score, StateGraphNode};
use crate::bubbles::index::{BubbleIndex, NodeBubbleMap};
use crate::graphs::{AlignableGraph, NodeIndexType};


pub struct ReachedBubbleExits<'a, O, C>
where
    O: OffsetType,
    C: AlignmentCosts,
{
    costs: &'a C,
    seq: &'a [u8],
    reached_exits: Vec<BTreeMap<O, Score>>
}

impl<'a, O, C> ReachedBubbleExits<'a, O, C>
where
    O: OffsetType,
    C: AlignmentCosts,
{
    pub fn new<G: AlignableGraph>(costs: &'a C, graph: &G, seq: &'a [u8]) -> Self {
        Self {
            costs,
            seq,
            reached_exits: vec![BTreeMap::default(); graph.node_count_with_start()]
        }
    }

    pub fn mark_reached<N: NodeIndexType>(&mut self, node: N, offset: O, score: Score) {
        self.reached_exits[node.index()]
            .insert(offset, score);
    }

    pub fn can_improve_alignment<S, N>(&self, bubble_index: &BubbleIndex<N, O>, state: &S, current_score: Score) -> bool
    where
        S: StateGraphNode<N, O>,
        N: NodeIndexType,
    {
        if !bubble_index.node_is_part_of_bubble(state.node()) {
            return true;
        }

        bubble_index.get_node_bubbles(state.node())
            .iter()
            .all(|bubble|
                self.can_improve_bubble(bubble, state, current_score))
    }

    pub fn can_improve_bubble<S, N>(&self, bubble: &NodeBubbleMap<N, O>, state: &S, current_score: Score) -> bool
    where
        S: StateGraphNode<N, O>,
        N: NodeIndexType
    {
        if self.reached_exits[bubble.bubble_exit.index()].is_empty() {
            return true;
        }

        let target_offset = state.offset() + bubble.dist_to_exit;

        if target_offset.as_usize() >= self.seq.len() {
            return false;
        }

        if state.node() == bubble.bubble_exit {
            return true;
        }

        // CASE 1: We could open an insertion from a reached bubble exit. Check if the most
        // optimal path from the current state (i.e., all matches, thus going diagonally and no
        // increase in alignment score) could improve over the score as reached when opening a
        // insertion from the bubble exit.
        if let Some((prev_offset, prev_score)) = self.reached_exits[bubble.bubble_exit.index()]
            .range(RangeTo { end: target_offset })
            .rev()
            .next()
        {
            let gap_length = target_offset - *prev_offset;
            let ins_open_score_from_exit = *prev_score
                + self.costs.gap_cost(AlignState::Match, gap_length.as_usize());

            if ins_open_score_from_exit <= current_score {
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
            let del_open_score_from_exit = *next_score
                + self.costs.gap_cost(AlignState::Match, gap_length.as_usize());

            if del_open_score_from_exit <= current_score {
                return false;
            }
        }

        true
    }

}
