use std::fmt::{Display, Formatter};

use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::io::state_tree::write_tree_gml;
use crate::aligner::offsets::OffsetType;
use crate::aligner::queue::AlignStateQueue;
use crate::aligner::scoring::{AlignmentCosts, AlignmentStateTree};
use crate::aligner::state::{AlignState, StateTree, StateTreeNode, Backtrace, TreeIndexType};
use crate::aligner::visited::{VisitedSet, VisitedSetPerNode};

#[derive(Clone, Copy, Debug)]
pub struct GapAffine {
    cost_mismatch: u8,
    cost_gap_extend: u8,
    cost_gap_open: u8
}

impl GapAffine {
    pub fn new(cost_mismatch: u8, cost_gap_extend: u8, cost_gap_open: u8) -> Self {
        Self { cost_mismatch, cost_gap_extend, cost_gap_open }
    }
}

impl AlignmentCosts for GapAffine {
    type StateTreeType<N, O, Ix> = GapAffineStateTree<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType;

    fn to_new_state_tree<G, N, O, Ix>(&self, graph: &G) -> Self::StateTreeType<N, O, Ix>
    where
        G: AlignableGraph<NodeIndex=N>,
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType
    {
        GapAffineStateTree::new(*self, graph)
    }

    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.cost_mismatch
    }

    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.cost_gap_open
    }

    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.cost_gap_extend
    }

    #[inline(always)]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        0
    }
}


pub struct GapAffineStateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    costs: GapAffine,
    tree: StateTree<N, O, Ix>,
    visited_m: VisitedSetPerNode<O>,
    visited_d: VisitedSetPerNode<O>,
    visited_i: VisitedSetPerNode<O>,
}

impl<N, O, Ix> GapAffineStateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    fn new<G: AlignableGraph<NodeIndex=N>>(costs: GapAffine, graph: &G) -> Self {
        Self {
            costs,
            tree: StateTree::new(),
            visited_m: VisitedSetPerNode::new(graph),
            visited_d: VisitedSetPerNode::new(graph),
            visited_i: VisitedSetPerNode::new(graph),
        }
    }

    #[inline]
    fn is_end_node<G>(&self, graph: &G, seq_len: usize, state: Ix) -> bool
    where
        G: AlignableGraph<NodeIndex=N>
    {
        let node = self.tree.get_node(state);
        node.offset().as_usize() == seq_len && graph.is_end(node.node())
    }
}

impl<N, O, Ix> AlignmentStateTree<N, O, Ix> for GapAffineStateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    #[inline(always)]
    fn add_node(&mut self, node: StateTreeNode<N, O, Ix>) -> Ix {
        self.tree.add_node(node)
    }

    #[inline(always)]
    fn get_node(&self, node_ix: Ix) -> &StateTreeNode<N, O, Ix> {
        self.tree.get_node(node_ix)
    }

    fn num_nodes(&self) -> usize {
        self.tree.num_nodes()
    }

    #[inline]
    fn visited(&self, node: N, offset: O, state: AlignState) -> bool {
        match state {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => self.visited_m.visited(node, offset),
            AlignState::Deletion => self.visited_d.visited(node, offset),
            AlignState::Insertion => self.visited_i.visited(node, offset),
            _ => panic!("Invalid alignment state for GapAffine!")
        }
    }

    fn mark_visited(&mut self, node: N, offset: O, state: AlignState) {
        match state {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => self.visited_m.mark_visited(node, offset),
            AlignState::Deletion => self.visited_d.mark_visited(node, offset),
            AlignState::Insertion => self.visited_i.mark_visited(node, offset),
            _ => panic!("Invalid alignment state for GapAffine!")
        };
    }

    fn close_indels_for(&mut self, node_indices: &[Ix]) -> Vec<Ix> {
        node_indices.iter().filter_map(|v| {
            let node = self.tree.get_node(*v);

            match node.state() {
                AlignState::Deletion | AlignState::Insertion => {
                    if !self.visited(node.node(), node.offset(), AlignState::Match) {
                        self.visited_m.mark_visited(node.node(), node.offset());

                        let new_state = StateTreeNode::new(
                            node.node(), node.offset(), AlignState::Match, Backtrace::ClosedIndel(*v));
                        let new_ix = self.tree.add_node(new_state);

                        Some(new_ix)
                    } else {
                        None
                    }
                },
                _ => None
            }
        }).collect()
    }

    fn generate_next<G>(
        &mut self,
        queue: &mut AlignStateQueue<Ix>,
        graph: &G,
        seq_len: usize,
        curr_score: usize,
        curr_ix: Ix
    ) -> Option<Ix>
    where
        G: AlignableGraph<NodeIndex=N>,
        Ix: TreeIndexType,
    {
        let seq_len_as_o = O::new(seq_len);
        assert!(curr_ix.index() < self.tree.num_nodes());

        match self.get_node(curr_ix).state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch  => {
                // For each successor we can either enter a mismatch state or open a deletion
                let new_score_mis = curr_score + self.costs.mismatch() as usize;
                for succ in graph.successors(self.get_node(curr_ix).node()) {
                    // Mismatch, only if there's still query sequence to match
                    if self.get_node(curr_ix).offset() < seq_len_as_o {
                        let new_offset = self.get_node(curr_ix).offset().increase_one();
                        if !self.visited(succ, new_offset, AlignState::Match) {
                            self.visited_m.mark_visited(succ, new_offset);

                            let new_state = StateTreeNode::new(
                                succ, new_offset,
                                AlignState::Mismatch, Backtrace::Step(curr_ix));
                            let new_ix = self.tree.add_node(new_state);

                            if self.is_end_node(graph, seq_len, new_ix) {
                                return Some(new_ix)
                            }

                            queue.queue_endpoint(self.costs.mismatch() - 1, new_ix);
                        }
                    }

                    // Open deletion
                    let offset = self.get_node(curr_ix).offset();
                    if !self.visited(succ, offset, AlignState::Deletion) {
                        self.visited_d.mark_visited(succ, offset);

                        let new_state = StateTreeNode::new(
                            succ, offset, AlignState::Deletion, Backtrace::Step(curr_ix));
                        let new_ix = self.tree.add_node(new_state);

                        queue.queue_endpoint(self.costs.gap_open() + self.costs.gap_extend() - 1, new_ix);
                    }
                }

                if self.get_node(curr_ix).offset() < seq_len_as_o {
                    // Open insertion
                    let curr_node = self.get_node(curr_ix).node();
                    let new_offset = self.get_node(curr_ix).offset().increase_one();
                    if !self.visited(curr_node, new_offset, AlignState::Insertion) {
                        self.visited_i.mark_visited(curr_node, new_offset);

                        let new_state = StateTreeNode::new(
                            curr_node, new_offset,
                            AlignState::Insertion, Backtrace::Step(curr_ix));
                        let new_ix = self.tree.add_node(new_state);

                        queue.queue_endpoint(self.costs.gap_open() + self.costs.gap_extend() - 1, new_ix);
                    }
                }
            },
            AlignState::Deletion => {
                // Extend deletion for each successor
                for succ in graph.successors(self.get_node(curr_ix).node()) {
                    let offset = self.get_node(curr_ix).offset();
                    if !self.visited(succ, offset, AlignState::Deletion) {
                        self.visited_d.mark_visited(succ, offset);

                        let new_state = StateTreeNode::new(
                            succ, offset,
                            AlignState::Deletion, Backtrace::Step(curr_ix));
                        let new_ix = self.tree.add_node(new_state);
                        queue.queue_endpoint(self.costs.gap_extend() - 1, new_ix);
                    }
                }
            },
            AlignState::Insertion => {
                if self.get_node(curr_ix).offset() < seq_len_as_o {
                    let curr_node = self.get_node(curr_ix).node();
                    let new_offset = self.get_node(curr_ix).offset().increase_one();
                    if !self.visited(curr_node, new_offset, AlignState::Insertion) {
                        self.visited_i.mark_visited(curr_node, new_offset);

                        let new_state = StateTreeNode::new(
                            curr_node, new_offset,
                            AlignState::Insertion, Backtrace::Step(curr_ix));
                        let new_ix = self.tree.add_node(new_state);
                        queue.queue_endpoint(self.costs.gap_extend() - 1, new_ix);
                    }
                }
            }
            _ => panic!("Invalid alignment state for gap affine!")
        };

        None
    }
}

impl<N, O, Ix> Display for GapAffineStateTree<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write_tree_gml(f, self)?;

        Ok(())
    }
}
