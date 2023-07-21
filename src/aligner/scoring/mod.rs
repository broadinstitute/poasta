pub mod gap_linear;
pub mod gap_affine;

use std::fmt::Display;
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, StateTreeNode, TreeIndexType, Backtrace};
use crate::graphs::{AlignableGraph, NodeIndexType};

pub use gap_affine::GapAffine;
use crate::aligner::queue::AlignStateQueue;


pub trait AlignmentCosts: Copy {
    type StateTreeType<N, O, Ix>: AlignmentStateTree<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType;

    fn to_new_state_tree<G, N, O, Ix>(&self, graph: &G) -> Self::StateTreeType<N, O, Ix>
    where
        G: AlignableGraph<NodeIndex=N>,
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType;
}

pub trait AlignmentCostsEdit: AlignmentCosts {
    fn mismatch(&self) -> u8;
}

pub trait AlignmentCostsLinear: AlignmentCostsEdit {
    fn gap_extend(&self) -> u8;
}

pub trait AlignmentCostsAffine: AlignmentCostsLinear {
    fn gap_open(&self) -> u8;
}

pub trait AlignmentCostsTwoPiece: AlignmentCostsAffine {
    fn gap_extend2(&self) -> u8;
    fn gap_open2(&self) -> u8;
}

/// A trait that defines operations on the alignment state tree, and thus provides an interface
/// to generate new alignment states based on current states.
pub trait AlignmentStateTree<N, O, Ix>: Display
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    fn add_node(&mut self, node: StateTreeNode<N, O, Ix>) -> Ix;
    fn get_node(&self, node_ix: Ix) -> &StateTreeNode<N, O, Ix>;
    fn num_nodes(&self) -> usize;

    fn visited(&self, node: N, offset: O, state: AlignState) -> bool;
    fn mark_visited(&mut self, node: N, offset: O, state: AlignState);

    fn close_indels_for(&mut self, node_indices: &[Ix]) -> Vec<Ix>;

    fn add_extended_path(&mut self, start_ix: Ix, path: impl Iterator<Item=(N, O)>) -> Ix {
        let mut parent_tree_node = start_ix;
        for (node, offset) in path {
            let new_state = StateTreeNode::new(node, offset, AlignState::Match,
                Backtrace::Step(parent_tree_node));
            parent_tree_node = self.add_node(new_state)
        }

        parent_tree_node
    }

    fn generate_next<G>(
        &mut self,
        queue: &mut AlignStateQueue<Ix>,
        graph: &G,
        seq_len: usize,
        curr_ix: Ix,
    ) -> Option<Ix>
    where
        G: AlignableGraph<NodeIndex=N>;
}
