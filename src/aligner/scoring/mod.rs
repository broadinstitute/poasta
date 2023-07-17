pub mod gap_linear;
pub mod gap_affine;

use std::fmt::Display;
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, StateTreeNode, TreeIndexType, Backtrace};
use crate::graphs::{AlignableGraph, NodeIndexType};

pub use gap_affine::GapAffine;


pub trait AlignmentCosts: Copy {
    type StateTreeType<N, O, Ix>: AlignmentStateTree<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType;

    fn to_new_state_tree<N, O, Ix>(&self) -> Self::StateTreeType<N, O, Ix>
    where
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

    fn add_extended_path(&mut self, start_ix: Ix, path: Vec<(N, O)>) -> Ix {
        let last = path.last().unwrap();

        // Store matching nodes other than the last with the edge to the parent state
        let edge_nodes: Vec<N> = path[..path.len()-1].iter()
            .map(|v| v.0)
            .collect();

        let new_state = StateTreeNode::new(last.0, last.1, AlignState::Match,
                                           Backtrace::ExtraMatches(start_ix, edge_nodes));
        self.add_node(new_state)
    }

    fn generate_next<G>(
        &mut self,
        graph: &G,
        seq_len: usize,
        curr_ix: Ix,
    ) -> Vec<(u8, Ix)>
    where
        G: AlignableGraph<NodeIndex=N>;
}
