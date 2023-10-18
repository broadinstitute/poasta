pub mod gap_linear;
pub mod gap_affine;

use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, StateGraph};
use crate::graphs::{AlignableGraph, NodeIndexType};

pub use gap_affine::GapAffine;


pub trait AlignmentCosts: Copy {
    type StateGraphType<N, O>: StateGraph<N, O>
    where
        N: NodeIndexType,
        O: OffsetType;

    fn new_state_graph<G: AlignableGraph, O: OffsetType>(&self, graph: &G) -> Self::StateGraphType<G::NodeIndex, O>;

    fn mismatch(&self) -> u8;
    
    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;
    
    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;

    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize;
}
