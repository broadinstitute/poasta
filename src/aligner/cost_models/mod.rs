use super::{astar::AstarAlignableGraph, fr_points::IndexType, AlignmentMode};

pub mod affine;

pub trait AlignmentCostModel {
    type AstarState<N, D>;
    
    fn initialize<G, D>(&self, graph: &G, alignment_mode: AlignmentMode) -> Self::AstarState<G::NodeType, D>
    where
        G: AstarAlignableGraph,
        D: IndexType;
    
    fn mismatch(&self) -> u8;
    
    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;
    
    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;
}