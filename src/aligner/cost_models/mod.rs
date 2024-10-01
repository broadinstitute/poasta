use super::{astar::{AlignableGraphRef, AstarState}, fr_points::{DiagType, OffsetType}, AlignmentMode};

pub mod affine;

pub trait AlignmentCostModel {
    type AstarStateType<G, D, O>: AstarState<G>
    where
        G: AlignableGraphRef,
        D: DiagType,
        O: OffsetType;
    
    fn initialize<G, D, O>(&self, graph: G, seq: &[u8], alignment_mode: AlignmentMode) -> Self::AstarStateType<G, D, O>
    where
        G: AlignableGraphRef,
        D: DiagType,
        O: OffsetType;
    
    fn initial_states<G, D, O>(&self, graph: G, seq: &[u8], alignment_mode: AlignmentMode) -> impl Iterator<
        Item=<Self::AstarStateType<G, D, O> as AstarState<G>>::AstarItem
    >
    where
        G: AlignableGraphRef,
        D: DiagType,
        O: OffsetType;
    
    fn mismatch(&self) -> u8;
    
    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;
    
    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;
}