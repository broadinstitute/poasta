use super::{
    astar::{AlignableGraph, AstarState},
    fr_points::{DiagType, PosType},
    AlignmentMode,
};

pub mod affine;

pub trait AlignmentCostModel {
    type AstarStateType<G, D, O>: AstarState<G>
    where
        G: AlignableGraph,
        D: DiagType,
        O: PosType;

    fn initialize<G, D, O, F>(
        &self,
        graph: &G,
        seq: &[u8],
        alignment_mode: AlignmentMode,
        heuristic: F,
    ) -> Self::AstarStateType<G, D, O>
    where
        G: AlignableGraph,
        D: DiagType,
        O: PosType,
        F: Fn(&<Self::AstarStateType<G, D, O> as AstarState<G>>::AstarItem) -> usize;

    fn mismatch(&self) -> u8;

    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;

    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;
}
