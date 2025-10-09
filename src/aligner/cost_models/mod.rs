use std::fmt::Debug;
use std::sync::Arc;

use super::astar::{AlignState, AstarState};
use super::fr_points::{DiagType, PosType};
use super::AlignmentMode;
use crate::aligner::traits::AlignableGraph;
use crate::graph::bubbles::index::BubbleIndex;

pub mod affine;

pub trait AlignmentCostModel {
    type DiagType: DiagType;
    type PosType: PosType;

    type Item: Debug + Clone;

    type AstarStateType<G>: AstarState<G, Self::Item>
    where
        G: AlignableGraph;

    fn init_astar<G>(
        &self,
        graph: &G,
        seq: &[u8],
        bubble_index: Arc<BubbleIndex>,
        mode: AlignmentMode,
    ) -> Self::AstarStateType<G>
    where
        G: AlignableGraph;

    fn mismatch(&self) -> u8;

    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;

    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;

    fn gap_cost(&self, current_state: AlignState, gap_length: usize) -> usize;
}
