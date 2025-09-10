
use std::sync::Arc;

use crate::aligner::traits::AlignableGraph;
use crate::graph::bubbles::index::BubbleIndex;
use super::astar::{AlignState, AstarState};
use super::fr_points::{Diag, DiagType, PosType, Score};
use super::AlignmentMode;

pub mod affine;


pub trait AstarItem<D>: Clone + std::fmt::Debug
where 
    D: DiagType,
{
    fn new(score: Score, node_rank: usize, node_diag: Diag<D>, aln_state: AlignState) -> Self;
    
    fn score(&self) -> Score;
    fn node_rank(&self) -> usize;
    fn node_diag(&self) -> Diag<D>;
    fn aln_state(&self) -> AlignState;
}


pub trait AlignmentCostModel {
    type DiagType: DiagType;
    type PosType: PosType;
    
    type Item: AstarItem<Self::DiagType>;
    
    type AstarStateType<G>: AstarState<G, Self::Item>
    where
        G: AlignableGraph;
    
    fn init_astar<G>(&self, 
        graph: &G, 
        seq: &[u8], 
        bubble_index: Arc<BubbleIndex<G::NodeType>>, 
        mode: AlignmentMode
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
