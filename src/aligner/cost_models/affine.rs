use std::marker::PhantomData;

use crate::aligner::{astar::{AlignState, AstarAlignableGraph, AstarState}, fr_points::{Diag, IndexType, NodeFrPoints}, AlignmentMode};
use super::AlignmentCostModel;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Affine {
    mismatch: u8,
    gap_open: u8,
    gap_extend: u8,
}

impl Affine {
    pub fn new(mismatch: u8, gap_open: u8, gap_extend: u8) -> Self {
        Self {
            mismatch,
            gap_open,
            gap_extend,
        }
    }
}

impl AlignmentCostModel for Affine {
    type AstarState<N, D> = AffineAstarState<N, D>;
    
    fn initialize<G, D>(&self, graph: &G, alignment_mode: AlignmentMode) -> Self::AstarState<G::NodeType, D>
    where
        G: AstarAlignableGraph,
        D: IndexType 
    {
        let mut state = AffineAstarState::new(graph.node_count());
        
        match alignment_mode {
            AlignmentMode::Global => {
                
            },
            AlignmentMode::EndsFree { qry_free_begin, qry_free_end, graph_free_begin, graph_free_end } => {
                
            }
        }
        
        state
    }
    
    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.mismatch
    }
    
    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.gap_open
    }
    
    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.gap_extend
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AffineAstarItem<N, D> {
    pub node: N,
    pub diag: Diag<D>,
    pub offset: D,
    pub state: AlignState,
}

#[derive(Clone, Debug, Default)]
struct AffineNodeDiagonals<D> {
    fr_points_m: NodeFrPoints<D>,
    fr_points_i: NodeFrPoints<D>,
    fr_points_d: NodeFrPoints<D>,
}


pub(crate) struct AffineAstarState<N, D> {
    fr_points: Vec<AffineNodeDiagonals<D>>,
    node_type: PhantomData<N>
}

impl<N, D> AffineAstarState<N, D>
    where D: IndexType
{
    fn new(node_count: usize) -> Self {
        Self {
            fr_points: vec![AffineNodeDiagonals::default(); node_count],
            node_type: PhantomData
        }
    }
}

impl<N, D> AstarState for AffineAstarState<N, D> {
    type QueueItem = AffineAstarItem<N, D>;
    
}


// struct AffineDiagonals<O, const B: usize> {
//     fr_points_m: [O; B],
//     fr_points_i: [O; B],
//     fr_points_d: [O; B],
// }

// struct BlockedAffineDiagonals<O, const B: usize=16> {
//     /// For each node in the POA graph, we store the query offsets per diagonal. For cache-efficiency,
//     /// we use a blocked storage scheme where we store the diagonals in blocks of size B.
//     fr_point_blocks: Vec<FxHashMap<usize, AffineDiagonals<O, B>>>,
// }

// impl<O, const B: usize> BlockedAffineDiagonals<O, B> {
//     fn new() -> Self {
//         if B & (B-1) != 0 {
//             panic!("BlockedAffineDiagonals block size must be a power of 2");
//         }
        
//         Self {
//             fr_point_blocks: vec![FxHashMap::default()]
//         }
//     }
// }


