use std::fmt;
use std::hash::Hash;

use crate::aligner::utils::print_alignment;
use crate::errors::PoastaError;
use heuristic::AstarHeuristic;
use super::utils::AlignedPair;
use super::{fr_points::Score, AlignmentMode};

pub mod queue;
pub mod heuristic;

/// Enum representing the alignment state of a particular cell in the alignment matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignState {
    Match,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
}


pub trait AlignableGraphNodeId: fmt::Debug + Clone + Copy + PartialEq + Eq + Hash {
    fn index(&self) -> usize;
}

pub trait AlignableGraphNodePos: fmt::Debug + Clone + Copy + PartialEq + Eq {
    type NodeType: AlignableGraphNodeId;
    
    fn new(node: Self::NodeType, pos: usize) -> Self;
    fn node(&self) -> Self::NodeType;
    fn pos(&self) -> usize;
}


pub trait AlignableGraphRef: Copy {
    type NodeType: AlignableGraphNodeId;
    type NodePosType: AlignableGraphNodePos<NodeType = Self::NodeType>;
    type Successors: Iterator<Item=Self::NodeType>;
    type Predecessors: Iterator<Item=Self::NodeType>;

    fn start_node(&self) -> Self::NodeType;
    fn end_node(&self) -> Self::NodeType;

    fn node_count(&self) -> usize;
    fn node_bound(&self) -> usize;

    fn node_seq(&self, node: Self::NodeType) -> &[u8];

    fn node_len(&self, node: Self::NodeType) -> usize {
        self.node_seq(node).len()
    }
    
    fn get_node_symbol(&self, p: Self::NodePosType) -> u8 {
        self.node_seq(p.node())[p.pos()]
    }

    fn successors(&self, node: Self::NodeType) -> Self::Successors;
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors;
}


pub trait AstarState<G> 
where
    G: AlignableGraphRef,
{
    type AstarItem: fmt::Debug;

    fn pop_front(&mut self) -> Option<Self::AstarItem>;

    fn is_further(&self, item: &Self::AstarItem, offset: usize) -> bool;
    fn is_end(&self, item: &Self::AstarItem) -> bool;
    
    fn get_score(&self, item: &Self::AstarItem) -> Score;
    fn get_offset(&self, item: &Self::AstarItem) -> usize;
    
    fn is_visited(&self, item: &Self::AstarItem) -> bool;
    fn set_visited(&mut self, item: &Self::AstarItem);

    fn update_if_further(&mut self, item: &Self::AstarItem, offset: usize) -> bool;

    fn queue_item(&mut self, item: Self::AstarItem, heuristic: usize);

    fn relax<F>(&mut self, seq: &[u8], item: &Self::AstarItem, heuristic: F)
    where
        F: Fn(&Self::AstarItem) -> usize;
    
    fn backtrace(&self, end: &Self::AstarItem) -> Vec<AlignedPair<G::NodePosType>>;
}


pub(crate) struct Astar<'a, 'b, G, H, AS> {
    graph: G,
    seq: &'a [u8],
    heuristic: &'b H,
    alignment_mode: AlignmentMode,
    state: AS,
}

impl<'a, 'b, G, H, AS> Astar<'a, 'b, G, H, AS>
where
    G: AlignableGraphRef,
    H: AstarHeuristic<
        G,
        <AS as AstarState<G>>::AstarItem
    >,
    AS: AstarState<G>,
{
    pub fn new(graph: G, seq: &'a [u8], heuristic: &'b H, alignment_mode: AlignmentMode, state: AS) -> Self {
        Self {
            graph,
            seq,
            heuristic,
            alignment_mode,
            state,
        }
    }

    pub fn run(&mut self) -> Result<AstarResult<G>, PoastaError> {
        let mut result = AstarResult::default();

        let (end_score, end_point) = loop {
            let Some(front) = self.state.pop_front() else {
                panic!("Empty queue before reaching end!")
            };
            
            if self.state.is_end(&front) {
                break (self.state.get_score(&front), front);
            }
            
            if self.state.is_visited(&front) {
                continue;
            }

            eprintln!("--- [FRONT]: {:?}", front);
            self.state.set_visited(&front);
            self.state.relax(self.seq, &front, |item| self.heuristic.h(item));
            
            result.num_visited += 1;
            
            eprintln!();
        };
        
        eprintln!("END: {end_score:?}, {end_point:?}");
        result.score = end_score;
        result.alignment = self.state.backtrace(&end_point);
        
        let aln_str = print_alignment(self.graph, self.seq, &result.alignment);
        eprintln!("{aln_str}");

        Ok(result)
    }
}

pub struct AstarResult<G>
where
    G: AlignableGraphRef,
{
    pub score: Score,
    pub alignment: Vec<AlignedPair<G::NodePosType>>,

    pub num_visited: usize,
    pub num_pruned: usize,
}

impl<G> Default for AstarResult<G>
where
    G: AlignableGraphRef,
{
    fn default() -> Self {
        Self {
            score: Score::default(),
            alignment: Vec::default(),

            num_visited: 0,
            num_pruned: 0
        }
    }
}
