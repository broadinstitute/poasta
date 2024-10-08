use std::fmt;
use std::hash::Hash;
use std::marker::PhantomData;

use crate::errors::PoastaError;
use heuristic::AstarHeuristic;
use tracing::{debug, debug_span, span, Level};
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


pub trait AlignableGraph {
    type NodeType: AlignableGraphNodeId;
    type NodePosType: AlignableGraphNodePos<NodeType = Self::NodeType>;
    
    type Successors<'a>: Iterator<Item=Self::NodeType> + 'a
        where Self: 'a;
    type Predecessors<'a>: Iterator<Item=Self::NodeType> + 'a
        where Self: 'a;

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

    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_>;
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_>;
}


pub trait AstarState<G> 
where
    G: AlignableGraph,
{
    type AstarItem: fmt::Debug;

    fn pop_front(&mut self) -> Option<Self::AstarItem>;

    fn is_further(&self, item: &Self::AstarItem, offset: usize) -> bool;
    fn is_end(&self, graph: &G, item: &Self::AstarItem) -> bool;
    
    fn get_score(&self, item: &Self::AstarItem) -> Score;
    fn get_offset(&self, item: &Self::AstarItem) -> usize;
    
    fn is_visited(&self, item: &Self::AstarItem) -> bool;
    fn set_visited(&mut self, item: &Self::AstarItem);

    fn update_if_further(&mut self, item: &Self::AstarItem, offset: usize) -> bool;

    fn queue_item(&mut self, item: Self::AstarItem, heuristic: usize);

    fn relax<F>(&mut self, graph: &G, seq: &[u8], item: &Self::AstarItem, heuristic: F)
    where
        F: Fn(&Self::AstarItem) -> usize;
    
    fn backtrace(&self, graph: &G, end: &Self::AstarItem) -> Vec<AlignedPair<G::NodePosType>>;
}


pub(crate) struct Astar<'a, 'b, 'c, G, H, AS> {
    graph: &'a G,
    seq: &'b [u8],
    heuristic: &'c H,
    alignment_mode: AlignmentMode,
    state: AS,
    dummy: PhantomData<G>,
}

impl<'a, 'b, 'c, G, H, AS> Astar<'a, 'b, 'c, G, H, AS>
where
    G: AlignableGraph,
    H: AstarHeuristic<
        G,
        <AS as AstarState<G>>::AstarItem
    >,
    AS: AstarState<G>,
{
    pub fn new(graph: &'a G, seq: &'b [u8], heuristic: &'c H, alignment_mode: AlignmentMode, state: AS) -> Self {
        Self {
            graph,
            seq,
            heuristic,
            alignment_mode,
            state,
            dummy: PhantomData,
        }
    }

    pub fn run(&mut self) -> Result<AstarResult<G>, PoastaError> {
        let span = span!(Level::INFO, "astar_run");
        let _enter = span.enter();
        
        let mut result = AstarResult::default();

        let (end_score, end_point) = loop {
            let Some(front) = self.state.pop_front() else {
                panic!("Empty queue before reaching end!")
            };
            
            if self.state.is_end(self.graph, &front) {
                self.state.set_visited(&front);
                break (self.state.get_score(&front), front);
            }
            
            if self.state.is_visited(&front) {
                continue;
            }
            
            debug!("--- FRONT {:?}", front);
            
            self.state.set_visited(&front);
            result.num_visited += 1;
            
            self.state.relax(self.graph, self.seq, &front, |item| self.heuristic.h(item));
        };
        
        debug!(score = end_score.as_usize(), ?end_point, "END");
        
        result.score = end_score;
        result.alignment = self.state.backtrace(self.graph, &end_point);

        Ok(result)
    }
}

pub struct AstarResult<G>
where
    G: AlignableGraph,
{
    pub score: Score,
    pub alignment: Vec<AlignedPair<G::NodePosType>>,

    pub num_visited: usize,
    pub num_pruned: usize,
}

impl<G> Default for AstarResult<G>
where
    G: AlignableGraph,
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
