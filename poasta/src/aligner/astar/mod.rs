use std::fmt;
use std::marker::PhantomData;

use crate::errors::PoastaError;
use heuristic::AstarHeuristic;
use tracing::{debug, span, Level};

use super::utils::AlignedPair;
use super::{fr_points::Score, AlignmentMode};

pub mod queue;
pub mod heuristic;

use crate::aligner::traits::AlignableGraph;

/// Enum representing the alignment state of a particular cell in the alignment matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignState {
    Match,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
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
