use super::fr_points::Score;
use super::utils::AlignedPair;

pub mod heuristic;
pub mod queue;
pub mod runnable;

use crate::align::traits::AlignableGraph;

pub trait AstarState<G, Item>
where
    G: AlignableGraph,
{
    fn pop_front(&mut self) -> Option<Item>;

    fn is_further(&self, item: &Item, offset: usize) -> bool;
    fn is_end(&self, graph: &G, item: &Item) -> bool;

    fn get_score(&self, item: &Item) -> Score;
    fn get_offset(&self, item: &Item) -> usize;

    fn update_if_further(&mut self, item: &Item, offset: usize) -> bool;

    fn queue_item(&mut self, item: Item, heuristic: usize);

    fn relax<F>(&mut self, graph: &G, seq: &[u8], item: &Item, heuristic: F)
    where
        F: Fn(&Self, &Item) -> usize;

    fn backtrace(&self, graph: &G, end: &Item) -> Vec<AlignedPair<G::NodePosType>>;
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
            num_pruned: 0,
        }
    }
}
