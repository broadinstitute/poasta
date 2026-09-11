//! Priority (bucket) queue for the A* engine, ported from `stable`'s
//! `aligner::queue` + the queue-layer half of `aligner::scoring::gap_affine`.

use std::collections::VecDeque;
use std::fmt::Debug;

use crate::align::AlignState;
use crate::align::engine::astar::aln_graph::AlignmentGraphNode;
use crate::align::offsets::OffsetType;
use crate::graph::traits::GraphNodeId;

pub trait QueueLayer: Default + Clone {
    type QueueItem: Debug;

    fn queue(&mut self, item: Self::QueueItem);
    fn pop(&mut self) -> Option<Self::QueueItem>;
    fn is_empty(&self) -> bool;
    fn capacity(&self) -> usize;
}

/// Layered (bucket) queue of items, keyed by priority `f = g + h`.
///
/// LIFO within a bucket (stack-like via each layer's `pop`) preserves the
/// depth-first flavour that pairs with the greedy match-extension DFA.
pub struct LayeredQueue<L> {
    layers: VecDeque<L>,
    layer_min: usize,
}

impl<L> LayeredQueue<L>
where
    L: QueueLayer,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn queue(&mut self, item: L::QueueItem, priority: usize) {
        if self.layers.is_empty() {
            self.layers.push_back(L::default());
            self.layer_min = priority;
        } else {
            let layer_max = self.layer_min + self.layers.len();

            if priority < self.layer_min {
                let diff = self.layer_min - priority;
                self.layers.reserve(diff);

                for _ in 0..diff {
                    self.layers.push_front(L::default())
                }

                self.layer_min = priority;
            } else if priority >= layer_max {
                self.layers.resize(priority - self.layer_min + 1, L::default());
            }
        }

        let ix = priority - self.layer_min;
        self.layers[ix].queue(item);
    }

    pub fn pop(&mut self) -> Option<L::QueueItem> {
        let popped_item = self.layers.front_mut()?.pop();

        while !self.layers.is_empty() {
            if self.layers[0].is_empty() {
                self.layers.pop_front();
                self.layer_min += 1;
            } else {
                break;
            }
        }

        popped_item
    }
}

impl<L> Default for LayeredQueue<L>
where
    L: QueueLayer,
{
    fn default() -> Self {
        Self {
            layers: VecDeque::with_capacity(64),
            layer_min: 0,
        }
    }
}

/// An item popped off the A* queue: its (stale-checked) score at queue time,
/// the alignment-graph node, and the state it was queued in.
#[derive(Clone, Copy, Debug)]
pub struct AstarQueuedItem<N, O>(pub u32, pub AlignmentGraphNode<N, O>, pub AlignState)
where
    N: GraphNodeId,
    O: OffsetType;

impl<N, O> AstarQueuedItem<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    #[inline]
    pub fn score(&self) -> u32 {
        self.0
    }

    #[inline]
    pub fn aln_node(&self) -> AlignmentGraphNode<N, O> {
        self.1
    }

    #[inline]
    pub fn aln_state(&self) -> AlignState {
        self.2
    }
}

pub trait AstarQueue<N, O>: Default
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>>;
    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32, h: usize);
}

/// A queue layer (bucket) for the gap-affine scoring model: keeps queued
/// alignment-graph nodes in different alignment states in separate vectors
/// to reduce branch mispredicts in the main A* loop.
#[derive(Clone)]
pub struct AffineQueueLayer<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    queued_states_m: Vec<(u32, AlignmentGraphNode<N, O>)>,
    queued_states_i: Vec<(u32, AlignmentGraphNode<N, O>)>,
    queued_states_d: Vec<(u32, AlignmentGraphNode<N, O>)>,
}

impl<N, O> QueueLayer for AffineQueueLayer<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    type QueueItem = AstarQueuedItem<N, O>;

    fn queue(&mut self, item: Self::QueueItem) {
        match item.aln_state() {
            AlignState::Match => self.queued_states_m.push((item.score(), item.aln_node())),
            AlignState::Insertion => self.queued_states_i.push((item.score(), item.aln_node())),
            AlignState::Deletion => self.queued_states_d.push((item.score(), item.aln_node())),
            AlignState::Insertion2 | AlignState::Deletion2 => {
                panic!("Invalid gap-affine state {:?}", item.aln_state())
            }
        }
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        self.queued_states_m
            .pop()
            .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Match))
            .or_else(|| {
                self.queued_states_d
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Deletion))
                    .or_else(|| {
                        self.queued_states_i
                            .pop()
                            .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Insertion))
                    })
            })
    }

    fn is_empty(&self) -> bool {
        self.queued_states_m.is_empty() && self.queued_states_d.is_empty() && self.queued_states_i.is_empty()
    }

    fn capacity(&self) -> usize {
        self.queued_states_m.capacity() + self.queued_states_d.capacity() + self.queued_states_i.capacity()
    }
}

impl<N, O> Default for AffineQueueLayer<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn default() -> Self {
        Self {
            queued_states_m: Vec::with_capacity(16),
            queued_states_i: Vec::with_capacity(4),
            queued_states_d: Vec::with_capacity(4),
        }
    }
}

pub type AffineLayeredQueue<N, O> = LayeredQueue<AffineQueueLayer<N, O>>;

impl<N, O> AstarQueue<N, O> for AffineLayeredQueue<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>> {
        self.pop()
    }

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32, h: usize) {
        let priority = score as usize + h;
        let item = AstarQueuedItem(score, node, aln_state);
        self.queue(item, priority)
    }
}

/// A queue layer (bucket) for the two-piece gap-affine scoring model: same
/// branch-reducing idea as [`AffineQueueLayer`], with two extra vectors for
/// the second piece's Insertion2/Deletion2 states.
#[derive(Clone)]
pub struct TwoPieceQueueLayer<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    queued_states_m: Vec<(u32, AlignmentGraphNode<N, O>)>,
    queued_states_i1: Vec<(u32, AlignmentGraphNode<N, O>)>,
    queued_states_i2: Vec<(u32, AlignmentGraphNode<N, O>)>,
    queued_states_d1: Vec<(u32, AlignmentGraphNode<N, O>)>,
    queued_states_d2: Vec<(u32, AlignmentGraphNode<N, O>)>,
}

impl<N, O> QueueLayer for TwoPieceQueueLayer<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    type QueueItem = AstarQueuedItem<N, O>;

    fn queue(&mut self, item: Self::QueueItem) {
        match item.aln_state() {
            AlignState::Match => self.queued_states_m.push((item.score(), item.aln_node())),
            AlignState::Insertion => self.queued_states_i1.push((item.score(), item.aln_node())),
            AlignState::Insertion2 => self.queued_states_i2.push((item.score(), item.aln_node())),
            AlignState::Deletion => self.queued_states_d1.push((item.score(), item.aln_node())),
            AlignState::Deletion2 => self.queued_states_d2.push((item.score(), item.aln_node())),
        }
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        self.queued_states_m
            .pop()
            .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Match))
            .or_else(|| {
                self.queued_states_d1
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Deletion))
            })
            .or_else(|| {
                self.queued_states_d2
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Deletion2))
            })
            .or_else(|| {
                self.queued_states_i1
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Insertion))
            })
            .or_else(|| {
                self.queued_states_i2
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Insertion2))
            })
    }

    fn is_empty(&self) -> bool {
        self.queued_states_m.is_empty()
            && self.queued_states_i1.is_empty()
            && self.queued_states_i2.is_empty()
            && self.queued_states_d1.is_empty()
            && self.queued_states_d2.is_empty()
    }

    fn capacity(&self) -> usize {
        self.queued_states_m.capacity()
            + self.queued_states_i1.capacity()
            + self.queued_states_i2.capacity()
            + self.queued_states_d1.capacity()
            + self.queued_states_d2.capacity()
    }
}

impl<N, O> Default for TwoPieceQueueLayer<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn default() -> Self {
        Self {
            queued_states_m: Vec::with_capacity(16),
            queued_states_i1: Vec::with_capacity(4),
            queued_states_i2: Vec::with_capacity(4),
            queued_states_d1: Vec::with_capacity(4),
            queued_states_d2: Vec::with_capacity(4),
        }
    }
}

pub type TwoPieceLayeredQueue<N, O> = LayeredQueue<TwoPieceQueueLayer<N, O>>;

impl<N, O> AstarQueue<N, O> for TwoPieceLayeredQueue<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>> {
        self.pop()
    }

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32, h: usize) {
        let priority = score as usize + h;
        let item = AstarQueuedItem(score, node, aln_state);
        self.queue(item, priority)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl QueueLayer for Vec<usize> {
        type QueueItem = usize;

        fn queue(&mut self, item: Self::QueueItem) {
            self.push(item);
        }

        fn pop(&mut self) -> Option<Self::QueueItem> {
            Vec::pop(self)
        }

        fn is_empty(&self) -> bool {
            <[usize]>::is_empty(self)
        }

        fn capacity(&self) -> usize {
            Vec::capacity(self)
        }
    }

    #[test]
    fn test_layered_queue() {
        let mut queue: LayeredQueue<Vec<usize>> = LayeredQueue::new();

        queue.queue(10, 10);
        assert_eq!(queue.layers.len(), 1);
        assert_eq!(queue.layer_min, 10);

        queue.queue(5, 5);
        assert_eq!(queue.layers.len(), 6);
        assert_eq!(queue.layer_min, 5);

        let popped = queue.pop();
        assert_eq!(popped, Some(5));

        assert_eq!(queue.layers.len(), 1);
        assert_eq!(queue.layer_min, 10);

        queue.queue(15, 15);
        assert_eq!(queue.layers.len(), 6);
        assert_eq!(queue.layer_min, 10);

        queue.queue(12, 12);
        assert_eq!(queue.layers.len(), 6);
        assert_eq!(queue.layer_min, 10);
        assert_eq!(queue.layers[2], vec![12]);
    }

    /// Popping must always return the lowest-priority item queued so far
    /// (monotone by priority), regardless of insertion order.
    #[test]
    fn pop_order_is_monotone_by_priority() {
        let mut queue: LayeredQueue<Vec<usize>> = LayeredQueue::new();

        let insertions = [(42, 7), (1, 3), (100, 20), (2, 3), (0, 1)];
        for (item, priority) in insertions {
            queue.queue(item, priority);
        }

        let mut popped_priorities = Vec::new();
        // Re-derive priority by tracking which bucket we popped last is hard
        // without extra bookkeeping, so instead assert the values pop out
        // grouped in non-decreasing priority order by construction: items at
        // priority 1, then 3 (LIFO: 2 before 1), then 20.
        while let Some(v) = queue.pop() {
            popped_priorities.push(v);
        }

        assert_eq!(popped_priorities, vec![0, 2, 1, 42, 100]);
    }
}
