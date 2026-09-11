//! Priority (bucket) queue for the A* engine.
//!
//! Priorities are small integers (`score + heuristic`), so the queue is a
//! deque of buckets rather than a heap: pushing and popping are O(1) and the
//! bucket for the current minimum is always at the front.

use std::collections::VecDeque;
use std::fmt::Debug;

use crate::align::AlignState;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::offsets::OffsetType;
use crate::graph::traits::GraphNodeId;

pub trait QueueLayer: Default + Clone {
    type QueueItem: Debug;

    fn queue(&mut self, item: Self::QueueItem);
    fn pop(&mut self) -> Option<Self::QueueItem>;
    fn is_empty(&self) -> bool;
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

impl<N, O, L> AstarQueue<N, O> for LayeredQueue<L>
where
    N: GraphNodeId,
    O: OffsetType,
    L: QueueLayer<QueueItem = AstarQueuedItem<N, O>>,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>> {
        self.pop()
    }

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32, h: usize) {
        let priority = score as usize + h;
        tracing::trace!(
            target: "poasta::astar::queue_item",
            score,
            node = node.node().index(),
            offset = node.offset().as_usize(),
            state = ?aln_state,
            priority,
        );
        self.queue(AstarQueuedItem(score, node, aln_state), priority)
    }
}

/// Pop order within one bucket: matches first, then deletions, then
/// insertions, piece one before piece two. States whose index is `>= S` are
/// skipped, so the same order serves every state count.
const POP_ORDER: [AlignState; 5] = [
    AlignState::Match,
    AlignState::Deletion,
    AlignState::Deletion2,
    AlignState::Insertion,
    AlignState::Insertion2,
];

/// One bucket of the layered queue. States in different alignment states go
/// in separate vectors, which keeps the main loop's state dispatch
/// predictable. `S` is the cost model's state count: 1 for linear, 3 for
/// affine, 5 for two-piece affine.
#[derive(Clone)]
pub struct StateQueueLayer<N, O, const S: usize>
where
    N: GraphNodeId,
    O: OffsetType,
{
    states: [Vec<(u32, AlignmentGraphNode<N, O>)>; S],
}

impl<N, O, const S: usize> Default for StateQueueLayer<N, O, S>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn default() -> Self {
        Self {
            states: std::array::from_fn(|i| Vec::with_capacity(if i == 0 { 16 } else { 4 })),
        }
    }
}

impl<N, O, const S: usize> QueueLayer for StateQueueLayer<N, O, S>
where
    N: GraphNodeId,
    O: OffsetType,
{
    type QueueItem = AstarQueuedItem<N, O>;

    fn queue(&mut self, item: Self::QueueItem) {
        self.states[item.aln_state() as usize].push((item.score(), item.aln_node()));
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        for state in POP_ORDER {
            let ix = state as usize;
            if ix >= S {
                continue;
            }
            if let Some((score, node)) = self.states[ix].pop() {
                return Some(AstarQueuedItem(score, node, state));
            }
        }
        None
    }

    fn is_empty(&self) -> bool {
        self.states.iter().all(Vec::is_empty)
    }
}

pub type LinearLayeredQueue<N, O> = LayeredQueue<StateQueueLayer<N, O, 1>>;
pub type AffineLayeredQueue<N, O> = LayeredQueue<StateQueueLayer<N, O, 3>>;
pub type TwoPieceLayeredQueue<N, O> = LayeredQueue<StateQueueLayer<N, O, 5>>;

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

    /// Popping must always return the lowest-priority item queued so far,
    /// regardless of insertion order; within a bucket the order is LIFO.
    #[test]
    fn pop_order_is_monotone_by_priority() {
        let mut queue: LayeredQueue<Vec<usize>> = LayeredQueue::new();

        let insertions = [(42, 7), (1, 3), (100, 20), (2, 3), (0, 1)];
        for (item, priority) in insertions {
            queue.queue(item, priority);
        }

        let mut popped = Vec::new();
        while let Some(v) = queue.pop() {
            popped.push(v);
        }

        assert_eq!(popped, vec![0, 2, 1, 42, 100]);
    }

    #[test]
    fn state_layer_pops_matches_before_deletions_before_insertions() {
        use crate::graph::mock::NIx;

        let mut layer: StateQueueLayer<NIx, u32, 5> = StateQueueLayer::default();
        let node = AlignmentGraphNode::new(NIx::new(0), 0u32);
        for state in [
            AlignState::Insertion2,
            AlignState::Insertion,
            AlignState::Deletion2,
            AlignState::Deletion,
            AlignState::Match,
        ] {
            layer.queue(AstarQueuedItem(0, node, state));
        }

        let order: Vec<AlignState> = std::iter::from_fn(|| layer.pop().map(|i| i.aln_state())).collect();
        assert_eq!(order, POP_ORDER.to_vec());
        assert!(layer.is_empty());
    }
}
