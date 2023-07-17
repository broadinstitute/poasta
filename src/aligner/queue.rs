use std::collections::VecDeque;
use crate::aligner::state::TreeIndexType;

/// Score-based queue of alignment state tree nodes to expand
pub struct AlignStateQueue<Ix: TreeIndexType> {
    queue: VecDeque<Vec<Ix>>
}

impl<Ix: TreeIndexType> AlignStateQueue<Ix> {
    pub fn new() -> Self {
        Self { queue: VecDeque::default() }
    }

    pub fn enqueue(&mut self, score_delta: u8, tree_node_ix: Ix) {
        if self.queue.len() <= score_delta as usize {
            self.queue.resize(score_delta as usize + 1, Vec::default());
        }

        self.queue[score_delta as usize].push(tree_node_ix);
    }

    pub fn pop_current(&mut self) -> Option<Vec<Ix>> {
        self.queue.pop_front()
    }

    pub fn next(&mut self) {
        self.queue.remove(0);
    }
}