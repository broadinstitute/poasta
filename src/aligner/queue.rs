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
            if let Some(last) = self.queue.back() {
                self.queue.resize(score_delta as usize + 1, Vec::with_capacity(last.len()));
            } else {
                self.queue.resize(score_delta as usize + 1, Vec::default());
            }
        }

        self.queue[score_delta as usize].push(tree_node_ix);
    }

    pub fn add_additional(&mut self, additional: Vec<Ix>) {
        if let Some(front) = self.queue.front_mut() {
            front.extend(additional.into_iter())
        } else {
            self.queue.push_back(additional)
        }
    }

    pub fn pop_current(&mut self) -> Option<Vec<Ix>> {
        self.queue.pop_front()
    }

    pub fn next(&mut self) {
        self.queue.pop_front();
    }
}