use std::collections::VecDeque;
use std::fmt::Debug;

pub trait QueueLayer : Default + Clone {
    type QueueItem: Debug;

    fn queue(&mut self, item: Self::QueueItem);

    fn pop(&mut self) -> Option<Self::QueueItem>;

    fn is_empty(&self) -> bool;

    fn capacity(&self) -> usize;
}

const MAX_REUSABLE_CAPACITY: usize = 4096;

/// Layered (bucket) queue of items
///
/// The layer number represents the priority of the queued item.
pub struct LayeredQueue<L> {
    layers: VecDeque<L>,
    layer_min: usize,
    reusable_layers: Vec<L>,
    reusable_capacity: usize,
}

impl<L> LayeredQueue<L>
    where L: QueueLayer
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

                let mut num_added = 0;
                for _ in 0..std::cmp::min(self.reusable_layers.len(), diff) {
                    let new_layer = self.reusable_layers.pop().unwrap();
                    self.reusable_capacity -= new_layer.capacity();
                    self.layers.push_front(new_layer);
                    num_added += 1;
                }

                for _ in 0..(diff-num_added) {
                    self.layers.push_front(L::default())
                }

                self.layer_min = priority;
            } else if priority >= layer_max {
                let diff = priority - layer_max + 1;
                for _ in 0..std::cmp::min(self.reusable_layers.len(), diff) {
                    let new_layer = self.reusable_layers.pop().unwrap();
                    self.reusable_capacity -= new_layer.capacity();
                    self.layers.push_back(new_layer);
                }

                self.layers.resize(priority - self.layer_min + 1, L::default());
            }
        }

        let ix = priority - self.layer_min;
        self.layers[ix].queue(item);
    }

    pub fn pop(&mut self) -> Option<L::QueueItem> {
        let popped_item = self.layers.get_mut(0)?
            .pop();

        while !self.layers.is_empty() {
            if self.layers[0].is_empty() {
                let popped_layer = self.layers.pop_front().unwrap();

                // Retain empty layers because we can reuse our hard-fought allocated memory
                // for new layers
                if popped_layer.capacity() > 32 && self.reusable_capacity < MAX_REUSABLE_CAPACITY {
                    self.reusable_capacity += popped_layer.capacity();
                    self.reusable_layers.push(popped_layer)
                }
                self.layer_min += 1;
            } else {
                break
            }
        }

        popped_item
    }
}

impl<L> Default for LayeredQueue<L>
    where L: QueueLayer,
{
    fn default() -> Self {
        Self {
            layers: VecDeque::with_capacity(1024),
            layer_min: 0,
            reusable_layers: Vec::with_capacity(32),
            reusable_capacity: 0,
        }
    }
}
