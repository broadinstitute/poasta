use std::collections::VecDeque;
use std::fmt::Debug;

pub trait QueueLayer : Default + Clone {
    type QueueItem: Debug;

    fn queue(&mut self, item: Self::QueueItem);

    fn pop(&mut self) -> Option<Self::QueueItem>;

    fn is_empty(&self) -> bool;

    fn capacity(&self) -> usize;
}


/// Layered (bucket) queue of items
///
/// The layer number represents the priority of the queued item.
pub struct LayeredQueue<L> {
    layers: VecDeque<L>,
    layer_min: usize,
    reusable_layers: Vec<L>,
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
                    self.layers.push_front(self.reusable_layers.pop().unwrap());
                    num_added += 1;
                }

                for _ in 0..(diff-num_added) {
                    self.layers.push_front(L::default())
                }

                self.layer_min = priority;
            } else if priority >= layer_max {
                let diff = priority - layer_max + 1;
                for _ in 0..std::cmp::min(self.reusable_layers.len(), diff) {
                    self.layers.push_back(self.reusable_layers.pop().unwrap());
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
                if popped_layer.capacity() > 32 {
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
        }
    }
}
