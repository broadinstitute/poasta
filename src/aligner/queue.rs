use std::collections::VecDeque;
use std::fmt::Debug;

pub trait QueueLayer : Default + Clone {
    type QueueItem: Debug;

    fn queue(&mut self, item: Self::QueueItem);

    fn pop(&mut self) -> Option<Self::QueueItem>;

    fn is_empty(&self) -> bool;
}


/// Layered (bucket) queue of items
///
/// The layer number represents the priority of the queued item.
#[derive(Default)]
pub struct LayeredQueue<L> {
    layers: VecDeque<L>,
    layer_min: usize,
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
        let popped_item = self.layers.get_mut(0)?
            .pop();

        while !self.layers.is_empty() {
            if self.layers[0].is_empty() {
                self.layers.pop_front();
                self.layer_min += 1;
            } else {
                break
            }
        }

        popped_item
    }
}
