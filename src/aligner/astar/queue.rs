use std::collections::VecDeque;
use std::fmt::Debug;

pub trait QueueLayer : Default + Clone {
    type QueueItem: Debug;

    fn queue(&mut self, item: Self::QueueItem);

    fn pop(&mut self) -> Option<Self::QueueItem>;

    fn is_empty(&self) -> bool;

    fn capacity(&self) -> usize;
}

impl<T> QueueLayer for Vec<T>
    where T: Clone + Debug
{
    type QueueItem = T;
    
    fn queue(&mut self, item: Self::QueueItem) {
        self.push(item);
    }
    
    fn pop(&mut self) -> Option<Self::QueueItem> {
        self.pop()
    }
    
    fn is_empty(&self) -> bool {
        self.is_empty()
    }
    
    fn capacity(&self) -> usize {
        self.capacity()
    }
}

/// Layered (bucket) queue of items
///
/// The layer number represents the priority of the queued item.
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

impl<L> Default for LayeredQueue<L>
    where L: QueueLayer,
{
    fn default() -> Self {
        Self {
            layers: VecDeque::with_capacity(64),
            layer_min: 0,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::{LayeredQueue, QueueLayer};

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
}