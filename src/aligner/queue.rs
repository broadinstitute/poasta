use std::collections::VecDeque;
use crate::aligner::state:: TreeIndexType;


#[derive(Copy, Clone, Debug)]
pub enum PathShift<O> {
    Right(O),
    Down(O),
}


#[derive(Clone, Debug, Default)]
pub struct ScoreLayer<Ix> {
    end_points: Vec<Ix>,
}

impl<Ix> ScoreLayer<Ix>
where
    Ix: TreeIndexType,
{
    pub fn new() -> Self {
        Self {
            end_points: Vec::default(),
        }
    }

    pub fn from_prev(other: &Self) -> Self {
        Self {
            end_points: Vec::with_capacity(other.end_points.len()),
        }
    }

    pub fn with_endpoints(endpoints: Vec<Ix>) -> Self {
        Self {
            end_points: endpoints,
        }
    }

    pub fn queue_endpoint(&mut self, end_point: Ix) {
        self.end_points.push(end_point);
    }

    pub fn queue_additional<T: IntoIterator<Item=Ix>>(&mut self, additional: T) {
        self.end_points.extend(additional.into_iter())
    }

    pub fn endpoints(&self) -> &[Ix] {
        &self.end_points
    }

    pub fn endpoints_mut(&mut self) -> &mut [Ix] {
        &mut self.end_points
    }

    pub fn is_empty(&self) -> bool {
        self.end_points.is_empty()
    }
}



/// Score-based queue of alignment state tree nodes to expand
pub struct AlignStateQueue<Ix>
where
    Ix: TreeIndexType,
{
    queue: VecDeque<ScoreLayer<Ix>>
}

impl<Ix> AlignStateQueue<Ix>
where
    Ix: TreeIndexType,
{
    pub fn new() -> Self {
        Self { queue: VecDeque::default() }
    }

    pub fn queue_endpoint(&mut self, score_delta: u8, tree_node_ix: Ix) {
        if self.queue.len() <= score_delta as usize {
            if let Some(last) = self.queue.back() {
                self.queue.resize(score_delta as usize + 1, ScoreLayer::from_prev(last));
            } else {
                self.queue.resize(score_delta as usize + 1, ScoreLayer::default());
            }
        }

        self.queue[score_delta as usize].queue_endpoint(tree_node_ix);
    }

    pub fn add_additional(&mut self, additional: Vec<Ix>) {
        if let Some(front) = self.queue.front_mut() {
            front.queue_additional(additional)
        } else {
            self.queue.push_back(ScoreLayer::with_endpoints(additional))
        }
    }

    pub fn pop_current(&mut self) -> Option<ScoreLayer<Ix>> {
        self.queue.pop_front()
    }

    pub fn next(&mut self) {
        self.queue.pop_front();
    }
}