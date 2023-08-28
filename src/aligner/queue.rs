use std::collections::VecDeque;
use crate::aligner::extend::{ExtendedPath, ExtendHit};
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, TreeIndexType};
use crate::graphs::NodeIndexType;


#[derive(Copy, Clone, Debug)]
pub enum PathShift<O> {
    Right(O),
    Down(O),
}


#[derive(Clone, Debug, Default)]
pub struct ScoreLayer<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
{
    end_points: Vec<Ix>,
    ext_hits_continued: Vec<(N, ExtendHit<O>)>,
}

impl<N, O, Ix> ScoreLayer<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    pub fn new() -> Self {
        Self {
            end_points: Vec::default(),
            ext_hits_continued: Vec::default(),
        }
    }

    pub fn from_prev(other: &Self) -> Self {
        Self {
            end_points: Vec::with_capacity(other.end_points.len()),
            ext_hits_continued: Vec::with_capacity(other.ext_hits_continued.len()),
        }
    }

    pub fn with_endpoints(endpoints: Vec<Ix>) -> Self {
        Self {
            end_points: endpoints,
            ext_hits_continued: Vec::default(),
        }
    }

    pub fn queue_endpoint(&mut self, end_point: Ix) {
        self.end_points.push(end_point);
    }

    pub fn queue_ext_hit(&mut self, node: N, ext_hit: ExtendHit<O>) {
        self.ext_hits_continued.push((node, ext_hit))
    }

    pub fn queue_additional(&mut self, additional: Vec<Ix>) {
        self.end_points.extend(additional.into_iter())
    }

    pub fn endpoints(&self) -> &[Ix] {
        &self.end_points
    }

    pub fn endpoints_mut(&mut self) -> &mut [Ix] {
        &mut self.end_points
    }

    pub fn extend_hits(&self) -> &[(N, ExtendHit<O>)] {
        &self.ext_hits_continued
    }

    pub fn is_empty(&self) -> bool {
        self.end_points.is_empty() && self.ext_hits_continued.is_empty()
    }
}



/// Score-based queue of alignment state tree nodes to expand
pub struct AlignStateQueue<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    queue: VecDeque<ScoreLayer<N, O, Ix>>
}

impl<N, O, Ix> AlignStateQueue<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
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

    pub fn queue_ext_hit(&mut self, score_delta: u8, node: N, ext_hit: ExtendHit<O>) {
        if self.queue.len() <= score_delta as usize {
            if let Some(last) = self.queue.back() {
                self.queue.resize(score_delta as usize + 1, ScoreLayer::from_prev(last));
            } else {
                self.queue.resize(score_delta as usize + 1, ScoreLayer::default());
            }
        }

        self.queue[score_delta as usize].queue_ext_hit(node, ext_hit);
    }

    pub fn pop_current(&mut self) -> Option<ScoreLayer<N, O, Ix>> {
        self.queue.pop_front()
    }

    pub fn next(&mut self) {
        self.queue.pop_front();
    }
}