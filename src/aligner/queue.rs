use std::collections::VecDeque;
use crate::aligner::extend::ExtendedPath;
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, TreeIndexType};
use crate::graphs::NodeIndexType;


#[derive(Copy, Clone, Debug)]
pub enum PathShift<O> {
    Right(O),
    Down(O),
}

#[derive(Clone, Debug)]
pub struct VisitedPath<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    pub start_node: N,
    pub start_offset: O,
    pub length: O,
    pub shift: PathShift<O>
}

impl<N, O> VisitedPath<N, O>
where
    N: NodeIndexType,
    O: OffsetType
{
    pub fn from_extended_path<Ix>(path: &ExtendedPath<N, O, Ix>, new_state: AlignState) -> Self
    where
        Ix: TreeIndexType,
    {
        let shift = match new_state {
            AlignState::Insertion | AlignState::Insertion2 =>
                PathShift::Right(O::one()),
            AlignState::Deletion | AlignState::Deletion2 =>
                PathShift::Down(O::one()),
            AlignState::Start | AlignState::Match | AlignState::Mismatch =>
                panic!("Invalid state!")
        };

        Self {
            start_node: path.start_node(),
            start_offset: path.start_offset(),
            length: path.len(),
            shift,
        }
    }

    pub fn to_next(&self) -> Self {
        let mut new = self.clone();

        new.shift = match self.shift {
            PathShift::Right(offset) => PathShift::Right(offset.increase_one()),
            PathShift::Down(offset) => PathShift::Down(offset.increase_one()),
        };

        new
    }
}


#[derive(Clone, Debug, Default)]
pub struct ScoreLayer<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
{
    end_points: Vec<Ix>,
    visited_paths: Vec<VisitedPath<N, O>>,
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
            visited_paths: Vec::default(),
        }
    }

    pub fn from_prev(other: &Self) -> Self {
        Self {
            end_points: Vec::with_capacity(other.end_points.len()),
            visited_paths: Vec::with_capacity(other.visited_paths.len()),
        }
    }

    pub fn with_endpoints(endpoints: Vec<Ix>) -> Self {
        Self {
            end_points: endpoints,
            visited_paths: Vec::default(),
        }
    }

    pub fn queue_endpoint(&mut self, end_point: Ix) {
        self.end_points.push(end_point);
    }

    pub fn queue_visited_path(&mut self, path: VisitedPath<N, O>) {
        self.visited_paths.push(path)
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

    pub fn visited_paths(&self) -> &[VisitedPath<N, O>] {
        &self.visited_paths
    }

    pub fn is_empty(&self) -> bool {
        self.end_points.is_empty() && self.visited_paths.is_empty()
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

    pub fn queue_visited_path(&mut self, score_delta: u8, visited_path: VisitedPath<N, O>) {
        if self.queue.len() <= score_delta as usize {
            if let Some(last) = self.queue.back() {
                self.queue.resize(score_delta as usize + 1, ScoreLayer::from_prev(last));
            } else {
                self.queue.resize(score_delta as usize + 1, ScoreLayer::default());
            }
        }

        self.queue[score_delta as usize].queue_visited_path(visited_path);
    }

    pub fn pop_current(&mut self) -> Option<ScoreLayer<N, O, Ix>> {
        self.queue.pop_front()
    }

    pub fn next(&mut self) {
        self.queue.pop_front();
    }
}