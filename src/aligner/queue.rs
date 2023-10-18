use std::collections::VecDeque;
use std::marker::PhantomData;
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::StateGraphNode;
use crate::graphs::NodeIndexType;



#[derive(Clone, Debug)]
pub struct ScoreLayer<S, N, O> {
    states: Vec<S>,
    dummy: PhantomData<(N, O)>
}

impl<S, N, O> ScoreLayer<S, N, O>
where
    S: StateGraphNode<N, O>,
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new() -> Self {
        Self {
            states: Vec::default(),
            dummy: PhantomData,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            states: Vec::with_capacity(capacity),
            dummy: PhantomData,
        }
    }

    pub fn with_states(states: Vec<S>) -> Self {
        Self {
            states,
            dummy: PhantomData
        }
    }

    pub fn queue_state(&mut self, state: S) {
        self.states.push(state);
    }

    pub fn queue_additional<T: IntoIterator<Item=S>>(&mut self, additional: T) {
        self.states.extend(additional.into_iter())
    }

    pub fn states(&self) -> &[S] {
        &self.states
    }

    pub fn states_mut(&mut self) -> &mut [S] {
        &mut self.states
    }

    pub fn iter(&self) -> std::slice::Iter<'_, S> {
        self.states.iter()
    }

    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }
}

impl<S, N, O> Default for ScoreLayer<S, N, O>
where
    S: StateGraphNode<N, O>,
    O: OffsetType,
    N: NodeIndexType,
{
    fn default() -> Self {
        Self {
            states: Vec::default(),
            dummy: PhantomData,
        }
    }
}

impl<S, N, O> IntoIterator for ScoreLayer<S, N, O>
where
    S: StateGraphNode<N, O>,
    N: NodeIndexType,
    O: OffsetType,
{
    type Item = S;
    type IntoIter = <Vec<S> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.states.into_iter()
    }
}

impl<'a, S, N, O> IntoIterator for &'a ScoreLayer<S, N, O>
where
    S: StateGraphNode<N, O>,
    N: NodeIndexType,
    O: OffsetType,
{
    type Item = &'a S;
    type IntoIter = std::slice::Iter<'a, S>;

    fn into_iter(self) -> Self::IntoIter {
        self.states.iter()
    }
}


/// Layered, score-based queue of alignment state tree nodes to explore
pub struct AlignStateQueue<S, N, O> {
    current: ScoreLayer<S, N, O>,
    queue: VecDeque<ScoreLayer<S, N, O>>
}

impl<S, N, O> AlignStateQueue<S, N, O>
where
    S: StateGraphNode<N, O>,
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new() -> Self {
        Self {
            current: ScoreLayer::default(),
            queue: VecDeque::default()
        }
    }

    pub fn queue_state(&mut self, state: S, score_delta: u8) {
        if score_delta == 0 {
            self.current.queue_state(state)
        } else {
            if self.queue.len() < score_delta as usize {
                let capacity = if let Some(last) = self.queue.back() {
                    last.len()
                } else {
                    1
                };

                self.queue
                    .resize(score_delta as usize, ScoreLayer::with_capacity(capacity));
            }

            self.queue[score_delta as usize - 1].queue_state(state)
        }
    }

    pub fn pop_front(&mut self) -> Option<ScoreLayer<S, N, O>> {
        if self.current.is_empty() && self.queue.is_empty() {
            None
        } else {
            // Replace current with an empty layer, to ensure queuing new states with a certain
            // score delta remains intuitive.
            Some(std::mem::take(&mut self.current))
        }
    }

    pub fn next_score(&mut self) {
        if !self.current.is_empty() {
            panic!("Can't move to next score: front queue is not empty!")
        }

        if let Some(front) = self.queue.pop_front() {
            let _ = std::mem::replace(&mut self.current, front);
        }
    }
}

impl<S, N, O> Default for AlignStateQueue<S, N, O>
where
    S: StateGraphNode<N, O>,
    N: NodeIndexType,
    O: OffsetType,
{
    fn default() -> Self {
        Self::new()
    }
}
