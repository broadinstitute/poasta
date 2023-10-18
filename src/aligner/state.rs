use super::offsets::OffsetType;
use crate::graphs::{AlignableGraph, NodeIndexType};
use std::cmp::Ordering;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::io::Write;
use std::ops::{Add, AddAssign};
use std::cell::{RefCell, RefMut};
use std::marker::PhantomData;
use crate::aligner::Alignment;
use crate::aligner::queue::AlignStateQueue;
use crate::bubbles::index::BubbleIndex;
use crate::bubbles::ReachedBubbleExits;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignState {
    Start,
    Match,
    Mismatch,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
}


#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Score {
    Score(usize),
    Unvisited
}

impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self {
            Self::Score(score) => match other {
                Self::Score(other_score) => score.partial_cmp(other_score),
                Self::Unvisited => Some(Ordering::Less),
            },
            Self::Unvisited => match other {
                Self::Score(_) => Some(Ordering::Greater),
                Self::Unvisited => Some(Ordering::Equal)
            }
        }
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Add<usize> for Score {
    type Output = Self;

    fn add(self, rhs: usize) -> Self::Output {
        match self {
            Self::Score(score) => Self::Score(score + rhs),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl AddAssign<usize> for Score {
    fn add_assign(&mut self, rhs: usize) {
        match self {
            Self::Score(score) => *score += rhs,
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl Add<u8> for Score {
    type Output = Self;

    fn add(self, rhs: u8) -> Self::Output {
        match self {
            Self::Score(score) => Self::Score(score + rhs as usize),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl AddAssign<u8> for Score {
    fn add_assign(&mut self, rhs: u8) {
        match self {
            Self::Score(score) => *score += rhs as usize,
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl From<Score> for usize {
    fn from(value: Score) -> Self {
        match value {
            Score::Score(score) => score,
            Score::Unvisited => panic!("Trying to convert Score::Unvisited!")
        }
    }
}

impl Display for Score {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Score(score) => Display::fmt(score, f),
            Self::Unvisited => Display::fmt("unvisited", f)
        }
    }
}

pub trait StateGraphNode<N, O>: Debug + Clone {
    fn new(node: N, offset: O, state: AlignState) -> Self;

    fn node(&self) -> N;
    fn offset(&self) -> O;
    fn state(&self) -> AlignState;
}

pub trait StateGraph<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    type StateNode: StateGraphNode<N, O>;
    type NewStatesContainer: IntoIterator<Item=(Self::StateNode, u8)>;

    fn get_score(&self, state: &Self::StateNode) -> Score;
    fn update_score(&mut self, state: &Self::StateNode, score: Score, prev: &Self::StateNode);

    fn get_prev(&self, state: &Self::StateNode) -> Option<&Self::StateNode>;

    fn new_match_state(&mut self, parent: &Self::StateNode, child: N, current_score: Score) -> Option<Self::StateNode>;
    fn new_mismatch_state(&mut self, parent: &Self::StateNode, child: N, current_score: Score) -> Option<(Self::StateNode, u8)>;
    fn open_or_extend_insertion(&mut self, parent: &Self::StateNode, score: Score) -> Self::NewStatesContainer;
    fn extend_insertion(&mut self, parent: &Self::StateNode, score: Score) -> Option<(Self::StateNode, u8)>;
    fn open_or_extend_deletion(&mut self, parent: &Self::StateNode, child_node: N, score: Score) -> Self::NewStatesContainer;
    fn extend_deletion(&mut self, parent: &Self::StateNode, child: N, score: Score) -> Option<(Self::StateNode, u8)>;

    fn backtrace(&self, end_state: &Self::StateNode) -> Alignment<N>;

    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>>;
}

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
struct StackNode<S, N, O, I>(S, RefCell<I>, PhantomData<(N, O)>);

impl<S, N, O, I> StackNode<S, N, O, I>
where
    S: StateGraphNode<N, O>,
    I: Iterator,
{
    fn new(state: S, iter: I)  -> Self {
        Self(state, RefCell::new(iter), PhantomData)
    }

    #[inline(always)]
    fn node(&self) -> N {
        self.0.node()
    }

    #[inline(always)]
    fn offset(&self) -> O {
        self.0.offset()
    }

    #[inline(always)]
    fn state(&self) -> AlignState {
        self.0.state()
    }

    #[inline(always)]
    fn state_graph_node(&self) -> &S {
        &self.0
    }

    #[inline]
    fn into_state_graph_node(self) -> S {
        self.0
    }

    #[inline]
    fn children_iter(&self) -> RefMut<I> {
        self.1.borrow_mut()
    }
}

pub enum ExtendResult<S, N> {
    Match(S),
    Mismatch(N),
    None,
}


/// A struct representing the depth-first extension state when
/// extending from a (node, query offset).
///
/// The algorithm is inspired by NetworkX's `dfs_labeled_edges` function, which in turn
/// is based on David Eppstein's depth-first search tree algorithm.
///
/// See also: https://www.ics.uci.edu/~eppstein/PADS/ or
/// https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.traversal.depth_first_search.dfs_labeled_edges.html
pub struct StateGraphSuccessors<'a, G, O, S>
where
    G: AlignableGraph,
    O: OffsetType,
    S: StateGraphNode<G::NodeIndex, O>,
{
    /// Borrow of the (reference) graph to align to
    graph: &'a G,

    /// The found bubbles in the graph
    bubble_index: &'a BubbleIndex<G::NodeIndex, O>,

    /// The query sequence to align
    seq: &'a [u8],

    /// The current alignment score
    score: Score,

    /// Stack for depth-first alignment of matches between query and graph
    stack: Vec<StackNode<S, G::NodeIndex, O, G::SuccessorIterator<'a>>>,

    /// Flag to check if we are at a leaf node in the DFS tree
    has_new_path: bool,
}

impl<'a, G, O, S> StateGraphSuccessors<'a, G, O, S>
where
    G: AlignableGraph,
    O: OffsetType,
    S: StateGraphNode<G::NodeIndex, O>
{
    pub fn new(
        graph: &'a G,
        bubble_index: &'a BubbleIndex<G::NodeIndex, O>,
        seq: &'a [u8],
        score: Score,
        start_state: &S,
    ) -> Self {
        Self {
            graph,
            bubble_index,
            seq,
            score,
            stack: vec![StackNode::new(start_state.clone(), graph.successors(start_state.node()))],
            has_new_path: false,
        }
    }

    pub fn queue_next<SG: StateGraph<G::NodeIndex, O, StateNode=S>>(
        &mut self,
        state_graph: &mut SG,
        queue: &mut AlignStateQueue<S, G::NodeIndex, O>,
        bubble_exits_reached: &mut ReachedBubbleExits<O>,
    ) -> Option<S> {
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();

            match self.next_successor(parent, state_graph, queue) {
                ExtendResult::Match(child_match) => {
                    self.has_new_path = true;

                    if self.bubble_index.is_exit(child_match.node()) {
                        // eprintln!("State {:?} reached bubble exit at score {}", child_match, self.score);

                        bubble_exits_reached[child_match.node().index()]
                            .push((child_match.offset(), self.score));
                    }

                    let child_succ = self.graph.successors(child_match.node());
                    self.stack.push(StackNode::new(child_match, child_succ));
                },
                ExtendResult::Mismatch(child) => {
                    // Queue the mismatch state if it will improve the score
                    if let Some((child_mismatch, score_delta)) = state_graph
                        .new_mismatch_state(parent.state_graph_node(), child, self.score)
                    {
                        // eprintln!("- EXPAND queue MIS {:?} (score: {}+{}={})", child_mismatch, self.score, score_delta, self.score + score_delta);
                        queue.queue_state(child_mismatch, score_delta);
                    }

                    // In case of a mismatch, also queue open/extend indels from parent
                    for (ins_state, score_delta) in state_graph
                        .open_or_extend_insertion(parent.state_graph_node(), self.score)
                    {
                        // eprintln!("- EXPAND queue INS {:?} (score: {}+{}={})", ins_state, self.score, score_delta, self.score + score_delta);
                        queue.queue_state(ins_state, score_delta);
                    }

                    for (del_state, score_delta) in state_graph
                        .open_or_extend_deletion(parent.state_graph_node(), child, self.score)
                    {
                        // eprintln!("- EXPAND queue DEL {:?} (score: {}+{}={})", del_state, self.score, score_delta, self.score + score_delta);
                        queue.queue_state(del_state, score_delta);
                    }

                    self.has_new_path = false;
                },
                ExtendResult::None => {
                    // We are about to move up on the stack
                    // Check if we are at a leaf node in the DFS tree, if yes,
                    // return leaf.
                    let popped = self.stack.pop().unwrap();

                    if self.has_new_path {
                        self.has_new_path = false;

                        return Some(popped.into_state_graph_node());
                    }

                }
            }
        }

        None
    }

    fn next_successor<SG: StateGraph<G::NodeIndex, O, StateNode=S>>(
        &self,
        parent: &StackNode<S, G::NodeIndex, O, G::SuccessorIterator<'a>>,
        // Separate mutable borrows for queue and state graph to prevent requiring mutable &self
        state_graph: &mut SG,
        queue: &mut AlignStateQueue<S, G::NodeIndex, O>,
    ) -> ExtendResult<S, G::NodeIndex> {
        if parent.offset().as_usize() > self.seq.len() {
            return ExtendResult::None;
        }

        while let Some(child) = parent.children_iter().next() {
            match parent.offset().as_usize().cmp(&self.seq.len()) {
                // Make sure we do not extend beyond the query sequence length
                Ordering::Less => {
                    let child_offset = parent.offset().increase_one();
                    let symbol_equal = self.graph
                        .is_symbol_equal(child, self.seq[child_offset.as_usize() - 1]);

                    if symbol_equal {
                        if let Some(child_state) = state_graph
                            .new_match_state(parent.state_graph_node(), child, self.score)
                        {
                            // On match, return to put the new state on the stack.
                            return ExtendResult::Match(child_state)
                        }
                    } else {
                        return ExtendResult::Mismatch(child)
                    }
                },
                Ordering::Equal => {
                    // Here, open/extend deletions from states in the last column, since these
                    // won't move the query offset, and thus will not move beyond the query length.
                    for (del_state, score_delta) in state_graph
                        .open_or_extend_deletion(parent.state_graph_node(), child, self.score)
                    {
                        // eprintln!("- EXPAND queue DEL {:?} (score: {}+{}={})", del_state, self.score, score_delta, self.score + score_delta);
                        queue.queue_state(del_state, score_delta);
                    }

                },
                Ordering::Greater => (),
            }
        }

        ExtendResult::None
    }

}
