use super::offsets::OffsetType;
use crate::graphs::NodeIndexType;
use core::fmt;
use num::{FromPrimitive, Unsigned};
use std::hash::Hash;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignState {
    Start,
    Match,
    Mismatch,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
}

/// Allow for various index types to refer to the alignment state tree nodes
pub trait TreeIndexType: Copy + Unsigned + FromPrimitive + Hash + fmt::Debug + Default {
    fn new(value: usize) -> Self;
    fn index(&self) -> usize;
}

impl TreeIndexType for usize {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value
    }

    fn index(&self) -> usize {
        *self
    }
}

impl TreeIndexType for u32 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as u32
    }

    #[inline(always)]
    fn index(&self) -> usize {
        *self as usize
    }
}

impl TreeIndexType for u64 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as u64
    }

    #[inline(always)]
    fn index(&self) -> usize {
        *self as usize
    }
}

/// The main tree data structure representing the alignment state
///
/// We use "arena pattern" to structure and define ownership of all nodes. This
/// main struct is the owner of all nodes in the tree. Through numeric indices
/// in each [`StateTreeNode`] object we define parent-child relationships.
///
/// The structure has the following generic type parameters:
///
/// * `N` - The type used for indexing nodes in the graph we are aligning to
/// * `O` - The integer used for storing query sequence offsets
/// * `S` - A type that generates next alignment states based on existing nodes in the
///   alignment state tree. This provides support for different alignment scoring schemes.
/// * `Ix` - The integer type used for node indices in the alignment state tree. Defaults to u32.
pub struct StateTree<N, O, Ix = u32>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    /// The vector storing all current nodes in the alignment state tree
    nodes: Vec<StateTreeNode<N, O, Ix>>,
}

impl<N, O, Ix> StateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    pub fn new() -> Self {
        Self {
            nodes: Vec::default(),
        }
    }

    pub fn reserve(&mut self, additional: usize) {
        self.nodes.reserve(additional)
    }

    pub fn node_indices(&self) -> impl Iterator<Item=Ix> + '_ {
        (0..self.nodes.len()).map(|ix| Ix::new(ix))
    }

    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    #[inline(always)]
    pub fn get_node(&self, ix: Ix) -> &StateTreeNode<N, O, Ix> {
        &self.nodes[ix.index()]
    }

    pub fn add_node(&mut self, node: StateTreeNode<N, O, Ix>) -> Ix {
        self.nodes.push(node);

        Ix::new(self.nodes.len() - 1)
    }
}

#[derive(Debug)]
pub struct StateTreeNode<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    graph_node: N,
    offset: O,
    state: AlignState,
    backtrace: Option<Backtrace<Ix>>,
}

impl<N, O, Ix> StateTreeNode<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    pub fn new_start(graph_node: N) -> Self {
        Self { graph_node, offset: O::zero(), state: AlignState::Start, backtrace: None }
    }

    pub fn new(graph_node: N, offset: O, state: AlignState, backtrace: Backtrace<Ix>) -> Self {
        Self { graph_node, offset, state, backtrace: Some(backtrace) }
    }

    pub fn node(&self) -> N {
        self.graph_node
    }

    pub fn offset(&self) -> O {
        self.offset
    }

    pub fn state(&self) -> AlignState {
        self.state
    }

    pub fn backtrace(&self) -> Option<&Backtrace<Ix>> {
        self.backtrace.as_ref()
    }
}

#[derive(Debug)]
pub enum Backtrace<Ix>
where
    Ix: TreeIndexType
{
    /// Represents a single alignment step. The only data stored in this
    /// variant is the index to the previous node in the alignment state tree.
    Step(Ix),

    /// Represents closing an indel, with the only data stored is the index of the align tree
    /// state node that represents the still open indel. It's a special backtrace state because it
    /// doesn't output any alignment characters.
    ClosedIndel(Ix),
}

impl<Ix> Backtrace<Ix>
where
    Ix: TreeIndexType,
{
    pub fn prev(&self) -> Ix {
        match *self {
            Self::Step(prev) => prev,
            Self::ClosedIndel(prev) => prev,
        }
    }
}
