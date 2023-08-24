pub mod poa;

use std::fmt::Debug;
use std::hash::Hash;

use petgraph::graph::IndexType;

pub trait NodeIndexType: Copy + Hash + PartialOrd + Ord + PartialEq + Eq + Debug + Default {
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
}

impl<T: IndexType> NodeIndexType for T {
    #[inline(always)]
    fn new(value: usize) -> Self {
        Self::new(value)
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        self.index()
    }
}

pub trait AlignableGraph {
    type NodeIndex: NodeIndexType;
    type NodeIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;

    type PredecessorIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;
    type SuccessorIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;

    fn all_nodes(&self) -> Self::NodeIterator<'_>;
    fn node_count(&self) -> usize;
    fn edge_count(&self) -> usize;

    fn start_nodes(&self) -> &Vec<Self::NodeIndex>;
    fn end_nodes(&self) -> &Vec<Self::NodeIndex>;
    fn node_ix_to_row(&self, node: Self::NodeIndex) -> Self::NodeIndex;
    fn row_to_node_ix<T: NodeIndexType>(&self, row: T) -> Self::NodeIndex;

    fn predecessors(&self, node: Self::NodeIndex) -> Self::PredecessorIterator<'_>;
    fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_>;

    fn is_successor(&self, row1: Self::NodeIndex, row2: Self::NodeIndex) -> bool;

    fn is_end(&self, node: Self::NodeIndex) -> bool;

    fn get_symbol(&self, node: Self::NodeIndex) -> char;
    fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool;
}
