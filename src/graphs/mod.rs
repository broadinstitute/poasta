pub mod poa;
pub mod tools;

#[cfg(test)]
pub(crate) mod mock;

use std::fmt::Debug;
use std::hash::Hash;

use petgraph::graph::IndexType;

pub trait NodeIndexType: Copy + Hash + PartialOrd + Ord + PartialEq + Eq + Debug + Default {
    fn index(&self) -> usize;
}

impl<T: IndexType> NodeIndexType for T {
    #[inline(always)]
    fn index(&self) -> usize {
        self.index()
    }
}

pub trait AlignableRefGraph {
    type NodeIndex: NodeIndexType;
    type NodeIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;
    type PredecessorIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;
    type SuccessorIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;

    fn all_nodes(&self) -> Self::NodeIterator<'_>;
    fn node_count(&self) -> usize;
    fn node_count_with_start_and_end(&self) -> usize;

    fn edge_count(&self) -> usize;

    fn start_node(&self) -> Self::NodeIndex;
    fn end_node(&self) -> Self::NodeIndex;

    fn predecessors(&self, node: Self::NodeIndex) -> Self::PredecessorIterator<'_>;
    fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_>;

    fn in_degree(&self, node: Self::NodeIndex) -> usize;
    fn out_degree(&self, node: Self::NodeIndex) -> usize;

    fn is_end(&self, node: Self::NodeIndex) -> bool;

    fn get_symbol(&self, node: Self::NodeIndex) -> char;
    fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool;

    fn get_node_ordering(&self) -> Vec<usize>;
}
