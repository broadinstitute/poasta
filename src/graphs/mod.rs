pub mod poa;

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

pub trait AlignableGraph {
    type NodeIndex: NodeIndexType;
    type NodeIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;
    type SuccessorIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;

    fn all_nodes(&self) -> Self::NodeIterator<'_>;
    fn node_count(&self) -> usize;
    fn start_nodes(&self) -> &Vec<Self::NodeIndex>;
    fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_>;

    fn is_end(&self, node: Self::NodeIndex) -> bool;

    fn get_symbol(&self, node: Self::NodeIndex) -> char;
    fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool;
}
