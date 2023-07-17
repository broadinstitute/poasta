pub mod poa;

use std::fmt::Debug;
use std::hash::Hash;

use petgraph::graph::IndexType;

pub trait NodeIndexType: Copy + Hash + PartialOrd + Ord + PartialEq + Eq + Debug + Default { }
impl<T: IndexType> NodeIndexType for T { }

pub trait AlignableGraph {
    type NodeIndex: NodeIndexType;
    type SuccessorIterator<'a>: Iterator<Item=Self::NodeIndex> + 'a
        where Self: 'a;

    fn start_nodes(&self) -> &Vec<Self::NodeIndex>;
    fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_>;

    fn is_end(&self, node: Self::NodeIndex) -> bool;

    fn get_symbol(&self, node: Self::NodeIndex) -> char;
    fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool;
}
