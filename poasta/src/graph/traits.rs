use std::fmt;
use std::hash::Hash;

pub trait GraphNodeId: fmt::Debug + Clone + Copy + PartialEq + Eq + Hash {
    fn index(&self) -> usize;
}

pub trait GraphBase {
    type NodeType: GraphNodeId;
    
    type NodeIter<'a>: Iterator<Item=Self::NodeType> + 'a
        where Self: 'a;
    type Successors<'a>: Iterator<Item=Self::NodeType> + 'a
        where Self: 'a;
    type Predecessors<'a>: Iterator<Item=Self::NodeType> + 'a
        where Self: 'a;

    fn all_nodes_iter(&self) -> Self::NodeIter<'_>;
    fn node_count(&self) -> usize;
    
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_>;
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_>;
}

/// Trait for graphs that have dedicated start and end nodes, i.e.,
/// a start node without any incoming edges and an end node without any outgoing edges.
pub trait GraphWithStartEnd: GraphBase {
    fn start_node(&self) -> <Self as GraphBase>::NodeType;
    fn end_node(&self) -> <Self as GraphBase>::NodeType;
}
