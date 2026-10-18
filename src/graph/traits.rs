use std::fmt;
use std::hash::Hash;

pub trait GraphNodeId: fmt::Debug + Clone + Copy + PartialEq + Eq + Hash + Ord {
    fn index(&self) -> usize;
}

/// The integer type a graph stores node and edge indices in. `max()` is
/// reserved as "none".
pub trait IndexType: Copy + Default + Hash + Ord + fmt::Debug + 'static {
    fn new(x: usize) -> Self;
    fn index(&self) -> usize;
    fn max() -> Self;
}

macro_rules! index_type {
    ($($t:ty),*) => {$(
        impl IndexType for $t {
            #[inline(always)]
            fn new(x: usize) -> Self {
                x as $t
            }
            #[inline(always)]
            fn index(&self) -> usize {
                *self as usize
            }
            #[inline(always)]
            fn max() -> Self {
                <$t>::MAX
            }
        }
    )*};
}

index_type!(u8, u16, u32, u64, usize);

impl<T> GraphNodeId for T
where
    T: IndexType,
{
    #[inline(always)]
    fn index(&self) -> usize {
        IndexType::index(self)
    }
}

pub trait GraphBase {
    type NodeType: GraphNodeId;

    type NodeIter<'a>: Iterator<Item = Self::NodeType> + 'a
    where
        Self: 'a;
    type Successors<'a>: Iterator<Item = Self::NodeType> + 'a
    where
        Self: 'a;
    type Predecessors<'a>: Iterator<Item = Self::NodeType> + 'a
    where
        Self: 'a;

    fn all_nodes_iter(&self) -> Self::NodeIter<'_>;
    fn node_count(&self) -> usize;
    fn is_empty(&self) -> bool;

    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_>;
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_>;

    fn out_degree(&self, node: Self::NodeType) -> usize;
    fn in_degree(&self, node: Self::NodeType) -> usize;
}

/// Trait for graphs that have dedicated start and end nodes, i.e.,
/// a start node without any incoming edges and an end node without any outgoing edges.
pub trait GraphWithNodeOrdering: GraphBase {
    fn start_node(&self) -> Self::NodeType;
    fn end_node(&self) -> Self::NodeType;

    fn node_rank(&self, node: Self::NodeType) -> usize;
    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType;
}
