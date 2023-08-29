use std::cell::{RefCell, RefMut};
use std::marker::PhantomData;
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentStateTree;
use crate::aligner::state::{AlignState, Backtrace, StateTreeNode, TreeIndexType};

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
pub struct StackNode<'a, N, O>(N, O, RefCell<Box<dyn Iterator<Item=N> + 'a>>);

impl<'a, N, O> StackNode<'a, N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    #[inline(always)]
    pub fn node(&self) -> N {
        self.0
    }

    #[inline(always)]
    pub fn offset(&self) -> O {
        self.1
    }

    #[inline]
    fn children_iter(&self) -> RefMut<Box<dyn Iterator<Item=N> + 'a>> {
        self.2.borrow_mut()
    }
}

/// A struct that represents the start node, end node, and length of an extended path
#[derive(Debug)]
pub enum ExtendedPath<N, O, Ix> {
    Matches(N, O, usize, Ix),
    Mismatch(N, O)
}

pub enum ExtendedChild<N, O> {
    None,
    Match(N, O),
    Mismatch(N, O),
}

pub trait VisitedProxy<N, O> {
    fn visited(&self, node: N, offset: O) -> bool;
    fn mark_visited(&mut self, node: N, offset: O);
}


/// A struct representing the depth-first extension state when
/// extending from a (node, query offset).
///
/// Extension occurs when a) the node symbol and the query symbol matches, and b) the path to the
/// current node brings the alignment further (i.e., the query offset is higher than before).
///
/// The algorithm is inspired by NetworkX's `dfs_labeled_edges` function, which in turn
/// is based on David Eppstein's depth-first search tree algorithm.
///
/// See also: https://www.ics.uci.edu/~eppstein/PADS/ or
/// https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.traversal.depth_first_search.dfs_labeled_edges.html
pub struct DepthFirstExtension<'a, G, N, O, V> {
    graph: &'a G,
    seq: &'a [u8],
    visited: V,
    stack: Vec<StackNode<'a, N, O>>
}

impl<'a, G, N, O, V> DepthFirstExtension<'a, G, N, O, V>
where
    G: AlignableGraph<NodeIndex=N>,
    N: NodeIndexType,
    O: OffsetType,
    V: VisitedProxy<N, O>,
{
    pub fn new(graph: &'a G, seq: &'a [u8], visited: V, start_node: N, start_offset: O) -> Self {
        let succ_iter = RefCell::new(Box::new(graph.successors(start_node)));
        Self {
            graph,
            seq,
            visited,
            stack: vec![StackNode(start_node, start_offset, succ_iter)]
        }
    }

    /// Find the next valid child of a given node
    ///
    /// A child is valid if it hasn't been visited before and we flag if the symbol in the query
    /// sequence matches the symbol of the node.
    fn next_valid_child(&self, parent: &StackNode<N, O>) -> ExtendedChild<N, O> {
        if parent.offset().as_usize() >= self.seq.len() {
            return ExtendedChild::None
        }

        while let Some(child) = parent.children_iter().next() {
            let child_offset = parent.offset().increase_one();

            let visited = self.visited.visited(child, child_offset);
            if !visited {
                return if self.graph.is_symbol_equal(child, self.seq[child_offset.as_usize() - 1]) {
                    ExtendedChild::Match(child, child_offset)
                } else {
                    ExtendedChild::Mismatch(child, child_offset)
                }
            }
        }

        ExtendedChild::None
    }

    pub fn get_visited_proxy(&self) -> &V {
        &self.visited
    }

    pub fn get_visited_proxy_mut(&mut self) -> &mut V {
        &mut self.visited
    }
}

impl<'a, G, N, O, V> Iterator for DepthFirstExtension<'a, G, N, O, V>
where
    G: AlignableGraph<NodeIndex=N>,
    N: NodeIndexType,
    O: OffsetType,
    V: VisitedProxy<N, O>,
{
    type Item = DFEEdge<N, O>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(parent) = self.stack.last() {
            let parent_node = parent.node();
            let parent_offset = parent.offset();

            return match self.next_valid_child(parent) {
                ExtendedChild::Match(child, child_offset) => {
                    self.visited.mark_visited(child, child_offset);

                    let child_succ = RefCell::new(Box::new(self.graph.successors(child)));
                    self.stack.push(StackNode(child, child_offset, child_succ));

                    Some(DFEEdge::ForwardMatch(parent_node, parent_offset, child, child_offset))
                },
                ExtendedChild::Mismatch(child, child_offset) => {
                    self.visited.mark_visited(child, child_offset);

                    Some(DFEEdge::ForwardMismatch(parent_node, parent_offset, child, child_offset))
                },
                ExtendedChild::None => {
                    let popped = self.stack.pop();

                    if let Some(new_top) = self.stack.last() {
                        popped.map(
                            |v| DFEEdge::Backward(v.node(), v.offset(), new_top.node(), new_top.offset())
                        )
                    } else {
                        None
                    }
                }
            }
        }

        None
    }
}


pub enum DFEEdge<N, O> {
    ForwardMatch(N, O, N, O),
    ForwardMismatch(N, O, N, O),
    Backward(N, O, N, O),
}

pub struct AlignStateTreeVisited<'a, T, N, O, Ix> {
    tree: &'a mut T,
    dummy: PhantomData<(N, O, Ix)>
}

impl<'a, T, N, O, Ix> AlignStateTreeVisited<'a, T, N, O, Ix>
where
    T: AlignmentStateTree<N, O, Ix>,
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    fn new(tree: &'a mut T) -> Self {
        Self {
            tree,
            dummy: PhantomData
        }
    }

    fn add_extend_node(&mut self, prev_state: Ix, node: N, offset: O) -> Ix
    where
        N: NodeIndexType,
        O: OffsetType,
    {
        let new_tree_node = StateTreeNode::new(
            node, offset,
            AlignState::Match,
            Backtrace::Step(prev_state)
        );

        self.tree.add_node(new_tree_node)
    }
}

impl<'a, T, N, O, Ix> VisitedProxy<N, O> for AlignStateTreeVisited<'a, T, N, O, Ix>
where
    T: AlignmentStateTree<N, O, Ix>,
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    fn visited(&self, node: N, offset: O) -> bool {
        self.tree.visited(node, offset, AlignState::Match)
    }

    fn mark_visited(&mut self, node: N, offset: O) {
        self.tree.mark_visited(node, offset, AlignState::Match);
    }
}


pub struct EndpointExtender<'a, G, N, O, T, Ix> {
    dfe: DepthFirstExtension<'a, G, N, O, AlignStateTreeVisited<'a, T, N, O, Ix>>,
    stack: Vec<Ix>,
    has_new_path: bool,
}

impl<'a, G, N, O, T, Ix> EndpointExtender<'a, G, N, O, T, Ix>
where
    G: AlignableGraph<NodeIndex=N>,
    N: NodeIndexType,
    O: OffsetType,
    T: AlignmentStateTree<N, O, Ix>,
    Ix: TreeIndexType,
{
    pub fn new(graph: &'a G, seq: &'a [u8], tree: &'a mut T, start_state: Ix) -> Self {
        let start_node = tree.get_node(start_state).node();
        let start_offset = tree.get_node(start_state).offset();
        let visited = AlignStateTreeVisited::new(tree);

        Self {
            dfe: DepthFirstExtension::new(graph, seq, visited, start_node, start_offset),
            stack: vec![start_state],
            has_new_path: false
        }
    }
}

pub struct NewExtendEndpoint<N, O, Ix>(AlignState, N, O, Ix, O);

impl<N, O, Ix> NewExtendEndpoint<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    pub fn state(&self) -> AlignState {
        self.0
    }

    pub fn node(&self) -> N {
        self.1
    }

    pub fn offset(&self) -> O {
        self.2
    }

    pub fn state_ix(&self) -> Ix {
        self.3
    }

    pub fn path_len(&self) -> O {
        self.4
    }
}

impl<'a, G, N, O, T, Ix> Iterator for EndpointExtender<'a, G, N, O, T, Ix>
where
    G: AlignableGraph<NodeIndex=N>,
    N: NodeIndexType,
    O: OffsetType,
    T: AlignmentStateTree<N, O, Ix>,
    Ix: TreeIndexType,
{
    type Item = NewExtendEndpoint<N, O, Ix>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(edge) = self.dfe.next() {
            match edge {
                DFEEdge::ForwardMatch(_, _, child, child_offset) => {
                    self.has_new_path = true;

                    let new_ix = self.dfe.get_visited_proxy_mut()
                        .add_extend_node(*self.stack.last().unwrap(), child, child_offset);
                    self.stack.push(new_ix)
                },
                DFEEdge::ForwardMismatch(_, _, child, child_offset) => {
                    let new_ix = self.dfe.get_visited_proxy_mut()
                        .add_extend_node(*self.stack.last().unwrap(), child, child_offset);
                    return Some(NewExtendEndpoint(
                        AlignState::Mismatch,
                        child, child_offset,
                        new_ix,
                        O::new(self.stack.len()))
                    );
                },
                DFEEdge::Backward(leaf_node, leaf_offset, _, _) => {
                    let leaf_state = self.stack.pop();

                    if self.has_new_path {
                        self.has_new_path = false;
                        return leaf_state.map(|v| NewExtendEndpoint(
                            AlignState::Match,
                            leaf_node, leaf_offset,
                            v,
                            O::new(self.stack.len() - 1))
                        )
                    }
                }
            }
        }

        None
    }
}
