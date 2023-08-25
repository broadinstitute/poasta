use std::cell::{RefCell, RefMut};
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentStateTree;
use crate::aligner::state::{AlignState, Backtrace, StateTreeNode, TreeIndexType};

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
pub struct ExtendNode<'a, N, O, Ix>(N, O, Ix, RefCell<Box<dyn Iterator<Item=N> + 'a>>);

impl<'a, N, O, Ix> ExtendNode<'a, N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    #[inline(always)]
    pub fn node(&self) -> N {
        self.0
    }

    #[inline(always)]
    pub fn offset(&self) -> O {
        self.1
    }

    #[inline(always)]
    pub fn tree_ix(&self) -> Ix {
        self.2
    }

    #[inline]
    fn children_iter(&self) -> RefMut<Box<dyn Iterator<Item=N> + 'a>> {
        self.3.borrow_mut()
    }
}

/// A struct that represents the start node, end node, and length of an extended path
#[derive(Debug)]
pub struct ExtendedPath<N, O, Ix>(Ix, Ix, N, O, O);

impl<N, O, Ix> ExtendedPath<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    #[inline(always)]
    pub fn start(&self) -> Ix {
        self.0
    }

    #[inline(always)]
    pub fn end(&self) -> Ix {
        self.1
    }

    #[inline(always)]
    pub fn start_node(&self) -> N {
        self.2
    }

    #[inline(always)]
    pub fn start_offset(&self) -> O {
        self.3
    }

    #[inline(always)]
    pub fn len(&self) -> O {
        self.4
    }
}

/// A struct representing the depth-first search state when
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
pub struct PathExtender<'a, G, T, N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
{
    graph: &'a G,
    seq: &'a [u8],
    tree: &'a mut T,
    stack: Vec<ExtendNode<'a, N, O, Ix>>,
    has_new_path: bool,
    new_path_start_ix: usize,
}

impl<'a, G, T, N, O, Ix> PathExtender<'a, G, T, N, O, Ix>
where
    G: AlignableGraph<NodeIndex=N>,
    T: AlignmentStateTree<N, O, Ix>,
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    pub fn new(graph: &'a G, seq: &'a [u8], tree: &'a mut T, start_ix: Ix) -> Self {
        let start_node = tree.get_node(start_ix).node();
        let offset = tree.get_node(start_ix).offset();
        let succ_iter = RefCell::new(Box::new(graph.successors(start_node)));

        Self {
            graph,
            seq,
            tree,
            stack: vec![ExtendNode(start_node, offset, start_ix, succ_iter)],
            has_new_path: false,
            new_path_start_ix: 0,
        }
    }

    /// Find the next valid child of a given node
    ///
    /// Valid children are those that have a matching symbol in the query and that
    /// have not been visited before.
    fn next_valid_child(&self, parent: &ExtendNode<N, O, Ix>) -> Option<(N, O)> {
        if parent.offset().as_usize() >= self.seq.len() {
            return None
        }

        while let Some(child) = parent.children_iter().next() {
            let child_offset = parent.offset().increase_one();

            if self.graph.is_symbol_equal(child, self.seq[child_offset.as_usize() - 1]) &&
                !self.tree.visited(child, child_offset, AlignState::Match)
            {
                return Some((child, child_offset))
            }
        }

        None
    }
}

impl<'a, G, T, N, O, Ix> Iterator for PathExtender<'a, G, T, N, O, Ix>
where
    G: AlignableGraph<NodeIndex=N>,
    T: AlignmentStateTree<N, O, Ix>,
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    type Item = ExtendedPath<N, O, Ix>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();
            if let Some((child, child_offset)) = self.next_valid_child(parent) {
                // We are going forward in the depth-first matches search tree, so set flag that
                // we have a new path
                if !self.has_new_path {
                    self.new_path_start_ix = self.stack.len() - 1;
                }

                self.has_new_path = true;

                // Add the reached (node, offset) to the state tree.
                let new_tree_ix = self.tree.add_node(
                    StateTreeNode::new(child, child_offset, AlignState::Match,
                                       Backtrace::Step(self.stack.last().unwrap().tree_ix()))
                );

                self.tree.mark_visited(child, child_offset, AlignState::Match);

                let child_succ = RefCell::new(Box::new(self.graph.successors(child)));
                self.stack.push(ExtendNode(child, child_offset, new_tree_ix, child_succ));
            } else {
                // if self.has_new_path {
                //     let path: Vec<_> = self.stack[self.new_path_start_ix..].iter()
                //         .map(|item| (item.0, item.1))
                //         .collect();
                //
                //     eprintln!("Path: {path:?}");
                // }

                // If we reach here, we are about to move up the depth-first matches search tree.
                // Pop item on top of the stack, and see if we can navigate to other children of nodes
                // on the stack. If at a leaf node (`has_new_path` = true), then return the tree
                // indices of the path.
                let popped = self.stack.pop();

                if self.has_new_path {
                    self.has_new_path = false;
                    return popped.map(|v| {
                        ExtendedPath(
                            self.stack[self.new_path_start_ix].tree_ix(),
                            v.tree_ix(),
                            self.stack[self.new_path_start_ix].node(),
                            self.stack[self.new_path_start_ix].offset(),
                            O::new(self.stack.len() - self.new_path_start_ix + 1)
                        )
                    });
                }
            }
        }

        None
    }
}
