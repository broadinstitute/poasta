use std::cell::{RefCell, RefMut};
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentStateTree;
use crate::aligner::state::{AlignState, TreeIndexType};

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
struct ExtendCandidate<'a, N, O>(N, O, RefCell<Box<dyn Iterator<Item=N> + 'a>>)
where
    N: NodeIndexType,
    O: OffsetType;

impl<'a, N, O> ExtendCandidate<'a, N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn offset(&self) -> O {
        self.1
    }

    pub fn children_iter(&self) -> RefMut<Box<dyn Iterator<Item=N> + 'a>> {
        self.2.borrow_mut()
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
pub struct PathExtender<'a, G, N, O>
where
    G: AlignableGraph<NodeIndex=N>,
    N: NodeIndexType,
    O: OffsetType,
{
    graph: &'a G,
    seq: &'a [u8],
    stack: Vec<ExtendCandidate<'a, N, O>>,
    curr_path: Vec<(N, O)>,
    has_new_path: bool,
}

impl<'a, G, N, O> PathExtender<'a, G, N, O>
where
    G: AlignableGraph<NodeIndex=N>,
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new(graph: &'a G, seq: &'a [u8], start_node: N, offset: O) -> Self {
        let succ_iter = RefCell::new(Box::new(graph.successors(start_node)));

        Self {
            graph,
            seq,
            stack: vec![ExtendCandidate(start_node, offset, succ_iter)],
            curr_path: vec![],
            has_new_path: false,
        }
    }

    /// Find the next valid child of a given node
    ///
    /// Valid children are those that have a matching symbol in the query and that 
    /// have not been visited before.
    fn next_valid_child<T, Ix>(&self, parent: &ExtendCandidate<N, O>, tree: &T) -> Option<(N, O)>
    where
        T: AlignmentStateTree<N, O, Ix>,
        Ix: TreeIndexType
    {
        if parent.offset().as_usize() >= self.seq.len() {
            return None
        }

        while let Some(child) = parent.children_iter().next() {
            let child_offset = parent.offset().increase_one();

            if self.graph.is_symbol_equal(child, self.seq[child_offset.as_usize() - 1]) &&
                !tree.visited(child, child_offset, AlignState::Match)
            {
                return Some((child, child_offset))
            }
        }

        None
    }

    pub fn next<T, Ix>(&mut self, tree: &mut T) -> Option<Vec<(N, O)>>
    where
        T: AlignmentStateTree<N, O, Ix>,
        Ix: TreeIndexType
    {
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();
            if let Some((child, child_offset)) = self.next_valid_child(parent, tree) {
                // Going forward in the depth-first matches search tree, set flag that we have a new path
                self.has_new_path = true;
                tree.mark_visited(child, child_offset, AlignState::Match);

                self.curr_path.push((child, child_offset));
                let child_succ = RefCell::new(Box::new(self.graph.successors(child)));
                self.stack.push(ExtendCandidate(child, child_offset, child_succ));
            } else {
                // If we reach here, we are about to move up the depth-first matches search tree,
                // but first, if we are at a leaf node, return the path and prepare for further
                // exploration of the graph
                let path_to_return = if self.has_new_path {
                    Some(self.curr_path.clone())
                } else {
                    None
                };

                // Exhausted children of current node, remove from stack and `curr_path`
                // to prepare for the next iteration
                self.stack.pop();
                self.curr_path.pop();

                if path_to_return.is_some() {
                    self.has_new_path = false;
                    return path_to_return
                }
            }
        }

        None
    }
}
