use std::cell::{RefCell, RefMut};
use crate::wavefront::offsets::OffsetPrimitive;
use crate::graph::POAGraph;
use crate::wavefront::compute::WFCompute;

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
struct ExtendCandidate<'a>(usize, usize, RefCell<Box<dyn Iterator<Item=usize> + 'a>>);

impl<'a> ExtendCandidate<'a> {
    pub fn node(&self) -> usize {
        self.0
    }

    pub fn offset(&self) -> usize {
        self.1
    }

    pub fn children_iter(&self) -> RefMut<Box<dyn Iterator<Item=usize> + 'a>> {
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
pub struct ExtendPaths<'a> {
    graph: &'a POAGraph,
    seq: &'a [u8],
    stack: Vec<ExtendCandidate<'a>>,
    curr_path: Vec<(usize, usize)>,
    has_new_path: bool
}

impl<'a> ExtendPaths<'a> {
    pub fn new<Offset: OffsetPrimitive>(graph: &'a POAGraph, seq: &'a [u8],
               start_node: usize, offset: Offset) -> Self {
        let offset_as_usize: usize = match offset.try_into() {
            Ok(v) => v,
            Err(_) => panic!("Could not obtain offset!")
        };

        let succ_iter = RefCell::new(graph.successors(start_node));

        Self {
            graph, seq,
            stack: vec![ExtendCandidate(start_node, offset_as_usize, succ_iter)],
            curr_path: vec![],
            has_new_path: false
        }
    }

    /// Find the next valid child of a given node
    ///
    /// Valid children are those that have a matching symbol in the query, and bring the alignment
    /// further (i.e., the query offset is equal or higher than the current highest value)
    fn next_valid_child<Compute: WFCompute>(&self, parent: &ExtendCandidate, compute: &Compute) -> Option<(usize, usize)> {
        if parent.offset() >= self.seq.len() - 1 {
            return None
        }

        if let Some(child) = parent.children_iter().next() {
            let child_offset = parent.offset() + 1;

            eprintln!("- checking rank {} vs child offset {} - 1 = {}, equal: {}", child, child_offset, child_offset-1,
                      self.graph.is_symbol_equal(child, self.seq[child_offset-1]));
            if self.graph.is_symbol_equal(child, self.seq[child_offset-1]) && compute.is_further(child, child_offset) {
                eprintln!("- is further too");
                Some((child, child_offset))
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn next<Compute: WFCompute>(&mut self, compute: &Compute) -> Option<Vec<(usize, usize)>> {
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();
            if let Some((child, child_offset)) = self.next_valid_child(parent, compute) {
                // Going forward in the depth-first matches search tree, set flag that we have a new path
                self.has_new_path = true;

                self.curr_path.push((child, child_offset));
                let child_succ = RefCell::new(self.graph.successors(child));
                self.stack.push(ExtendCandidate(child, child_offset, child_succ));
            } else {
                // If we reach here, we are about to move up the depth-first matches search tree,
                // if we are at a leaf node, return the path and prepare for further exploration of the graph
                let path_to_return = if self.has_new_path {
                    eprintln!("At leaf node, path: {:?}", self.curr_path);
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
