use std::cell::{RefCell, RefMut};
use std::collections::HashMap;
use crate::aligner::layers::Layer;
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::{Diag, DiagonalPoint, MatrixPoint, OffsetType};
use crate::aligner::visited::{AlignState, VisitedIntervalType};
use crate::aligner::visited::{Backtrace, VisitedInterval, VisitedIntervalData, VisitedIntervalTree};

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
pub(crate) struct ExtendNode<'a, N, O>(N, O, RefCell<Box<dyn Iterator<Item=N> + 'a>>)
where
    N: NodeIndexType,
    O: OffsetType;

impl<'a, N, O> ExtendNode<'a, N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn node(&self) -> N {
        self.0
    }

    pub fn offset(&self) -> O {
        self.1
    }

    pub fn diag<G: AlignableGraph<NodeIndex=N>>(&self, graph: &G) -> Diag {
        (graph.node_ix_to_row(self.0), self.1).diag()
    }

    fn children_iter(&self) -> RefMut<Box<dyn Iterator<Item=N> + 'a>> {
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
pub struct PathExtender<'a, G, O>
where
    O: OffsetType,
{
    graph: &'a G,
    seq: &'a [u8],
    layer: &'a Layer<O>,
    new_intervals: HashMap<Diag, VisitedIntervalTree<O>>,
}

impl<'a, G, O> PathExtender<'a, G, O>
where
    G: AlignableGraph,
    O: OffsetType,
{
    pub(crate) fn new(graph: &'a G, seq: &'a [u8], layer: &'a Layer<O>) -> Self {
        Self {
            graph,
            seq,
            layer,
            new_intervals: HashMap::new(),
        }
    }

    fn next_interval_start_on_diag(&self, diag: Diag, curr_offset: O) -> Option<O> {
        [
            self.new_intervals.get(&diag)
                .and_then(|ivals| ivals.find(curr_offset, O::max_value())
                    .next()
                    .map(|ival| ival.start)),
            self.layer.get_visited(diag)
                .and_then(|ivals| ivals.find(curr_offset, O::max_value())
                    .next()
                    .map(|ival| ival.start)),
        ].into_iter()
            .flatten()
            .min()
    }

    pub fn extend(&mut self, start_diag: Diag, start_ival: &'a VisitedIntervalType<O>) {
        let row: G::NodeIndex = (start_diag, start_ival.end() - O::one()).row();
        let start_node = self.graph.row_to_node_ix(row);
        eprintln!("Start ival: {:?}, node: {:?} (row: {:?})", start_ival, start_node, row);

        let succ_iter = RefCell::new(Box::new(self.graph.successors(start_node)));
        let mut node_stack = vec![ExtendNode(start_node, start_ival.end() - O::one(), succ_iter)];
        let start_interval = start_ival.clone();

        let mut ival_stack = vec![(start_diag, start_interval)];
        let mut prev_diag = start_diag;
        let mut next_ival_start = self.next_interval_start_on_diag(start_diag, start_ival.end());
        let mut has_new_path = false;

        while !node_stack.is_empty() {
            let parent = node_stack.last().unwrap();
            let parent_node = parent.node();
            let parent_diag = parent.diag(self.graph);

            eprintln!("Curr diag: {:?}, offset: {:?}, row: {:?}, next interval on diag: {:?}", parent_diag, parent.offset(),
                      self.graph.node_ix_to_row(parent.node()), next_ival_start);
            eprintln!("Ival stack: {:?}", ival_stack);

            // See if we can find a visited interval on the same diagonal further downstream
            // of our current position, which aids in determining how far we can extend along
            // this diagonal.
            if parent_diag != prev_diag {
                next_ival_start = self.next_interval_start_on_diag(parent_diag, parent.offset().increase_one());
            }

            if let Some((child, child_diag, child_offset)) = self.next_valid_child(parent, next_ival_start) {
                eprintln!("Forward {:?}, offset: {:?}, row: {:?}, {}-{}",
                          child_diag, child_offset, self.graph.node_ix_to_row(child),
                          self.graph.get_symbol(child), char::from(self.seq[child_offset.as_usize()-1]));

                // Going forward in the depth-first matches search tree, set flag that we have a new path
                has_new_path = true;

                let child_succ = RefCell::new(Box::new(self.graph.successors(child)));
                node_stack.push(ExtendNode(child, child_offset, child_succ));

                let last_interval = &mut ival_stack.last_mut().unwrap().1;
                if child_diag == parent_diag {
                    // When staying on the same diagonal, we can simply increase interval length
                    *last_interval.end_mut() = last_interval.end() + O::one();
                    last_interval.clipped_mut().saturating_sub(O::one());
                    last_interval.set_extendable(true);
                } else {
                    // Changing diagonals, add interval on previous diagonal to our newly
                    // created intervals hash map, and start a new one
                    last_interval.set_extendable(false);
                    eprintln!("Changing to diagonal {:?}, adding {:?} to newly added intervals (diag: {:?})", child_diag, last_interval, parent_diag);

                    self.new_intervals.entry(parent_diag)
                        .and_modify(|itree| itree.insert_and_index(last_interval.clone()))
                        .or_insert_with(|| {
                            VisitedIntervalTree::from_iter(vec![last_interval.clone()])
                        });

                    let new_interval = VisitedIntervalType::new_size_one(
                        child_offset, VisitedIntervalData::new(true, O::zero(),
                                                               Backtrace::new(parent_diag, AlignState::Extended)));
                    eprintln!("Starting new interval: {:?}", new_interval);

                    ival_stack.push((child_diag, new_interval));
                }

                prev_diag = child_diag;
            } else {
                eprintln!("Up");
                // If we reach here, we are about to move up the depth-first matches search tree.
                // If we are at a leaf node (has_new_path = true), then add last interval to newly
                // added interval hashmap. Then, remove node from stack, and also update the interval stack.
                if has_new_path {
                    let to_insert = &ival_stack.last().unwrap().1;
                    eprintln!("At leaf, adding {:?} to newly added intervals", to_insert);

                    self.new_intervals.entry(parent_diag)
                        .and_modify(|itree| itree.insert_and_index(to_insert.clone()))
                        .or_insert_with(|| {
                            VisitedIntervalTree::from_iter(vec![to_insert.clone()])
                        });

                    has_new_path = false;
                }

                node_stack.pop();

                let mut remove_top_interval = false;
                if let Some(top_of_stack) = ival_stack.last_mut() {
                    let last_interval = &mut top_of_stack.1;
                    *last_interval.end_mut() = last_interval.end() - O::one();
                    remove_top_interval = last_interval.start == last_interval.end;
                    eprintln!("Changed last interval: {:?}", last_interval);
                    prev_diag = top_of_stack.0;
                }

                if remove_top_interval {
                    eprintln!("Remove interval from stack");
                    ival_stack.pop();

                    if let Some(top_of_stack) = ival_stack.last() {
                        prev_diag = top_of_stack.0;
                    }
                }
            }
        }
    }

    /// Check both our newly added intervals and the existing intervals in a layer to determine
    /// whether a (diag, offset) has been visited before.
    fn visited(&self, diag: Diag, offset: O) -> bool {
        self.new_intervals.get(&diag)
            .map_or(false, |ivals| ivals.contains(offset))
            || self.layer.visited(diag, offset)
    }

    /// Find the next valid child of a given node, or any child whose symbol matching
    /// that in the query and which hasn't been visited yet.
    ///
    /// If a child is on the same diagonal, we can take a shortcut to check if a (node, offset)
    /// combination is already visited. We use the pre-computed start location of the next visited
    /// interval on the same diagonal, and check if the child offset is lower than that. This saves
    /// a number of queries to the interval trees, that have logarithmic time complexity to check
    /// if a offset is already visited.
    ///
    /// # Arguments
    /// * parent - A reference to an [`ExtendedNode`], which holds the node to explore for
    ///   successor nodes, and the current query offset.
    /// * next_ival_start - The query offset of the next visited interval the diagonal of `parent`,
    ///   if any.
    fn next_valid_child<N: NodeIndexType>(&self, parent: &ExtendNode<N, O>, next_ival_start: Option<O>) -> Option<(N, Diag, O)>
    where
        G: AlignableGraph<NodeIndex=N>
    {
        if parent.offset().as_usize() >= self.seq.len() {
            return None
        }

        while let Some(child) = parent.children_iter().next() {
            let child_offset = parent.offset().increase_one();

            if self.graph.is_symbol_equal(child, self.seq[child_offset.as_usize() - 1]) {
                let child_row = self.graph.node_ix_to_row(child);
                let child_diag = (child_row, child_offset).diag();

                if child_diag == parent.diag(self.graph) {
                    if next_ival_start.map_or(true, |max_offset| child_offset < max_offset) {
                        return Some((child, child_diag, child_offset));
                    }
                } else if !self.visited(child_diag, child_offset) {
                    return Some((child, child_diag, child_offset));
                }
            }
        }

        None
    }

    pub fn get_new_intervals(&self) -> &HashMap<Diag, VisitedIntervalTree<O>> {
        &self.new_intervals
    }

    pub fn into_new_intervals(self) -> HashMap<Diag, VisitedIntervalTree<O>> {
        self.new_intervals
    }
}
