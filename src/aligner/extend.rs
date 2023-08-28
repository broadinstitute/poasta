use std::cell::{RefCell, RefMut};
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentStateTree};
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
pub struct ExtendedPath<N, O, Ix>(Ix, N, O, O);

impl<N, O, Ix> ExtendedPath<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    #[inline(always)]
    pub fn end(&self) -> Ix {
        self.0
    }

    #[inline(always)]
    pub fn end_node(&self) -> N {
        self.1
    }

    #[inline(always)]
    pub fn end_offset(&self) -> O {
        self.2
    }

    #[inline(always)]
    pub fn len(&self) -> O {
        self.3
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
    curr_score: usize,
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
    pub fn new(graph: &'a G, seq: &'a [u8], tree: &'a mut T, curr_score: usize, start_ix: Ix) -> Self {
        let start_node = tree.get_node(start_ix).node();
        let offset = tree.get_node(start_ix).offset();
        let succ_iter = RefCell::new(Box::new(graph.successors(start_node)));

        Self {
            graph,
            seq,
            tree,
            curr_score,
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
                !self.tree.visited(child, child_offset, AlignState::Match, self.curr_score)
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
    type Item = (ExtendNode<'a, N, O, Ix>, O);

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
                let new_tree_node = StateTreeNode::new(child, child_offset, AlignState::Match,
                                                       Backtrace::Step(self.stack.last().unwrap().tree_ix()));
                if self.stack.len() > 1 {
                    let last = self.stack.last().unwrap();
                    self.tree.add_intermediate_extend_node(self.curr_score, last.node(), last.offset(),  O::new(self.stack.len() - 1));
                }

                let new_tree_ix = self.tree.add_node(new_tree_node);
                self.tree.mark_visited(child, child_offset, AlignState::Match);

                let child_succ = RefCell::new(Box::new(self.graph.successors(child)));
                self.stack.push(ExtendNode(child, child_offset, new_tree_ix, child_succ));
            } else {
                // if self.has_new_path {
                //     let path: Vec<_> = self.stack[self.new_path_start_ix..].iter()
                //         .map(|item| (item.0, item.1))
                //         .collect();
                //
                //     eprintln!("Path: {path:?} (path len: {})", (self.stack.len() - self.new_path_start_ix));
                // }

                // If we reach here, we are about to move up the depth-first matches search tree.
                // Pop item on top of the stack, and see if we can navigate to other children of nodes
                // on the stack. If at a leaf node (`has_new_path` = true), then return the tree
                // indices of the path.
                let popped = self.stack.pop();

                if self.has_new_path {
                    self.has_new_path = false;
                    return popped.map(|v| (v, O::new(self.stack.len() - self.new_path_start_ix)));
                }
            }
        }

        None
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ExtendHitType {
    Full,
    DeletionExtension
}

#[derive(Copy, Clone, Debug)]
pub struct ExtendHit<O>(usize, O, O, ExtendHitType)
where
    O: OffsetType;

impl<O> ExtendHit<O>
where
    O: OffsetType,
{
    pub fn new(score: usize, offset: O, path_offset: O, extend_type: ExtendHitType) -> Self {
        Self(score, offset, path_offset, extend_type)
    }

    pub fn score(&self) -> usize {
        self.0
    }

    pub fn offset(&self) -> O {
        self.1
    }

    pub fn path_offset(&self) -> O {
        self.2
    }

    pub fn extend_type(&self) -> ExtendHitType {
        self.3
    }

    /// If node N is reached at offset x during the extension step, and it is the y'th node in that
    /// path, we can make some inferences whether other alignment states are reachable within
    /// a given score.
    ///
    /// # Example
    /// ```
    /// use poasta::aligner::offsets::OffsetType;
    /// use poasta::aligner::scoring::GapAffine;
    /// use poasta::aligner::extend::{ExtendHit, ExtendHitType};
    /// use poasta::aligner::state::AlignState;
    ///
    /// let costs = GapAffine::new(4, 2, 6);
    /// let hit = ExtendHit::<u32>::new(12, 20, 5, ExtendHitType::Full);
    ///
    /// // It's own offset is reachable
    /// assert_eq!(hit.reachable_from(costs, 20, 20), true);
    ///
    /// // Cases below are examples of deletions opened from previous nodes in the path
    /// assert_eq!(hit.reachable_from(costs, 19, 14, AlignState::Match), false); // too low score
    /// assert_eq!(hit.reachable_from(costs, 19, 20, AlignState::Match), true);
    /// assert_eq!(hit.reachable_from(costs, 18, 20, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 18, 22, AlignState::Match), true);
    /// assert_eq!(hit.reachable_from(costs, 17, 22, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 17, 24, AlignState::Match), true);
    /// assert_eq!(hit.reachable_from(costs, 16, 24, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 16, 26, AlignState::Match), true);
    /// assert_eq!(hit.reachable_from(costs, 15, 26, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 15, 28, AlignState::Match), true);
    ///
    /// // Going beyond the extended path should all return false
    /// assert_eq!(hit.reachable_from(costs, 14, 30, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 14, 32, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 12, 40, AlignState::Match), false);
    ///
    /// // Examples where we would open an insertion from the node reached by extension
    /// assert_eq!(hit.reachable_from(costs, 21, 14, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 21, 20, AlignState::Match), true);
    /// assert_eq!(hit.reachable_from(costs, 22, 20, AlignState::Match), false);
    /// assert_eq!(hit.reachable_from(costs, 22, 22, AlignState::Match), true);
    ///
    /// // Deletion extensions only check for deletions from previous nodes, not insertions
    /// let hit2 = ExtendHit::<u32>::new(12, 20, 5, ExtendHitType::DeletionExtension);
    /// assert_eq!(hit2.reachable_from(costs, 19, 14, AlignState::Match), false); // too low score
    /// assert_eq!(hit2.reachable_from(costs, 19, 20, AlignState::Match), true);
    /// assert_eq!(hit2.reachable_from(costs, 18, 20, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 18, 22, AlignState::Match), true);
    /// assert_eq!(hit2.reachable_from(costs, 17, 22, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 17, 24, AlignState::Match), true);
    /// assert_eq!(hit2.reachable_from(costs, 16, 24, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 16, 26, AlignState::Match), true);
    /// assert_eq!(hit2.reachable_from(costs, 15, 26, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 15, 28, AlignState::Match), true);
    ///
    /// assert_eq!(hit2.reachable_from(costs, 21, 14, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 21, 20, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 22, 20, AlignState::Match), false);
    /// assert_eq!(hit2.reachable_from(costs, 22, 22, AlignState::Match), false);
    ///
    /// ```
    pub fn reachable_from<C: AlignmentCosts>(&self, costs: C, offset: O, score: usize, state: AlignState) -> bool {
        // eprintln!("{self:?}");
        if offset == self.offset() {
            return true;
        }

        if (state == AlignState::Insertion || state == AlignState::Insertion2) && offset < self.offset() {
            return false;
        }

        if (state == AlignState::Deletion || state == AlignState::Deletion2) && offset > self.offset() {
            return false;
        }

        if offset < self.offset() - self.path_offset() {
            // eprintln!("outside deletion range.");
            return false;
        }

        if self.extend_type() == ExtendHitType::DeletionExtension && offset >= self.offset() {
            // eprintln!("Type = deletion, and offset indicates insertion.");
            return false;
        }

        let gap_length = (offset.as_isize() - self.offset().as_isize()).unsigned_abs();
        // eprintln!("gap length from hit to offset: {:?}-{:?}: {gap_length}", self.offset(), offset);

        let score_from_hit = self.score() + costs.gap_score(gap_length);
        // eprintln!("Score from hit: {score_from_hit} <= {score}: {}", score_from_hit <= score);
        score_from_hit <= score
    }

    pub fn to_next(&self) -> Self {
        Self(
            self.score(),
            self.offset(),
            self.path_offset(),
            ExtendHitType::DeletionExtension
        )
    }

}
