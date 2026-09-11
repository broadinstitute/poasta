//! Depth-first greedy match extension.
//!
//! Match edges cost 0, so once a `Match` state pops off the queue the search
//! walks matching `(node, offset+1)` successors depth first for free and
//! re-enters the queue only at a mismatch, the graph end or the query end.
//! This is the biggest constant-factor win in the engine.

use std::cell::{RefCell, RefMut};

use crate::align::AlignState;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::astar::visited::AstarVisited;
use crate::align::offsets::OffsetType;
use crate::align::traits::AlignableGraph;
use crate::graph::traits::GraphNodeId;

/// A node in the reference graph aligned to a certain query offset that we
/// will try to extend from, together with its (not-yet-exhausted) successor
/// iterator.
struct StackNode<N, O, I>(AlignmentGraphNode<N, O>, RefCell<I>)
where
    N: GraphNodeId,
    O: OffsetType;

impl<N, O, I> StackNode<N, O, I>
where
    N: GraphNodeId,
    O: OffsetType,
    I: Iterator,
{
    fn new(state: AlignmentGraphNode<N, O>, iter: I) -> Self {
        Self(state, RefCell::new(iter))
    }

    #[inline(always)]
    fn offset(&self) -> O {
        self.0.offset()
    }

    #[inline(always)]
    fn aln_graph_node(&self) -> &AlignmentGraphNode<N, O> {
        &self.0
    }

    #[inline]
    fn children_iter(&self) -> RefMut<'_, I> {
        self.1.borrow_mut()
    }
}

#[derive(Debug)]
pub enum ExtendResult<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    /// Reached the reference-graph end node. Holds the alignment-graph edge
    /// traversed: the parent node and the just-reached end-sentinel node.
    RefGraphEnd(AlignmentGraphNode<N, O>, AlignmentGraphNode<N, O>),

    /// Could not extend further because all query sequence has been aligned.
    /// Holds the parent alignment-graph node and the graph node we tried
    /// (but couldn't) align next.
    QueryEnd(AlignmentGraphNode<N, O>, N),

    /// Encountered a reference-graph node with a different symbol than the
    /// corresponding query symbol. Holds the last matching alignment-graph
    /// node (so callers can still open indel states from it) and the
    /// mismatching alignment-graph node.
    Mismatch(AlignmentGraphNode<N, O>, AlignmentGraphNode<N, O>),

    /// A just-matched state satisfies the caller-supplied `is_end` predicate
    /// (used by `EndsFree`/semi-global alignment, where a valid end state can
    /// sit in the middle of the graph, not just at the true end sentinel,
    /// so it must be checked on every match, not only when a state is popped
    /// off the priority queue).
    EndReached(AlignmentGraphNode<N, O>),
}

enum Successor<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    RefGraphEnd(AlignmentGraphNode<N, O>),
    QueryEnd(N),
    Match(AlignmentGraphNode<N, O>),
    Mismatch(AlignmentGraphNode<N, O>),
    SuccessorsExhausted,
}

pub struct DepthFirstGreedyAlignment<'a, G, O>
where
    G: AlignableGraph + 'a,
    O: OffsetType,
{
    ref_graph: &'a G,
    seq: &'a [u8],
    score: u32,
    num_visited: usize,
    stack: Vec<StackNode<G::Node, O, G::Successors<'a>>>,
}

impl<'a, G, O> DepthFirstGreedyAlignment<'a, G, O>
where
    G: AlignableGraph,
    O: OffsetType,
{
    pub fn new(ref_graph: &'a G, seq: &'a [u8], score: u32, start_node: &AlignmentGraphNode<G::Node, O>) -> Self {
        Self {
            ref_graph,
            seq,
            score,
            num_visited: 0,
            stack: vec![StackNode::new(*start_node, ref_graph.successors(start_node.node()))],
        }
    }

    pub fn get_num_visited(&self) -> usize {
        self.num_visited
    }

    /// `is_end` is checked against every freshly matched state before it is
    /// pushed for further extension. This is needed so `EndsFree` end states in the
    /// middle of the graph aren't silently walked past (see
    /// `ExtendResult::EndReached`). For `Global` mode this is essentially free:
    /// it's only ever true at the true end sentinel, which the DFA never
    /// reaches via a `Match` step in the first place.
    pub fn extend<V, F>(&mut self, astar_visited: &mut V, is_end: F) -> Option<ExtendResult<G::Node, O>>
    where
        V: AstarVisited<G::Node, O>,
        F: Fn(&AlignmentGraphNode<G::Node, O>) -> bool,
    {
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();

            match self.next_valid_successor(parent, astar_visited) {
                Successor::RefGraphEnd(aln_node) => {
                    return Some(ExtendResult::RefGraphEnd(*parent.aln_graph_node(), aln_node));
                }
                Successor::QueryEnd(ref_node) => {
                    return Some(ExtendResult::QueryEnd(*parent.aln_graph_node(), ref_node));
                }
                Successor::Match(child) => {
                    if astar_visited.prune(self.score, &child, AlignState::Match) {
                        continue;
                    }

                    astar_visited.dfa_match(self.score, parent.aln_graph_node(), &child);
                    self.num_visited += 1;
                    tracing::trace!(
                        target: "poasta::astar::extend",
                        score = self.score,
                        node = child.node().index(),
                        offset = parent.offset().as_usize(),
                        extended_offset = child.offset().as_usize(),
                        state = "Match",
                    );

                    if is_end(&child) {
                        return Some(ExtendResult::EndReached(child));
                    }

                    let child_succ = self.ref_graph.successors(child.node());
                    self.stack.push(StackNode::new(child, child_succ));
                }
                Successor::Mismatch(child) => {
                    return Some(ExtendResult::Mismatch(*parent.aln_graph_node(), child));
                }
                Successor::SuccessorsExhausted => {
                    self.stack.pop();
                }
            }
        }

        None
    }

    fn next_valid_successor<V>(
        &self,
        parent: &StackNode<G::Node, O, G::Successors<'a>>,
        astar_visited: &mut V,
    ) -> Successor<G::Node, O>
    where
        V: AstarVisited<G::Node, O>,
    {
        while let Some(child) = parent.children_iter().next() {
            if child == self.ref_graph.end_node() {
                let aln_termination = AlignmentGraphNode::new(child, parent.offset());
                astar_visited.update_score_if_lower(
                    &aln_termination,
                    AlignState::Match,
                    parent.aln_graph_node(),
                    AlignState::Match,
                    self.score,
                );

                return Successor::RefGraphEnd(aln_termination);
            }

            if parent.offset().as_usize() >= self.seq.len() {
                return Successor::QueryEnd(child);
            }

            let child_offset = parent.offset().increase_one();
            let child_node = AlignmentGraphNode::new(child, child_offset);

            if self.ref_graph.is_symbol_equal(child, self.seq[child_offset.as_usize() - 1]) {
                if astar_visited.update_score_if_lower(
                    &child_node,
                    AlignState::Match,
                    parent.aln_graph_node(),
                    AlignState::Match,
                    self.score,
                ) {
                    return Successor::Match(child_node);
                }
            } else {
                return Successor::Mismatch(child_node);
            }
        }

        Successor::SuccessorsExhausted
    }
}
