use std::cell::{RefCell, RefMut};
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::astar::AstarVisited;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::Score;
use crate::graphs::{AlignableRefGraph, NodeIndexType};

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
struct StackNode<N, O, I>(AlignmentGraphNode<N, O>, RefCell<I>)
where
    N: NodeIndexType,
    O: OffsetType;

impl<N, O, I> StackNode<N, O, I>
where
    N: NodeIndexType,
    O: OffsetType,
    I: Iterator,
{
    fn new(state: AlignmentGraphNode<N, O>, iter: I)  -> Self {
        Self(state, RefCell::new(iter))
    }

    #[inline(always)]
    fn node(&self) -> N {
        self.0.node()
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
    fn children_iter(&self) -> RefMut<I> {
        self.1.borrow_mut()
    }
}

#[derive(Debug)]
pub enum ExtendResult<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    /// Variant indicating that we reached the end node in the reference graph
    ///
    /// The associated values represent the edge traversed in the alignment graph
    /// and thus represents the parent alignment graph node and the
    /// just reached alignment graph node with N being the reference graph end node.
    RefGraphEnd(AlignmentGraphNode<N, O>, AlignmentGraphNode<N, O>),

    /// Variant indicating that we could not extend any further because all query sequence
    /// has been aligned.
    ///
    /// The associated value is the parent alignment graph node and the
    /// child reference graph node that we tried to align, but couldn't.
    QueryEnd(AlignmentGraphNode<N, O>, N),


    /// During extension, we encountered a reference graph node
    /// with a different symbol than the corresponding query symbol.
    ///
    /// The associated values indicate the edge in the alignment graph traversed, with the
    /// first ['AlignmentGraphNode'] the last node with matching symbols between the
    /// reference and the query, and the second ['AlignmentGraphNode'] the mismatching
    /// alignment graph node.
    ///
    /// We return the edge such that we can open indel states from the last matching
    /// alignment graph node (the first value), in addition to the mismatch state
    /// (the second value).
    Mismatch(AlignmentGraphNode<N, O>, AlignmentGraphNode<N, O>),
}

/// This struct represents the state for POASTA's
/// depth-first greedy alignment algorithm.
///
/// It traverses the reference graph in a depth-first
/// manner, and greedily aligns visited nodes with
/// the corresponding query offset if the characters
/// match.
///
/// This works because we set the alignment cost of
/// matching characters to 0.
pub struct DepthFirstGreedyAlignment<'a, G, O>
    where G: AlignableRefGraph + 'a,
          O: OffsetType,
{
    /// The alignment graph
    ref_graph: &'a G,

    /// Query sequence to align
    seq: &'a [u8],

    /// The current alignment score
    score: Score,

    /// Number of alignment states visited
    num_visited: usize,

    /// Stack for depth-first alignment of matches between query and graph
    stack: Vec<StackNode<G::NodeIndex, O, G::SuccessorIterator<'a>>>,
}

impl<'a, G, O> DepthFirstGreedyAlignment<'a, G, O>
where
    G: AlignableRefGraph,
    O: OffsetType,
{
    pub fn new(
        ref_graph: &'a G,
        seq: &'a [u8],
        score: Score,
        start_node: &AlignmentGraphNode<G::NodeIndex, O>,
    ) -> Self {
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

    pub fn extend<V>(
        &mut self,
        astar_visited: &mut V,
    ) -> Option<ExtendResult<G::NodeIndex, O>>
    where
        V: AstarVisited<G::NodeIndex, O>,
    {
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();

            match self.next_valid_successor(parent, astar_visited) {
                Successor::RefGraphEnd(aln_node) => return Some(ExtendResult::RefGraphEnd(
                    *parent.aln_graph_node(),
                    aln_node
                )),
                Successor::QueryEnd(ref_node) => return Some(ExtendResult::QueryEnd(
                    *parent.aln_graph_node(),
                    ref_node
                )),
                Successor::Match(child) => {
                    if astar_visited.prune(self.score, &child, AlignState::Match) {
                        continue;
                    }

                    astar_visited.visit(self.score, &child, AlignState::Match);
                    self.num_visited += 1;

                    let child_succ = self.ref_graph.successors(child.node());
                    self.stack.push(StackNode::new(child, child_succ));
                },
                Successor::Mismatch(child) => {
                    return Some(ExtendResult::Mismatch(
                        *parent.aln_graph_node(),
                        child
                    ))
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
        parent: &StackNode<G::NodeIndex, O, G::SuccessorIterator<'a>>,
        // Separate mutable borrows for  visited data and reached bubble exits
        // to prevent requiring mutable &self
        astar_visited: &mut V,
    ) -> Successor<G::NodeIndex, O>
    where
        V: AstarVisited<G::NodeIndex, O>,
    {
        while let Some(child) = parent.children_iter().next() {
            if child == self.ref_graph.end_node() {
                let aln_termination = AlignmentGraphNode::new(child, parent.offset());
                astar_visited.update_score_if_lower(&aln_termination, AlignState::Match,
                                                    parent.aln_graph_node(), AlignState::Match, self.score);

                return Successor::RefGraphEnd(aln_termination);
            }

            if parent.offset().as_usize() >= self.seq.len() {
                return Successor::QueryEnd(child)
            }

            let child_offset = parent.offset().increase_one();
            let child_node = AlignmentGraphNode::new(child, child_offset);

            if self.ref_graph.is_symbol_equal(child, self.seq[child_offset.as_usize()-1]) {
                if astar_visited
                    .update_score_if_lower(&child_node, AlignState::Match,
                                           parent.aln_graph_node(), AlignState::Match, self.score)
                {
                    return Successor::Match(child_node)
                }
            } else {
                return Successor::Mismatch(child_node)
            }
        }

        Successor::SuccessorsExhausted
    }
}

enum Successor<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    RefGraphEnd(AlignmentGraphNode<N, O>),
    QueryEnd(N),
    Match(AlignmentGraphNode<N, O>),
    Mismatch(AlignmentGraphNode<N, O>),
    SuccessorsExhausted,
}

