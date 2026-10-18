//! Alignment-graph state and expansion: the state machines A* searches over,
//! one per cost-model shape.

use std::ops::Bound;

use crate::align::astar::visited::AstarVisited;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::offsets::OffsetType;
use crate::align::traits::AlignableGraph;
use crate::align::{AlignState, AlignmentMode};
use crate::graph::traits::GraphNodeId;

/// A node in the alignment graph: a cursor into the POA graph paired with a
/// query offset. Does not carry the alignment state (Match/Insertion/...);
/// that is tracked separately by callers (queue items, visited cells).
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub struct AlignmentGraphNode<N, O> {
    node: N,
    offset: O,
}

impl<N, O> AlignmentGraphNode<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    pub fn new(node: N, offset: O) -> Self {
        Self { node, offset }
    }

    #[inline(always)]
    pub fn node(&self) -> N {
        self.node
    }

    #[inline(always)]
    pub fn offset(&self) -> O {
        self.offset
    }
}

/// Expands alignment-graph states according to a particular cost model and
/// alignment mode.
///
/// Every expansion method calls `f` with `(score_delta, successor_node,
/// successor_state)` for each successor whose score improved (that is, when
/// `AstarVisited::update_score_if_lower` returned `true`), so the caller can
/// queue exactly those. The main loop only expands `Match` states through
/// the greedy match DFA, which reports where extension stopped and why; the
/// three `expand_*` methods below handle each of those stops. Gap states are
/// expanded through `expand_all`.
pub trait AlignmentGraph {
    fn mode(&self) -> AlignmentMode;

    fn initial_states<G, O>(
        &self,
        ref_graph: &G,
        region: &EndsFreeRegion<G::Node>,
    ) -> Vec<AlignmentGraphNode<G::Node, O>>
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        match self.mode() {
            AlignmentMode::Global => {
                vec![AlignmentGraphNode::new(ref_graph.start_node(), O::zero())]
            }
            AlignmentMode::EndsFree { .. } => region.initial_states(),
        }
    }

    fn is_end<G, O>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        region: &EndsFreeRegion<G::Node>,
        node: &AlignmentGraphNode<G::Node, O>,
        aln_state: AlignState,
    ) -> bool
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        match self.mode() {
            AlignmentMode::Global => {
                aln_state == AlignState::Match
                    && node.node() == ref_graph.end_node()
                    && node.offset().as_usize() == seq.len()
            }
            AlignmentMode::EndsFree { .. } => region.is_end(node, aln_state),
        }
    }

    /// Expand all successor states of `node` in gap state `state`.
    #[allow(clippy::too_many_arguments)]
    fn expand_all<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: u32,
        node: &AlignmentGraphNode<G::Node, O>,
        state: AlignState,
        f: F,
    ) where
        V: AstarVisited<G::Node, O>,
        G: AlignableGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::Node, O>, AlignState);

    /// Expand from a node at the reference-graph end sentinel: only an
    /// insertion (extra query bases) is possible.
    fn expand_ref_graph_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        score: u32,
        f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);

    /// Expand from a node at the query end: only a deletion (remaining graph
    /// nodes with no query bases left) is possible.
    fn expand_query_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: N,
        score: u32,
        f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);

    /// Expand a mismatch found during greedy match extension: queues the
    /// mismatch (Match state) plus the indels opened from `parent`.
    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: u32,
        f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);
}

// ─── EndsFree (semi-global) helpers, shared across cost-model shapes ────────

/// Resolve a `Bound<usize>` (as used in `AlignmentMode::EndsFree`) into a
/// concrete maximum free length, capped at `max`.
pub(crate) fn resolve_bound(bound: Bound<usize>, max: usize) -> usize {
    match bound {
        Bound::Unbounded => max,
        Bound::Included(k) => k.min(max),
        Bound::Excluded(k) => k.saturating_sub(1).min(max),
    }
}

/// The free regions of an `EndsFree` alignment, computed once per
/// alignment. Empty in `Global` mode.
pub struct EndsFreeRegion<N> {
    /// Nodes within the free leading-graph bound of the start, in rank
    /// order, so the queue order does not depend on a hasher.
    start_nodes: Vec<N>,
    /// Largest free leading query offset.
    max_qry_free: usize,
    /// Lowest query offset that counts as an end.
    end_offset_min: usize,
    /// Nodes that reach the end within the free trailing-graph bound, by
    /// node index.
    end_node: Vec<bool>,
}

impl<N: GraphNodeId> EndsFreeRegion<N> {
    pub fn new<G>(ref_graph: &G, seq: &[u8], mode: AlignmentMode) -> Self
    where
        G: AlignableGraph<Node = N>,
    {
        let AlignmentMode::EndsFree {
            qry_free_begin,
            qry_free_end,
            graph_free_begin,
            graph_free_end,
        } = mode
        else {
            return Self {
                start_nodes: Vec::new(),
                max_qry_free: 0,
                end_offset_min: 0,
                end_node: Vec::new(),
            };
        };
        let n = ref_graph.node_count();
        let (start, end) = (ref_graph.start_node(), ref_graph.end_node());

        // Fewest edges from the start, then to the end, by one pass over the
        // ranks each. Successors have higher ranks.
        let mut hops = vec![usize::MAX; n];
        hops[start.index()] = 0;
        let max_hops = resolve_bound(graph_free_begin, n);
        let mut start_nodes = Vec::new();
        for rank in 0..n {
            let v = ref_graph.rank_to_node(rank);
            let d = hops[v.index()];
            if d == usize::MAX || v == end {
                continue;
            }
            if d <= max_hops {
                start_nodes.push(v);
            }
            for succ in ref_graph.successors(v) {
                hops[succ.index()] = hops[succ.index()].min(d + 1);
            }
        }

        let mut to_end = vec![usize::MAX; n];
        to_end[end.index()] = 0;
        let max_hops = resolve_bound(graph_free_end, n);
        let mut end_node = vec![false; n];
        for rank in (0..n).rev() {
            let v = ref_graph.rank_to_node(rank);
            let d = to_end[v.index()];
            if d == usize::MAX || v == start {
                continue;
            }
            end_node[v.index()] = d <= max_hops;
            for pred in ref_graph.predecessors(v) {
                to_end[pred.index()] = to_end[pred.index()].min(d + 1);
            }
        }

        Self {
            start_nodes,
            max_qry_free: resolve_bound(qry_free_begin, seq.len()),
            end_offset_min: seq
                .len()
                .saturating_sub(resolve_bound(qry_free_end, seq.len())),
            end_node,
        }
    }

    /// Every free start node at every free leading query offset, score 0.
    fn initial_states<O: OffsetType>(&self) -> Vec<AlignmentGraphNode<N, O>> {
        let mut states = Vec::with_capacity(self.start_nodes.len() * (self.max_qry_free + 1));
        for &node in &self.start_nodes {
            for q in 0..=self.max_qry_free {
                states.push(AlignmentGraphNode::new(node, O::new(q)));
            }
        }
        states
    }

    /// A `Match` state is an end once its query offset is within the free
    /// trailing query bound and its node reaches the graph end within the
    /// free trailing graph bound.
    fn is_end<O: OffsetType>(
        &self,
        node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> bool {
        aln_state == AlignState::Match
            && node.offset().as_usize() >= self.end_offset_min
            && self.end_node[node.node().index()]
    }
}

// ─── Linear: one state ───────────────────────────────────────────────────────

/// Linear-gap alignment graph. With no gap-open cost there is nothing to
/// remember between steps, so every state is `Match` and a gap is just an
/// edge that costs `gap_extend` and moves along one axis.
pub struct LinearAlignmentGraph {
    costs: Linear,
    mode: AlignmentMode,
}

impl LinearAlignmentGraph {
    pub fn new(costs: Linear, mode: AlignmentMode) -> Self {
        Self { costs, mode }
    }
}

impl AlignmentGraph for LinearAlignmentGraph {
    fn mode(&self) -> AlignmentMode {
        self.mode
    }

    fn expand_all<V, G, O, F>(
        &self,
        _visited_data: &mut V,
        _ref_graph: &G,
        _seq: &[u8],
        _score: u32,
        _node: &AlignmentGraphNode<G::Node, O>,
        state: AlignState,
        _f: F,
    ) where
        V: AstarVisited<G::Node, O>,
        G: AlignableGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::Node, O>, AlignState),
    {
        unreachable!("linear model has no gap states to expand, got {state:?}")
    }

    fn expand_ref_graph_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        score: u32,
        mut f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState),
    {
        let ge = self.costs.gap_extend();
        let ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        if visited_data.update_score_if_lower(
            &ins,
            AlignState::Match,
            parent,
            AlignState::Match,
            score + ge as u32,
        ) {
            f(ge, ins, AlignState::Match);
        }
    }

    fn expand_query_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: N,
        score: u32,
        mut f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState),
    {
        let ge = self.costs.gap_extend();
        let del = AlignmentGraphNode::new(child, parent.offset());
        if visited_data.update_score_if_lower(
            &del,
            AlignState::Match,
            parent,
            AlignState::Match,
            score + ge as u32,
        ) {
            f(ge, del, AlignState::Match);
        }
    }

    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: u32,
        mut f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState),
    {
        let mm = self.costs.mismatch();
        if visited_data.update_score_if_lower(
            child,
            AlignState::Match,
            parent,
            AlignState::Match,
            score + mm as u32,
        ) {
            f(mm, *child, AlignState::Match);
        }

        let ge = self.costs.gap_extend();
        let ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        if visited_data.update_score_if_lower(
            &ins,
            AlignState::Match,
            parent,
            AlignState::Match,
            score + ge as u32,
        ) {
            f(ge, ins, AlignState::Match);
        }

        let del = AlignmentGraphNode::new(child.node(), parent.offset());
        if visited_data.update_score_if_lower(
            &del,
            AlignState::Match,
            parent,
            AlignState::Match,
            score + ge as u32,
        ) {
            f(ge, del, AlignState::Match);
        }
    }
}

// ─── Gap-affine: one or two independent gap pieces ───────────────────────────

/// A gap piece has its own costs and its own two states, deletion and
/// insertion.
#[derive(Clone, Copy)]
struct GapPiece {
    open_extend: u8,
    extend: u8,
    del: AlignState,
    ins: AlignState,
}

/// Gap-affine alignment graph with `Match` plus one deletion and one
/// insertion state per gap piece. Each piece opens directly from `Match`
/// under its own open+extend cost and only ever extends within itself; the
/// zero-cost edge back to `Match` closes it.
///
/// With two pieces (two-piece affine) this means no piece-to-piece
/// transitions, matching the "two independent tracks" recurrence
/// `D2 = min(M + go2 + ge2, D2 + ge2)` rather than a model where piece two is
/// only reachable from piece one mid-gap.
pub struct GapAlignmentGraph<C> {
    costs: C,
    mode: AlignmentMode,
    pieces: [GapPiece; 2],
    n_pieces: usize,
}

pub type AffineAlignmentGraph = GapAlignmentGraph<Affine>;
pub type TwoPieceAlignmentGraph = GapAlignmentGraph<TwoPieceAffine>;

impl<C> GapAlignmentGraph<C>
where
    C: AlignmentCostModel,
{
    pub fn new(costs: C, mode: AlignmentMode, n_pieces: usize) -> Self {
        assert!(
            n_pieces == 1 || n_pieces == 2,
            "gap-affine graphs have one or two pieces"
        );
        let pieces = [
            GapPiece {
                open_extend: costs.gap_open() + costs.gap_extend(),
                extend: costs.gap_extend(),
                del: AlignState::Deletion,
                ins: AlignState::Insertion,
            },
            GapPiece {
                open_extend: costs.gap_open2() + costs.gap_extend2(),
                extend: costs.gap_extend2(),
                del: AlignState::Deletion2,
                ins: AlignState::Insertion2,
            },
        ];
        Self {
            costs,
            mode,
            pieces,
            n_pieces,
        }
    }

    #[inline(always)]
    fn pieces(&self) -> &[GapPiece] {
        &self.pieces[..self.n_pieces]
    }
}

impl<C> AlignmentGraph for GapAlignmentGraph<C>
where
    C: AlignmentCostModel,
{
    fn mode(&self) -> AlignmentMode {
        self.mode
    }

    fn expand_all<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: u32,
        node: &AlignmentGraphNode<G::Node, O>,
        state: AlignState,
        mut f: F,
    ) where
        V: AstarVisited<G::Node, O>,
        G: AlignableGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::Node, O>, AlignState),
    {
        let piece = match state {
            AlignState::Deletion | AlignState::Insertion => self.pieces[0],
            AlignState::Deletion2 | AlignState::Insertion2 => self.pieces[1],
            AlignState::Match => unreachable!("Match states are expanded through the match DFA"),
        };

        // Closing the gap is free.
        if visited_data.update_score_if_lower(node, AlignState::Match, node, state, score) {
            f(0, *node, AlignState::Match);
        }

        let extended = score + piece.extend as u32;
        if state == piece.del {
            for ref_succ in ref_graph.successors(node.node()) {
                let del = AlignmentGraphNode::new(ref_succ, node.offset());
                if visited_data.update_score_if_lower(&del, state, node, state, extended) {
                    f(piece.extend, del, state);
                }
            }
        } else if node.offset().as_usize() < seq.len() {
            let ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
            if visited_data.update_score_if_lower(&ins, state, node, state, extended) {
                f(piece.extend, ins, state);
            }
        }
    }

    fn expand_ref_graph_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        score: u32,
        mut f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState),
    {
        let ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        for piece in self.pieces() {
            let opened = score + piece.open_extend as u32;
            if visited_data.update_score_if_lower(
                &ins,
                piece.ins,
                parent,
                AlignState::Match,
                opened,
            ) {
                f(piece.open_extend, ins, piece.ins);
            }
        }
    }

    fn expand_query_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: N,
        score: u32,
        mut f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState),
    {
        let del = AlignmentGraphNode::new(child, parent.offset());
        for piece in self.pieces() {
            let opened = score + piece.open_extend as u32;
            if visited_data.update_score_if_lower(
                &del,
                piece.del,
                parent,
                AlignState::Match,
                opened,
            ) {
                f(piece.open_extend, del, piece.del);
            }
        }
    }

    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: u32,
        mut f: F,
    ) where
        V: AstarVisited<N, O>,
        N: GraphNodeId,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState),
    {
        let mm = self.costs.mismatch();
        if visited_data.update_score_if_lower(
            child,
            AlignState::Match,
            parent,
            AlignState::Match,
            score + mm as u32,
        ) {
            f(mm, *child, AlignState::Match);
        }

        let ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        for piece in self.pieces() {
            let opened = score + piece.open_extend as u32;
            if visited_data.update_score_if_lower(
                &ins,
                piece.ins,
                parent,
                AlignState::Match,
                opened,
            ) {
                f(piece.open_extend, ins, piece.ins);
            }
        }

        let del = AlignmentGraphNode::new(child.node(), parent.offset());
        for piece in self.pieces() {
            let opened = score + piece.open_extend as u32;
            if visited_data.update_score_if_lower(
                &del,
                piece.del,
                parent,
                AlignState::Match,
                opened,
            ) {
                f(piece.open_extend, del, piece.del);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::alignment::AddAlignment;
    use crate::graph::poa::POAGraph;
    use crate::graph::traits::GraphWithNodeOrdering;

    fn chain(seq: &[u8]) -> POAGraph<u32> {
        let mut g = POAGraph::new();
        g.add_alignment("s", seq, None, &vec![1; seq.len()])
            .unwrap();
        g
    }

    fn ends_free(begin: Bound<usize>, end: Bound<usize>) -> AlignmentMode {
        AlignmentMode::EndsFree {
            qry_free_begin: Bound::Excluded(0),
            qry_free_end: Bound::Excluded(0),
            graph_free_begin: begin,
            graph_free_end: end,
        }
    }

    #[test]
    fn free_start_nodes_come_in_rank_order() {
        let g = chain(b"ACGTACGT");
        let region = EndsFreeRegion::new(&g, b"ACG", ends_free(Bound::Unbounded, Bound::Unbounded));
        let ranks: Vec<usize> = region.start_nodes.iter().map(|&v| g.node_rank(v)).collect();
        assert!(ranks.windows(2).all(|w| w[0] < w[1]));
        // The start sentinel and all eight real nodes, not the end sentinel.
        assert_eq!(ranks.len(), 9);
        assert_eq!(region.end_node.iter().filter(|&&b| b).count(), 9);
    }

    #[test]
    fn hop_bounds_limit_both_ends() {
        let g = chain(b"ACGTACGT");
        let region = EndsFreeRegion::new(
            &g,
            b"ACG",
            ends_free(Bound::Included(2), Bound::Included(1)),
        );
        // The start sentinel and the first two nodes.
        assert_eq!(region.start_nodes.len(), 3);
        // The end sentinel and the last node.
        assert_eq!(region.end_node.iter().filter(|&&b| b).count(), 2);
        assert_eq!(region.initial_states::<u32>().len(), 3);
    }
}
