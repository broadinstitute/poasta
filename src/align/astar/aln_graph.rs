//! Alignment-graph state and expansion: the state machines A* searches over,
//! one per cost-model shape.

use std::collections::HashSet;
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

    fn initial_states<G, O>(&self, ref_graph: &G, seq: &[u8]) -> Vec<AlignmentGraphNode<G::Node, O>>
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        match self.mode() {
            AlignmentMode::Global => vec![AlignmentGraphNode::new(ref_graph.start_node(), O::zero())],
            AlignmentMode::EndsFree { qry_free_begin, graph_free_begin, .. } => {
                ends_free_initial_states(ref_graph, seq, qry_free_begin, graph_free_begin)
            }
        }
    }

    fn is_end<G, O>(
        &self,
        ref_graph: &G,
        seq: &[u8],
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
            AlignmentMode::EndsFree { qry_free_end, graph_free_end, .. } => {
                ends_free_is_end(ref_graph, seq, node, aln_state, qry_free_end, graph_free_end)
            }
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

/// Nodes reachable from the graph start within `max_hops` edges (inclusive of
/// the start node itself), excluding the end sentinel. Used to seed
/// `EndsFree`'s free leading-graph-skip initial states.
///
/// ponytail: recomputed on every `initial_states`/`is_end` call rather than
/// cached on the `AlignmentGraph` struct (which isn't generic over the graph
/// type, so it can't hold a `HashSet<G::Node>` field). Fine for now; bump if
/// EndsFree alignments on large graphs show up in profiles. The fix is to
/// precompute this once in `align()` and thread it through.
fn free_prefix_nodes<G: AlignableGraph>(ref_graph: &G, max_hops: usize) -> HashSet<G::Node> {
    let start = ref_graph.start_node();
    let end = ref_graph.end_node();
    let mut visited = HashSet::new();
    visited.insert(start);
    let mut frontier = vec![start];

    for _ in 0..max_hops {
        let mut next = Vec::new();
        for n in &frontier {
            for succ in ref_graph.successors(*n) {
                if succ != end && visited.insert(succ) {
                    next.push(succ);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }

    visited
}

/// Nodes that can reach the graph end within `max_hops` edges (inclusive of
/// the end node itself), excluding the start sentinel. Used for `EndsFree`'s
/// free trailing-graph-skip end-state check. See `free_prefix_nodes` for the
/// same recompute-per-call caveat.
fn free_suffix_nodes<G: AlignableGraph>(ref_graph: &G, max_hops: usize) -> HashSet<G::Node> {
    let start = ref_graph.start_node();
    let end = ref_graph.end_node();
    let mut visited = HashSet::new();
    visited.insert(end);
    let mut frontier = vec![end];

    for _ in 0..max_hops {
        let mut next = Vec::new();
        for n in &frontier {
            for pred in ref_graph.predecessors(*n) {
                if pred != start && visited.insert(pred) {
                    next.push(pred);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }

    visited
}

/// Shared `EndsFree` `initial_states`: seed every `(node, offset)` combination
/// of {graph nodes within the free leading-graph-skip bound} × {query offsets
/// within the free leading-query-skip bound}, all at score 0.
fn ends_free_initial_states<G, O>(
    ref_graph: &G,
    seq: &[u8],
    qry_free_begin: Bound<usize>,
    graph_free_begin: Bound<usize>,
) -> Vec<AlignmentGraphNode<G::Node, O>>
where
    G: AlignableGraph,
    O: OffsetType,
{
    let max_qry_free = resolve_bound(qry_free_begin, seq.len());
    let max_graph_hops = resolve_bound(graph_free_begin, ref_graph.node_count());
    let free_nodes = free_prefix_nodes(ref_graph, max_graph_hops);

    let mut states = Vec::with_capacity(free_nodes.len() * (max_qry_free + 1));
    for &node in &free_nodes {
        for q in 0..=max_qry_free {
            states.push(AlignmentGraphNode::new(node, O::new(q)));
        }
    }
    states
}

/// Shared `EndsFree` `is_end`: a `Match` state counts as a valid alignment end
/// once its query offset is within the free trailing-query-skip bound of
/// `seq.len()` *and* its graph node can reach the true end within the free
/// trailing-graph-skip bound.
fn ends_free_is_end<G, O>(
    ref_graph: &G,
    seq: &[u8],
    node: &AlignmentGraphNode<G::Node, O>,
    aln_state: AlignState,
    qry_free_end: Bound<usize>,
    graph_free_end: Bound<usize>,
) -> bool
where
    G: AlignableGraph,
    O: OffsetType,
{
    if aln_state != AlignState::Match {
        return false;
    }

    let max_qry_free = resolve_bound(qry_free_end, seq.len());
    let max_graph_hops = resolve_bound(graph_free_end, ref_graph.node_count());

    let qry_ok = node.offset().as_usize() >= seq.len().saturating_sub(max_qry_free);
    let graph_ok = free_suffix_nodes(ref_graph, max_graph_hops).contains(&node.node());

    qry_ok && graph_ok
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
        if visited_data.update_score_if_lower(&ins, AlignState::Match, parent, AlignState::Match, score + ge as u32) {
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
        if visited_data.update_score_if_lower(&del, AlignState::Match, parent, AlignState::Match, score + ge as u32) {
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
        if visited_data.update_score_if_lower(child, AlignState::Match, parent, AlignState::Match, score + mm as u32) {
            f(mm, *child, AlignState::Match);
        }

        let ge = self.costs.gap_extend();
        let ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        if visited_data.update_score_if_lower(&ins, AlignState::Match, parent, AlignState::Match, score + ge as u32) {
            f(ge, ins, AlignState::Match);
        }

        let del = AlignmentGraphNode::new(child.node(), parent.offset());
        if visited_data.update_score_if_lower(&del, AlignState::Match, parent, AlignState::Match, score + ge as u32) {
            f(ge, del, AlignState::Match);
        }
    }
}

// ─── Gap-affine: one or two independent gap pieces ───────────────────────────

/// One gap piece: its costs and the two states (deletion, insertion) it owns.
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
        assert!(n_pieces == 1 || n_pieces == 2, "gap-affine graphs have one or two pieces");
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
        Self { costs, mode, pieces, n_pieces }
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
            if visited_data.update_score_if_lower(&ins, piece.ins, parent, AlignState::Match, opened) {
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
            if visited_data.update_score_if_lower(&del, piece.del, parent, AlignState::Match, opened) {
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
        if visited_data.update_score_if_lower(child, AlignState::Match, parent, AlignState::Match, score + mm as u32) {
            f(mm, *child, AlignState::Match);
        }

        let ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        for piece in self.pieces() {
            let opened = score + piece.open_extend as u32;
            if visited_data.update_score_if_lower(&ins, piece.ins, parent, AlignState::Match, opened) {
                f(piece.open_extend, ins, piece.ins);
            }
        }

        let del = AlignmentGraphNode::new(child.node(), parent.offset());
        for piece in self.pieces() {
            let opened = score + piece.open_extend as u32;
            if visited_data.update_score_if_lower(&del, piece.del, parent, AlignState::Match, opened) {
                f(piece.open_extend, del, piece.del);
            }
        }
    }
}
