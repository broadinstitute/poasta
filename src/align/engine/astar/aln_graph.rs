//! Alignment-graph state and expansion, ported from `stable`'s `aligner::aln_graph`
//! + `aligner::scoring::gap_affine`.

use std::collections::HashSet;
use std::ops::Bound;

use crate::align::cost_models::AlignmentCostModel;
use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::engine::astar::visited::AstarVisited;
use crate::align::offsets::OffsetType;
use crate::align::traits::AlignableGraph;
use crate::align::{AlignState, AlignmentMode};
use crate::graph::traits::GraphNodeId;

/// A node in the alignment graph: a cursor into the POA graph paired with a
/// query offset. Does not carry the alignment state (Match/Insertion/...) —
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
/// alignment mode. Implementations have no DP-kernel analogue — this is the
/// A*-specific counterpart to `DPKernel`.
pub trait AlignmentGraph {
    type CostModel: AlignmentCostModel;

    fn get_costs(&self) -> &Self::CostModel;

    fn initial_states<G, O>(&self, ref_graph: &G, seq: &[u8]) -> Vec<AlignmentGraphNode<G::Node, O>>
    where
        G: AlignableGraph,
        O: OffsetType;

    fn is_end<G, O>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        node: &AlignmentGraphNode<G::Node, O>,
        aln_state: AlignState,
    ) -> bool
    where
        G: AlignableGraph,
        O: OffsetType;

    /// Expand all successor states of `node` in state `state`. Calls `f` with
    /// `(score_delta, successor_node, successor_state)` for each successor
    /// whose score was actually improved (`AstarVisited::update_score_if_lower`
    /// returned `true`).
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
    /// mismatch (Match state) plus the indel states opened from `parent`.
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
/// type, so it can't hold a `HashSet<G::Node>` field). Fine for now — bump if
/// EndsFree alignments on large graphs show up in profiles; the fix is to
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

/// Three-state (`Match`/`Deletion`/`Insertion`) alignment graph, generic over
/// any cost model shaped like gap-affine (linear is the `gap_open() == 0`
/// special case of the same recurrence — no separate implementation needed).
pub struct ThreeStateAlignmentGraph<C> {
    costs: C,
    mode: AlignmentMode,
}

impl<C> ThreeStateAlignmentGraph<C>
where
    C: AlignmentCostModel,
{
    pub fn new(costs: C, mode: AlignmentMode) -> Self {
        Self { costs, mode }
    }
}

/// Gap-affine alignment graph (three-state graph specialised to [`Affine`]).
pub type AffineAlignmentGraph = ThreeStateAlignmentGraph<Affine>;

/// Linear-gap alignment graph (three-state graph specialised to [`Linear`]);
/// `Linear::gap_open() == 0` degenerates the affine recurrence into a plain
/// per-base gap cost.
pub type LinearAlignmentGraph = ThreeStateAlignmentGraph<Linear>;

impl<C> AlignmentGraph for ThreeStateAlignmentGraph<C>
where
    C: AlignmentCostModel,
{
    type CostModel = C;

    fn get_costs(&self) -> &Self::CostModel {
        &self.costs
    }

    fn initial_states<G, O>(&self, ref_graph: &G, seq: &[u8]) -> Vec<AlignmentGraphNode<G::Node, O>>
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        match self.mode {
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
        match self.mode {
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
        let go_ge = self.costs.gap_open() + self.costs.gap_extend();

        match state {
            AlignState::Match => {
                let child_offset = node.offset().increase_one();
                for ref_succ in ref_graph.successors(node.node()) {
                    let new_node_mis = AlignmentGraphNode::new(ref_succ, child_offset);

                    let score_delta = if ref_graph.is_symbol_equal(ref_succ, seq[child_offset.as_usize() - 1]) {
                        0u8
                    } else {
                        self.costs.mismatch()
                    };
                    let new_score_mis = score + score_delta as u32;

                    if visited_data.update_score_if_lower(&new_node_mis, AlignState::Match, node, state, new_score_mis) {
                        f(score_delta, new_node_mis, AlignState::Match);
                    }

                    // Open deletion
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + go_ge as u32;
                    if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del) {
                        f(go_ge, new_node_del, AlignState::Deletion);
                    }
                }

                // Open insertion
                let new_node_ins = AlignmentGraphNode::new(node.node(), child_offset);
                let new_score_ins = score + go_ge as u32;
                if child_offset.as_usize() <= seq.len()
                    && visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins)
                {
                    f(go_ge, new_node_ins, AlignState::Insertion);
                }
            }
            AlignState::Insertion => {
                // I -> M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Insertion, score) {
                    f(0, *node, AlignState::Match);
                }

                let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
                let new_score_ins = score + self.costs.gap_extend() as u32;
                if node.offset().as_usize() < seq.len()
                    && visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins)
                {
                    f(self.costs.gap_extend(), new_node_ins, AlignState::Insertion);
                }
            }
            AlignState::Deletion => {
                // D -> M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Deletion, score) {
                    f(0, *node, AlignState::Match);
                }

                for ref_succ in ref_graph.successors(node.node()) {
                    let score_delta = self.costs.gap_extend();
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + score_delta as u32;
                    if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del) {
                        f(score_delta, new_node_del, AlignState::Deletion);
                    }
                }
            }
            AlignState::Insertion2 | AlignState::Deletion2 => {
                panic!("Invalid gap-affine state {state:?}!")
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
        let go_ge = self.costs.gap_open() + self.costs.gap_extend();
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + go_ge as u32;
        if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins) {
            f(go_ge, new_node_ins, AlignState::Insertion);
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
        let go_ge = self.costs.gap_open() + self.costs.gap_extend();
        let new_node_del = AlignmentGraphNode::new(child, parent.offset());
        let new_score_del = score + go_ge as u32;
        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del) {
            f(go_ge, new_node_del, AlignState::Deletion);
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
        let new_score_mis = score + self.costs.mismatch() as u32;
        if visited_data.update_score_if_lower(child, AlignState::Match, parent, AlignState::Match, new_score_mis) {
            f(self.costs.mismatch(), *child, AlignState::Match);
        }

        let go_ge = self.costs.gap_open() + self.costs.gap_extend();

        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + go_ge as u32;
        if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins) {
            f(go_ge, new_node_ins, AlignState::Insertion);
        }

        let new_node_del = AlignmentGraphNode::new(child.node(), parent.offset());
        let new_score_del = score + go_ge as u32;
        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del) {
            f(go_ge, new_node_del, AlignState::Deletion);
        }
    }
}

/// Two-piece gap-affine alignment graph: two fully independent gap tracks,
/// `Deletion`/`Insertion` under piece 1 (go1, ge1) and `Deletion2`/`Insertion2`
/// under piece 2 (go2, ge2). Each track opens directly from `Match` under its
/// own open+extend cost and only ever extends within itself — no piece1 <->
/// piece2 transition edges.
///
/// This intentionally diverges from `stable`'s `Affine2PieceAlignmentGraph`,
/// which only ever opens gaps under piece 1 and reaches piece 2 by
/// transitioning mid-gap. That shape doesn't match this crate's
/// `TwoPieceAffineKernel` DP recurrence (`kernels/two_piece.rs`), which models
/// the two pieces as fully independent tracks, each directly reachable from
/// `M` (`D2[v][q] = min(M[u][q] + go2+ge2, ...)`) with no cross-talk. Since
/// `CanonicalDP` is this crate's oracle (per `CLAUDE.md`), the A* state
/// machine is built to match that recurrence instead of stable's.
pub struct TwoPieceAlignmentGraph {
    costs: TwoPieceAffine,
    mode: AlignmentMode,
}

impl TwoPieceAlignmentGraph {
    pub fn new(costs: TwoPieceAffine, mode: AlignmentMode) -> Self {
        Self { costs, mode }
    }
}

impl AlignmentGraph for TwoPieceAlignmentGraph {
    type CostModel = TwoPieceAffine;

    fn get_costs(&self) -> &Self::CostModel {
        &self.costs
    }

    fn initial_states<G, O>(&self, ref_graph: &G, seq: &[u8]) -> Vec<AlignmentGraphNode<G::Node, O>>
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        match self.mode {
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
        match self.mode {
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
        let go1_ge1 = self.costs.gap_open() + self.costs.gap_extend();
        let go2_ge2 = self.costs.gap_open2() + self.costs.gap_extend2();

        match state {
            AlignState::Match => {
                let child_offset = node.offset().increase_one();
                for ref_succ in ref_graph.successors(node.node()) {
                    let new_node_mis = AlignmentGraphNode::new(ref_succ, child_offset);

                    let score_delta = if ref_graph.is_symbol_equal(ref_succ, seq[child_offset.as_usize() - 1]) {
                        0u8
                    } else {
                        self.costs.mismatch()
                    };
                    let new_score_mis = score + score_delta as u32;

                    if visited_data.update_score_if_lower(&new_node_mis, AlignState::Match, node, state, new_score_mis) {
                        f(score_delta, new_node_mis, AlignState::Match);
                    }

                    // Open deletion, piece 1
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del1 = score + go1_ge1 as u32;
                    if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del1) {
                        f(go1_ge1, new_node_del, AlignState::Deletion);
                    }

                    // Open deletion, piece 2
                    let new_score_del2 = score + go2_ge2 as u32;
                    if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion2, node, state, new_score_del2) {
                        f(go2_ge2, new_node_del, AlignState::Deletion2);
                    }
                }

                if child_offset.as_usize() <= seq.len() {
                    // Open insertion, piece 1
                    let new_node_ins = AlignmentGraphNode::new(node.node(), child_offset);
                    let new_score_ins1 = score + go1_ge1 as u32;
                    if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins1) {
                        f(go1_ge1, new_node_ins, AlignState::Insertion);
                    }

                    // Open insertion, piece 2
                    let new_score_ins2 = score + go2_ge2 as u32;
                    if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion2, node, state, new_score_ins2) {
                        f(go2_ge2, new_node_ins, AlignState::Insertion2);
                    }
                }
            }
            AlignState::Insertion => {
                // I1 -> M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Insertion, score) {
                    f(0, *node, AlignState::Match);
                }

                let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
                let new_score_ins = score + self.costs.gap_extend() as u32;
                if node.offset().as_usize() < seq.len()
                    && visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins)
                {
                    f(self.costs.gap_extend(), new_node_ins, AlignState::Insertion);
                }
            }
            AlignState::Insertion2 => {
                // I2 -> M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Insertion2, score) {
                    f(0, *node, AlignState::Match);
                }

                let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
                let new_score_ins = score + self.costs.gap_extend2() as u32;
                if node.offset().as_usize() < seq.len()
                    && visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion2, node, state, new_score_ins)
                {
                    f(self.costs.gap_extend2(), new_node_ins, AlignState::Insertion2);
                }
            }
            AlignState::Deletion => {
                // D1 -> M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Deletion, score) {
                    f(0, *node, AlignState::Match);
                }

                for ref_succ in ref_graph.successors(node.node()) {
                    let score_delta = self.costs.gap_extend();
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + score_delta as u32;
                    if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del) {
                        f(score_delta, new_node_del, AlignState::Deletion);
                    }
                }
            }
            AlignState::Deletion2 => {
                // D2 -> M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Deletion2, score) {
                    f(0, *node, AlignState::Match);
                }

                for ref_succ in ref_graph.successors(node.node()) {
                    let score_delta = self.costs.gap_extend2();
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + score_delta as u32;
                    if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion2, node, state, new_score_del) {
                        f(score_delta, new_node_del, AlignState::Deletion2);
                    }
                }
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
        let go1_ge1 = self.costs.gap_open() + self.costs.gap_extend();
        let go2_ge2 = self.costs.gap_open2() + self.costs.gap_extend2();

        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins1 = score + go1_ge1 as u32;
        if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins1) {
            f(go1_ge1, new_node_ins, AlignState::Insertion);
        }

        let new_score_ins2 = score + go2_ge2 as u32;
        if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion2, parent, AlignState::Match, new_score_ins2) {
            f(go2_ge2, new_node_ins, AlignState::Insertion2);
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
        let go1_ge1 = self.costs.gap_open() + self.costs.gap_extend();
        let go2_ge2 = self.costs.gap_open2() + self.costs.gap_extend2();

        let new_node_del = AlignmentGraphNode::new(child, parent.offset());
        let new_score_del1 = score + go1_ge1 as u32;
        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del1) {
            f(go1_ge1, new_node_del, AlignState::Deletion);
        }

        let new_score_del2 = score + go2_ge2 as u32;
        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion2, parent, AlignState::Match, new_score_del2) {
            f(go2_ge2, new_node_del, AlignState::Deletion2);
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
        let new_score_mis = score + self.costs.mismatch() as u32;
        if visited_data.update_score_if_lower(child, AlignState::Match, parent, AlignState::Match, new_score_mis) {
            f(self.costs.mismatch(), *child, AlignState::Match);
        }

        let go1_ge1 = self.costs.gap_open() + self.costs.gap_extend();
        let go2_ge2 = self.costs.gap_open2() + self.costs.gap_extend2();

        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins1 = score + go1_ge1 as u32;
        if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins1) {
            f(go1_ge1, new_node_ins, AlignState::Insertion);
        }
        let new_score_ins2 = score + go2_ge2 as u32;
        if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion2, parent, AlignState::Match, new_score_ins2) {
            f(go2_ge2, new_node_ins, AlignState::Insertion2);
        }

        let new_node_del = AlignmentGraphNode::new(child.node(), parent.offset());
        let new_score_del1 = score + go1_ge1 as u32;
        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del1) {
            f(go1_ge1, new_node_del, AlignState::Deletion);
        }
        let new_score_del2 = score + go2_ge2 as u32;
        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion2, parent, AlignState::Match, new_score_del2) {
            f(go2_ge2, new_node_del, AlignState::Deletion2);
        }
    }
}
