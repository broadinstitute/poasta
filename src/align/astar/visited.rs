//! Visited-state storage for the A* engine: the best score found so far for
//! every `(graph node, query offset, alignment state)` triple, plus the
//! bookkeeping bubble-exit pruning needs, plus the backtrace that turns the
//! scores back into an alignment.

use std::marker::PhantomData;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::align::AlignState;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::astar::pruning::ReachedBubbleExitsMatch;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::offsets::OffsetType;
use crate::align::score::INF;
use crate::align::traits::{AlignableGraph, AlignedPair};
use crate::graph::bubbles::index::BubbleIndex;
use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};

/// Insert `offset` into a per-node "reached bubble exit" list, keeping it sorted
/// ascending so `ReachedBubbleExitsMatch` can binary-search it.
#[inline]
fn insert_reached<O: OffsetType>(reached: &mut Vec<O>, offset: O) {
    let idx = reached.partition_point(|v| *v < offset);
    reached.insert(idx, offset);
}

/// State tracked by the A* engine for every `(alignment-graph node, state)`
/// pair it has queued or visited.
pub trait AstarVisited<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> u32;
    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32);

    fn mark_reached(&mut self, score: u32, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState);
    fn dfa_match(&mut self, score: u32, parent: &AlignmentGraphNode<N, O>, child: &AlignmentGraphNode<N, O>);

    fn prune(&self, score: u32, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool;

    /// Update the visited score for `aln_node`/`aln_state` if `score` improves
    /// on what's stored, returning whether it did (i.e. whether the caller
    /// should queue this successor).
    fn update_score_if_lower(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
        parent: &AlignmentGraphNode<N, O>,
        parent_state: AlignState,
        score: u32,
    ) -> bool;

    fn backtrace<G>(&self, ref_graph: &G, seq: &[u8], aln_node: &AlignmentGraphNode<N, O>) -> Vec<AlignedPair<N>>
    where
        G: AlignableGraph<Node = N>;
}

/// Node-rank × query-offset space tiled into `B×B` blocks, allocated on
/// demand. A cell is `[u32; S]`, indexed by `AlignState as usize`, so `S` is
/// `1` for linear, `3` for affine and `5` for two-piece costs.
pub struct BlockedVisitedStorage<N, O, const S: usize, const B: usize = 8> {
    node_blocks: Vec<FxHashMap<O, [[[u32; S]; B]; B]>>,
    node_ranks: Vec<usize>,
    _marker: PhantomData<N>,
}

impl<N, O, const S: usize, const B: usize> BlockedVisitedStorage<N, O, S, B>
where
    N: GraphNodeId,
    O: OffsetType,
{
    pub fn new<G>(graph: &G) -> Self
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        assert!(B.is_power_of_two(), "Block size B should be a power of 2!");

        let node_count = graph.node_count();
        let num_blocks_nodes = (node_count / B) + 1;
        let mut node_ranks = vec![0usize; node_count];
        for node in graph.all_nodes_iter() {
            node_ranks[node.index()] = graph.node_rank(node);
        }

        Self {
            node_blocks: vec![FxHashMap::default(); num_blocks_nodes],
            node_ranks,
            _marker: PhantomData,
        }
    }

    #[inline(always)]
    fn calc_block_ix(&self, aln_node: &AlignmentGraphNode<N, O>) -> (usize, O, usize, usize) {
        let node_rank = self.node_ranks[aln_node.node().index()];
        let node_block = node_rank >> B.ilog2();

        let offset = aln_node.offset().as_usize();
        let offset_block = O::new(offset >> B.ilog2());

        let within_block_node = node_rank & (B - 1);
        let within_block_qry = offset & (B - 1);

        (node_block, offset_block, within_block_node, within_block_qry)
    }

    #[inline]
    pub fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> u32 {
        let (node_block, offset_block, wn, wq) = self.calc_block_ix(aln_node);
        self.node_blocks[node_block]
            .get(&offset_block)
            .map(|v| v[wn][wq][aln_state as usize])
            .unwrap_or(INF)
    }

    #[inline]
    pub fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32) {
        let (node_block, offset_block, wn, wq) = self.calc_block_ix(aln_node);
        let map = &mut self.node_blocks[node_block];
        if map.is_empty() {
            map.reserve(4);
        }
        let v = map.entry(offset_block).or_insert_with(|| [[[INF; S]; B]; B]);
        v[wn][wq][aln_state as usize] = score;
    }

    #[inline]
    pub fn update_score_if_lower_block(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32) -> bool {
        let (node_block, offset_block, wn, wq) = self.calc_block_ix(aln_node);
        let map = &mut self.node_blocks[node_block];
        if map.is_empty() {
            map.reserve(4);
        }
        let v = map.entry(offset_block).or_insert_with(|| [[[INF; S]; B]; B]);
        let cell = &mut v[wn][wq][aln_state as usize];

        if score < *cell {
            *cell = score;
            true
        } else {
            false
        }
    }
}

/// Everything the A* search records per state, for any cost model: scores in
/// blocked storage plus, per bubble exit, the query offsets at which the exit
/// was reached in `Match` (what bubble-exit pruning checks against).
///
/// `S` is the cost model's state count and also selects the backtrace state
/// machine: `1` is the one-state linear model, `3` one gap piece, `5` two.
pub struct AstarData<N, O, C, const S: usize>
where
    N: GraphNodeId,
    O: OffsetType,
{
    costs: C,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    visited: BlockedVisitedStorage<N, O, S>,
    bubbles_reached_m: Vec<Vec<O>>,
}

pub type LinearAstarData<N, O> = AstarData<N, O, Linear, 1>;
pub type AffineAstarData<N, O> = AstarData<N, O, Affine, 3>;
pub type TwoPieceAstarData<N, O> = AstarData<N, O, TwoPieceAffine, 5>;

impl<N, O, C, const S: usize> AstarData<N, O, C, S>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    pub fn new<G>(costs: C, ref_graph: &G, seq: &[u8], bubble_index: Arc<BubbleIndex<N>>) -> Self
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        Self {
            costs,
            seq_len: seq.len(),
            visited: BlockedVisitedStorage::new(ref_graph),
            bubbles_reached_m: vec![Vec::new(); ref_graph.node_count()],
            bubble_index,
        }
    }
}

impl<N, O, C, const S: usize> AstarVisited<N, O> for AstarData<N, O, C, S>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    #[inline]
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> u32 {
        self.visited.get_score(aln_node, aln_state)
    }

    #[inline]
    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32) {
        self.visited.set_score(aln_node, aln_state, score)
    }

    fn mark_reached(&mut self, _score: u32, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) {
        if aln_state == AlignState::Match && self.bubble_index.is_exit(aln_node.node()) {
            insert_reached(&mut self.bubbles_reached_m[aln_node.node().index()], aln_node.offset());
        }
    }

    fn dfa_match(&mut self, score: u32, _parent: &AlignmentGraphNode<N, O>, child: &AlignmentGraphNode<N, O>) {
        self.mark_reached(score, child, AlignState::Match);
    }

    #[inline]
    fn prune(&self, score: u32, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool {
        if !self.bubble_index.node_is_part_of_bubble(aln_node.node()) {
            return false;
        }

        !self.bubble_index.get_node_bubbles(aln_node.node()).iter().all(|bubble| {
            let reached = &self.bubbles_reached_m[bubble.bubble_exit.index()];
            let get_match_score = |node: N, offset: O| self.visited.get_score(&AlignmentGraphNode::new(node, offset), AlignState::Match);

            ReachedBubbleExitsMatch::new(get_match_score, &self.costs, reached, self.seq_len)
                .can_improve_bubble(&self.bubble_index, bubble, aln_node, aln_state, score)
        })
    }

    #[inline]
    fn update_score_if_lower(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
        _parent: &AlignmentGraphNode<N, O>,
        _parent_state: AlignState,
        score: u32,
    ) -> bool {
        self.visited.update_score_if_lower_block(aln_node, aln_state, score)
    }

    fn backtrace<G>(&self, ref_graph: &G, seq: &[u8], aln_node: &AlignmentGraphNode<N, O>) -> Vec<AlignedPair<N>>
    where
        G: AlignableGraph<Node = N>,
    {
        if S == 1 {
            linear_backtrace(&self.visited, &self.costs, ref_graph, seq, aln_node)
        } else {
            gap_backtrace(&self.visited, &self.costs, ref_graph, seq, aln_node)
        }
    }
}

// ─── Backtrace: one-state linear model ───────────────────────────────────────

/// One backtrace step for the linear model. Returns the predecessor and the
/// kind of move that led from it to `curr` (`Match` for the diagonal), or
/// `None` at the start of the alignment.
fn linear_step<N, O, C, G, const S: usize>(
    visited: &BlockedVisitedStorage<N, O, S>,
    costs: &C,
    ref_graph: &G,
    seq: &[u8],
    curr: &AlignmentGraphNode<N, O>,
) -> Option<(AlignmentGraphNode<N, O>, AlignState)>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
    G: AlignableGraph<Node = N>,
{
    let score = visited.get_score(curr, AlignState::Match);
    let get = |node: N, offset: O| visited.get_score(&AlignmentGraphNode::new(node, offset), AlignState::Match);
    let ge = costs.gap_extend() as u32;
    let at_end = curr.node() == ref_graph.end_node();

    if at_end || curr.offset() > O::zero() {
        let (pred_offset, step) = if at_end {
            (curr.offset(), 0)
        } else if ref_graph.is_symbol_equal(curr.node(), seq[curr.offset().as_usize() - 1]) {
            (curr.offset() - O::one(), 0)
        } else {
            (curr.offset() - O::one(), costs.mismatch() as u32)
        };
        for p in ref_graph.predecessors(curr.node()) {
            if get(p, pred_offset) == score.wrapping_sub(step) {
                return Some((AlignmentGraphNode::new(p, pred_offset), AlignState::Match));
            }
        }
    }

    for p in ref_graph.predecessors(curr.node()) {
        if get(p, curr.offset()) == score.wrapping_sub(ge) {
            return Some((AlignmentGraphNode::new(p, curr.offset()), AlignState::Deletion));
        }
    }

    if curr.offset() > O::zero() {
        let pred_offset = curr.offset() - O::one();
        if get(curr.node(), pred_offset) == score.wrapping_sub(ge) {
            return Some((AlignmentGraphNode::new(curr.node(), pred_offset), AlignState::Insertion));
        }
    }

    None
}

fn linear_backtrace<N, O, C, G, const S: usize>(
    visited: &BlockedVisitedStorage<N, O, S>,
    costs: &C,
    ref_graph: &G,
    seq: &[u8],
    end: &AlignmentGraphNode<N, O>,
) -> Vec<AlignedPair<N>>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
    G: AlignableGraph<Node = N>,
{
    let mut alignment = Vec::new();
    let mut curr = *end;

    while let Some((pred, mv)) = linear_step(visited, costs, ref_graph, seq, &curr) {
        match mv {
            // The move into the end sentinel consumes nothing.
            AlignState::Match if curr.node() == ref_graph.end_node() => {}
            AlignState::Match => {
                alignment.push(AlignedPair::new(Some(curr.node()), Some(curr.offset().as_usize() - 1)));
            }
            AlignState::Deletion => alignment.push(AlignedPair::new(Some(curr.node()), None)),
            AlignState::Insertion => alignment.push(AlignedPair::new(None, Some(curr.offset().as_usize() - 1))),
            AlignState::Deletion2 | AlignState::Insertion2 => unreachable!("linear model has one gap piece"),
        }
        curr = pred;
    }

    alignment.reverse();
    alignment
}

// ─── Backtrace: gap-affine models (one or two pieces) ────────────────────────

#[derive(Clone, Copy)]
struct GapPiece {
    open_extend: u32,
    extend: u32,
    del: AlignState,
    ins: AlignState,
}

/// The gap pieces a storage with `S` slots can hold: one for `S == 3`, two
/// for `S == 5`.
fn gap_pieces<C: AlignmentCostModel>(costs: &C, n_states: usize) -> Vec<GapPiece> {
    let mut pieces = vec![GapPiece {
        open_extend: (costs.gap_open() + costs.gap_extend()) as u32,
        extend: costs.gap_extend() as u32,
        del: AlignState::Deletion,
        ins: AlignState::Insertion,
    }];
    if n_states == 5 {
        pieces.push(GapPiece {
            open_extend: (costs.gap_open2() + costs.gap_extend2()) as u32,
            extend: costs.gap_extend2() as u32,
            del: AlignState::Deletion2,
            ins: AlignState::Insertion2,
        });
    }
    pieces
}

/// Find the predecessor `(node, state)` whose score is
/// consistent with a transition into `aln_node`/`aln_state`, or `None` at the
/// start of the alignment.
fn gap_step<N, O, C, G, const S: usize>(
    visited: &BlockedVisitedStorage<N, O, S>,
    costs: &C,
    pieces: &[GapPiece],
    ref_graph: &G,
    seq: &[u8],
    aln_node: &AlignmentGraphNode<N, O>,
    aln_state: AlignState,
) -> Option<(AlignmentGraphNode<N, O>, AlignState)>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
    G: AlignableGraph<Node = N>,
{
    let curr_score = visited.get_score(aln_node, aln_state);
    let get = |node: N, offset: O, state: AlignState| visited.get_score(&AlignmentGraphNode::new(node, offset), state);

    if aln_state == AlignState::Match {
        if aln_node.offset() > O::zero() {
            let at_end = aln_node.node() == ref_graph.end_node();
            let is_match_or_end =
                at_end || ref_graph.is_symbol_equal(aln_node.node(), seq[aln_node.offset().as_usize() - 1]);
            let pred_offset = if at_end { aln_node.offset() } else { aln_node.offset() - O::one() };
            let expected = if is_match_or_end {
                curr_score
            } else {
                curr_score.wrapping_sub(costs.mismatch() as u32)
            };

            for p in ref_graph.predecessors(aln_node.node()) {
                if get(p, pred_offset, AlignState::Match) == expected {
                    return Some((AlignmentGraphNode::new(p, pred_offset), AlignState::Match));
                }
            }
        }

        // Closing a gap into Match is free, so an equal gap score means the
        // gap state is where this Match came from.
        for piece in pieces {
            if visited.get_score(aln_node, piece.del) == curr_score {
                return Some((*aln_node, piece.del));
            }
        }
        for piece in pieces {
            if visited.get_score(aln_node, piece.ins) == curr_score {
                return Some((*aln_node, piece.ins));
            }
        }
        return None;
    }

    let piece = pieces
        .iter()
        .find(|p| p.del == aln_state || p.ins == aln_state)
        .expect("gap state without a matching piece");

    if aln_state == piece.del {
        for p in ref_graph.predecessors(aln_node.node()) {
            if get(p, aln_node.offset(), AlignState::Match) == curr_score.wrapping_sub(piece.open_extend) {
                return Some((AlignmentGraphNode::new(p, aln_node.offset()), AlignState::Match));
            }
        }
        for p in ref_graph.predecessors(aln_node.node()) {
            if get(p, aln_node.offset(), piece.del) == curr_score.wrapping_sub(piece.extend) {
                return Some((AlignmentGraphNode::new(p, aln_node.offset()), piece.del));
            }
        }
    } else if aln_node.offset() > O::zero() {
        let pred = AlignmentGraphNode::new(aln_node.node(), aln_node.offset() - O::one());
        if visited.get_score(&pred, AlignState::Match) == curr_score.wrapping_sub(piece.open_extend) {
            return Some((pred, AlignState::Match));
        }
        if visited.get_score(&pred, piece.ins) == curr_score.wrapping_sub(piece.extend) {
            return Some((pred, piece.ins));
        }
    }

    None
}

fn gap_backtrace<N, O, C, G, const S: usize>(
    visited: &BlockedVisitedStorage<N, O, S>,
    costs: &C,
    ref_graph: &G,
    seq: &[u8],
    aln_node: &AlignmentGraphNode<N, O>,
) -> Vec<AlignedPair<N>>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
    G: AlignableGraph<Node = N>,
{
    let pieces = gap_pieces(costs, S);
    let step = |node: &AlignmentGraphNode<N, O>, state: AlignState| {
        gap_step(visited, costs, &pieces, ref_graph, seq, node, state)
    };

    let Some((mut curr, mut curr_state)) = step(aln_node, AlignState::Match) else {
        panic!("No backtrace for alignment end state?");
    };

    let mut alignment = Vec::new();

    // `aln_node` (the alignment end state) is itself never pushed below: the
    // loop only records `curr` once its predecessor has been found. For
    // `Global` mode that's fine, since the end state is always the virtual
    // end sentinel, which has no base to record. But `EndsFree` alignments
    // can terminate on a real graph node (see `ends_free_is_end`), whose own
    // base would otherwise be silently dropped. Record it here, unless it's
    // just an alias for a just-closed indel at the same (node, offset), which
    // the loop handles once it reaches it.
    let is_indel_alias = curr.node() == aln_node.node() && curr.offset() == aln_node.offset();
    if aln_node.node() != ref_graph.end_node() && !is_indel_alias {
        alignment.push(AlignedPair::new(Some(aln_node.node()), Some(aln_node.offset().as_usize() - 1)));
    }

    while let Some((bt_node, bt_state)) = step(&curr, curr_state) {
        // A Match that came from a gap state at the same cell is the free
        // gap-closing edge; skip it so the cell is recorded once, as the gap.
        if curr_state == AlignState::Match && bt_state != AlignState::Match {
            curr = bt_node;
            curr_state = bt_state;
            continue;
        }

        match curr_state {
            AlignState::Match => {
                alignment.push(AlignedPair::new(Some(curr.node()), Some(curr.offset().as_usize() - 1)));
            }
            AlignState::Insertion | AlignState::Insertion2 => {
                alignment.push(AlignedPair::new(None, Some(curr.offset().as_usize() - 1)));
            }
            AlignState::Deletion | AlignState::Deletion2 => {
                alignment.push(AlignedPair::new(Some(curr.node()), None));
            }
        }

        if bt_node.node() == ref_graph.start_node() {
            break;
        }

        curr = bt_node;
        curr_state = bt_state;
    }

    alignment.reverse();
    alignment
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::mock::create_test_graph1;
    use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

    type O = u32;

    #[test]
    fn get_set_round_trip_on_small_grid() {
        let g = create_test_graph1();
        let mut storage: BlockedVisitedStorage<_, O, 3> = BlockedVisitedStorage::new(&g);

        let n0 = g.rank_to_node(0);
        let node = AlignmentGraphNode::new(n0, 5u32);

        assert_eq!(storage.get_score(&node, AlignState::Match), INF, "unset cell should read INF");

        storage.set_score(&node, AlignState::Match, 7);
        assert_eq!(storage.get_score(&node, AlignState::Match), 7);
        // Other states at the same (node, offset) stay untouched.
        assert_eq!(storage.get_score(&node, AlignState::Deletion), INF);
    }

    #[test]
    fn update_score_if_lower_only_updates_on_improvement() {
        let g = create_test_graph1();
        let mut storage: BlockedVisitedStorage<_, O, 3> = BlockedVisitedStorage::new(&g);

        let n0 = g.rank_to_node(2);
        let node = AlignmentGraphNode::new(n0, 3u32);

        assert!(storage.update_score_if_lower_block(&node, AlignState::Match, 10));
        assert_eq!(storage.get_score(&node, AlignState::Match), 10);

        // Higher score: no update.
        assert!(!storage.update_score_if_lower_block(&node, AlignState::Match, 20));
        assert_eq!(storage.get_score(&node, AlignState::Match), 10);

        // Lower score: updates.
        assert!(storage.update_score_if_lower_block(&node, AlignState::Match, 4));
        assert_eq!(storage.get_score(&node, AlignState::Match), 4);
    }

    #[test]
    fn blocks_spanning_node_and_offset_do_not_alias() {
        let g = create_test_graph1();
        let mut storage: BlockedVisitedStorage<_, O, 3> = BlockedVisitedStorage::new(&g);

        // Exercise several nodes/offsets spanning multiple B=8 blocks.
        let n_ranks = g.node_count().min(10);
        for rank in 0..n_ranks {
            let node = g.rank_to_node(rank);
            for offset in [0u32, 1, 7, 8, 9, 17] {
                let aln_node = AlignmentGraphNode::new(node, offset);
                let score = (rank as u32) * 100 + offset;
                storage.set_score(&aln_node, AlignState::Match, score);
            }
        }

        for rank in 0..n_ranks {
            let node = g.rank_to_node(rank);
            for offset in [0u32, 1, 7, 8, 9, 17] {
                let aln_node = AlignmentGraphNode::new(node, offset);
                let expected = (rank as u32) * 100 + offset;
                assert_eq!(storage.get_score(&aln_node, AlignState::Match), expected);
            }
        }
    }
}
