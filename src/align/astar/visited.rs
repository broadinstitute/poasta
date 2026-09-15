//! Visited-state storage for the A* engine: the best score found so far for
//! every `(graph node, query offset, alignment state)` triple, the reached
//! bits bubble-exit pruning reads, and the backtrace that turns the scores
//! back into an alignment.
//!
//! Scores live in `B×B` blocks, one window of blocks over the query offsets
//! per group of `B` ranks. Reached `Match` states are one bit per query
//! offset per rank.

use std::marker::PhantomData;
use std::sync::Arc;

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

    /// Forget the previous alignment and get ready for `seq` on `graph`,
    /// keeping the allocations.
    fn reset<G>(&mut self, graph: &G, seq: &[u8], bubbles: Arc<BubbleIndex<N>>)
    where
        G: GraphWithNodeOrdering<NodeType = N>;
}

/// The blocks of one group of `B` ranks over a window of query offsets.
///
/// A window instead of a hash map because the offsets the search touches
/// for one rank form a band around the diagonal. A lookup is then a
/// subtraction, a bounds check and a load, where the hashed blocks probed
/// and often rehashed on every read and write. Memory is still
/// proportional to the band.
///
/// The blocks stay square. A* expands neighbours: the successors of a
/// state sit one rank down and at most one offset over, so with `B×B`
/// blocks the cells one pop touches share a block and a few cache lines.
pub struct BlockRow<const S: usize, const B: usize> {
    /// Offset block index of `blocks[0]`. Block `i` covers offsets
    /// `(base + i) * B .. (base + i + 1) * B`.
    base: u32,
    /// Rank-major inside a block, `[u32; S]` per cell, `INF` when unreached.
    blocks: Vec<[[[u32; S]; B]; B]>,
}

impl<const S: usize, const B: usize> BlockRow<S, B> {
    const EMPTY: [[[u32; S]; B]; B] = [[[INF; S]; B]; B];

    fn new() -> Self {
        Self { base: 0, blocks: Vec::new() }
    }

    #[inline]
    fn get(&self, rank_in_block: usize, offset: u32, state: usize) -> u32 {
        // `wrapping_sub` folds "below the window" into "past the end", so one
        // bounds check covers both.
        let block = ((offset >> B.ilog2()) as usize).wrapping_sub(self.base as usize);
        match self.blocks.get(block) {
            Some(b) => b[rank_in_block][offset as usize & (B - 1)][state],
            None => INF,
        }
    }

    #[inline]
    fn cell_mut(&mut self, rank_in_block: usize, offset: u32) -> &mut [u32; S] {
        let block = (offset >> B.ilog2()) as usize;
        let base = self.base as usize;
        if block < base || block >= base + self.blocks.len() {
            self.grow_to(block);
        }
        &mut self.blocks[block - self.base as usize][rank_in_block][offset as usize & (B - 1)]
    }

    // The first window is small on purpose. On similar sequences a rank is
    // touched at one to three offsets, and a graph with tens of thousands of
    // ranks would pay a large fill per alignment for a wide first window.
    // Growth doubles in either direction, so the total copy cost stays
    // linear in the final window. The base is a whole unit (a block, a
    // word), so growing downward is a prepend and never shifts anything
    // inside a unit.
    #[cold]
    #[inline(never)]
    fn grow_to(&mut self, block: usize) {
        let base = self.base as usize;
        if self.blocks.is_empty() {
            self.base = block as u32;
            self.blocks.resize(2, Self::EMPTY);
        } else if block < base {
            let grow = (base - block).max(self.blocks.len()).min(base);
            self.blocks.splice(0..0, std::iter::repeat_n(Self::EMPTY, grow));
            self.base -= grow as u32;
        } else {
            let needed = block - base + 1;
            self.blocks.resize(needed.max(2 * self.blocks.len()), Self::EMPTY);
        }
    }
}

/// One bit per query offset at which the `Match` state of one rank was
/// reached.
///
/// This replaces a sorted `Vec` of offsets per bubble exit. The pruning
/// check needs the reached offsets just below, inside and just above a
/// window. On high-error reads an exit is reached at hundreds to
/// thousands of offsets, so the sorted vector cost two binary searches
/// per check and a memmove per reached exit. That was half the run. With
/// bits each lookup is a masked word scan.
///
/// The bits are set for every reached `Match` state, not only at bubble
/// exits: setting a bit is cheaper than testing whether the node is an
/// exit, and a uniform meaning is what a unipath extension needs later.
/// Only exit ranks are read.
#[derive(Default)]
pub struct ReachedBits {
    /// Query offset of bit 0 of `words[0]`, always a multiple of 64.
    base: u32,
    words: Vec<u64>,
    /// Set bits in `words`. The pruning check returns early while it is
    /// zero, which on similar sequences is most checks: an exit is reached at
    /// one offset and most exits have not been reached yet when a state
    /// inside their bubble is expanded.
    n_reached: u32,
}

// Every scan below is over absolute offsets. `i = offset - base` is the bit
// index in the window, `i >> 6` the word and `i & 63` the bit. The masks:
// `u64::MAX >> (63 - k)` keeps bits `0..=k` and `u64::MAX << k` keeps bits
// `k..=63`, both with a shift in `0..=63`.
impl ReachedBits {
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n_reached == 0
    }

    /// One past the last offset the window covers.
    fn end(&self) -> u32 {
        self.base + 64 * self.words.len() as u32
    }

    pub fn mark(&mut self, offset: u32) {
        let word = (offset >> 6) as usize;
        let base_word = (self.base >> 6) as usize;
        // Grows like `BlockRow::grow_to`, one word at first, doubling after.
        if self.words.is_empty() {
            self.base = offset & !63;
            self.words.push(0);
        } else if word < base_word {
            let grow = (base_word - word).max(self.words.len()).min(base_word);
            self.words.splice(0..0, std::iter::repeat_n(0, grow));
            self.base -= 64 * grow as u32;
        } else if word >= base_word + self.words.len() {
            let needed = word - base_word + 1;
            self.words.resize(needed.max(2 * self.words.len()), 0);
        }
        let slot = &mut self.words[word - (self.base >> 6) as usize];
        let bit = 1u64 << (offset & 63);
        if *slot & bit == 0 {
            *slot |= bit;
            self.n_reached += 1;
        }
    }

    /// Largest reached offset `<= offset`.
    pub fn last_reached_at_or_below(&self, offset: u32) -> Option<u32> {
        if self.words.is_empty() || offset < self.base {
            return None;
        }
        let i = (offset.min(self.end() - 1) - self.base) as usize;
        let mut w = i >> 6;
        let mut word = self.words[w] & (u64::MAX >> (63 - (i & 63)));
        loop {
            if word != 0 {
                return Some(self.base + 64 * w as u32 + 63 - word.leading_zeros());
            }
            if w == 0 {
                return None;
            }
            w -= 1;
            word = self.words[w];
        }
    }

    /// Smallest reached offset `>= offset`.
    pub fn first_reached_at_or_above(&self, offset: u32) -> Option<u32> {
        if self.words.is_empty() || offset >= self.end() {
            return None;
        }
        let i = (offset.max(self.base) - self.base) as usize;
        let mut w = i >> 6;
        let mut word = self.words[w] & (u64::MAX << (i & 63));
        loop {
            if word != 0 {
                return Some(self.base + 64 * w as u32 + word.trailing_zeros());
            }
            w += 1;
            if w == self.words.len() {
                return None;
            }
            word = self.words[w];
        }
    }

    /// Largest reached offset `< offset`.
    #[inline]
    pub fn last_reached_below(&self, offset: u32) -> Option<u32> {
        self.last_reached_at_or_below(offset.checked_sub(1)?)
    }

    /// Smallest reached offset `> offset`.
    #[inline]
    pub fn first_reached_above(&self, offset: u32) -> Option<u32> {
        self.first_reached_at_or_above(offset.checked_add(1)?)
    }

    /// Reached offsets in `lo..=hi`, ascending.
    pub fn reached_in(&self, lo: u32, hi: u32) -> impl Iterator<Item = u32> + '_ {
        std::iter::successors(self.first_reached_at_or_above(lo), move |&o| self.first_reached_above(o))
            .take_while(move |&o| o <= hi)
    }

    fn reset(&mut self) {
        self.words.clear();
        self.base = 0;
        self.n_reached = 0;
    }
}

/// Scores in block windows per group of `B` ranks, and reached bits per
/// rank. A cell is `[u32; S]`, indexed by `AlignState as usize`, so `S` is
/// `1` for linear, `3` for affine and `5` for two-piece costs.
pub struct BlockStorage<N, O, const S: usize, const B: usize = 8> {
    /// By `rank / B`.
    block_rows: Vec<BlockRow<S, B>>,
    /// By rank. Empty until the first bit is set.
    reached: Vec<ReachedBits>,
    /// By node index.
    node_ranks: Vec<u32>,
    /// Block rows with at least one block, for `reset`.
    touched_rows: Vec<u32>,
    /// Ranks with at least one bit set, for `reset`.
    touched_reached: Vec<u32>,
    _marker: PhantomData<(N, O)>,
}

impl<N, O, const S: usize, const B: usize> BlockStorage<N, O, S, B>
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
        let mut node_ranks = vec![0u32; node_count];
        for node in graph.all_nodes_iter() {
            node_ranks[node.index()] = graph.node_rank(node) as u32;
        }
        Self {
            block_rows: std::iter::repeat_with(BlockRow::new).take(node_count.div_ceil(B)).collect(),
            reached: std::iter::repeat_with(ReachedBits::default).take(node_count).collect(),
            node_ranks,
            touched_rows: Vec::new(),
            touched_reached: Vec::new(),
            _marker: PhantomData,
        }
    }

    /// Clear what the previous alignment wrote and take the ranks of
    /// `graph`, which may have changed since. `blocks.clear()` keeps each
    /// row's capacity, which is the point. A row's `base` is stale after the
    /// clear, but `grow_to` resets it on the first touch because the row is
    /// empty. Every row the previous alignment wrote is cleared here, so no
    /// stale score can be read under a new rank.
    pub fn reset<G>(&mut self, graph: &G)
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        for &row in &self.touched_rows {
            self.block_rows[row as usize].blocks.clear();
        }
        self.touched_rows.clear();
        for &rank in &self.touched_reached {
            self.reached[rank as usize].reset();
        }
        self.touched_reached.clear();

        let node_count = graph.node_count();
        self.block_rows.resize_with(node_count.div_ceil(B), BlockRow::new);
        self.reached.resize_with(node_count, ReachedBits::default);
        self.node_ranks.clear();
        self.node_ranks.resize(node_count, 0);
        for node in graph.all_nodes_iter() {
            self.node_ranks[node.index()] = graph.node_rank(node) as u32;
        }
    }

    #[inline(always)]
    fn rank(&self, aln_node: &AlignmentGraphNode<N, O>) -> usize {
        self.node_ranks[aln_node.node().index()] as usize
    }

    #[inline]
    pub fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> u32 {
        let rank = self.rank(aln_node);
        self.block_rows[rank >> B.ilog2()].get(rank & (B - 1), aln_node.offset().as_usize() as u32, aln_state as usize)
    }

    #[inline]
    fn cell_mut(&mut self, aln_node: &AlignmentGraphNode<N, O>) -> &mut [u32; S] {
        let rank = self.rank(aln_node);
        let row_index = rank >> B.ilog2();
        let row = &mut self.block_rows[row_index];
        if row.blocks.is_empty() {
            self.touched_rows.push(row_index as u32);
        }
        row.cell_mut(rank & (B - 1), aln_node.offset().as_usize() as u32)
    }

    #[inline]
    pub fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32) {
        self.cell_mut(aln_node)[aln_state as usize] = score;
    }

    #[inline]
    pub fn update_score_if_lower_block(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: u32) -> bool {
        let cell = &mut self.cell_mut(aln_node)[aln_state as usize];
        if score < *cell {
            *cell = score;
            true
        } else {
            false
        }
    }

    #[inline]
    pub fn mark_reached(&mut self, aln_node: &AlignmentGraphNode<N, O>) {
        let rank = self.rank(aln_node);
        let bits = &mut self.reached[rank];
        if bits.words.is_empty() {
            self.touched_reached.push(rank as u32);
        }
        bits.mark(aln_node.offset().as_usize() as u32);
    }

    #[inline]
    pub fn reached(&self, node: N) -> &ReachedBits {
        &self.reached[self.node_ranks[node.index()] as usize]
    }
}

/// Everything the A* search records per state, for any cost model: scores in
/// block windows per rank group, and per rank the query offsets at which its
/// `Match` state was reached (what bubble-exit pruning checks against).
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
    visited: BlockStorage<N, O, S>,
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
            visited: BlockStorage::new(ref_graph),
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
        if aln_state == AlignState::Match {
            self.visited.mark_reached(aln_node);
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
            let reached = self.visited.reached(bubble.bubble_exit);
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

    fn reset<G>(&mut self, graph: &G, seq: &[u8], bubbles: Arc<BubbleIndex<N>>)
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        self.seq_len = seq.len();
        self.bubble_index = bubbles;
        self.visited.reset(graph);
    }
}

// ─── Backtrace: one-state linear model ───────────────────────────────────────

/// One backtrace step for the linear model. Returns the predecessor and the
/// kind of move that led from it to `curr` (`Match` for the diagonal), or
/// `None` at the start of the alignment.
fn linear_step<N, O, C, G, const S: usize>(
    visited: &BlockStorage<N, O, S>,
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
    visited: &BlockStorage<N, O, S>,
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
    visited: &BlockStorage<N, O, S>,
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
    visited: &BlockStorage<N, O, S>,
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
    use std::collections::BTreeSet;

    use super::*;
    use crate::graph::mock::create_test_graph1;
    use crate::graph::traits::GraphWithNodeOrdering;

    type O = u32;

    #[test]
    fn get_set_round_trip() {
        let g = create_test_graph1();
        let mut storage: BlockStorage<_, O, 3> = BlockStorage::new(&g);

        let node = AlignmentGraphNode::new(g.rank_to_node(0), 5u32);
        let next_rank = AlignmentGraphNode::new(g.rank_to_node(1), 5u32);

        assert_eq!(storage.get_score(&node, AlignState::Match), INF, "unset cell should read INF");

        storage.set_score(&node, AlignState::Match, 7);
        assert_eq!(storage.get_score(&node, AlignState::Match), 7);
        // Other states in the same cell stay untouched.
        assert_eq!(storage.get_score(&node, AlignState::Deletion), INF);
        // The next rank at the same offset shares the block but not the cell.
        assert_eq!(storage.get_score(&next_rank, AlignState::Match), INF);
    }

    #[test]
    fn grows_in_both_directions() {
        let g = create_test_graph1();
        let mut storage: BlockStorage<_, O, 3> = BlockStorage::new(&g);
        let n = g.rank_to_node(2);
        let at = |offset: u32| AlignmentGraphNode::new(n, offset);

        for (offset, score) in [(100, 1), (20, 2), (900, 3)] {
            storage.set_score(&at(offset), AlignState::Match, score);
        }

        assert_eq!(storage.get_score(&at(100), AlignState::Match), 1);
        assert_eq!(storage.get_score(&at(20), AlignState::Match), 2);
        assert_eq!(storage.get_score(&at(900), AlignState::Match), 3);
        for offset in [19, 21, 899] {
            assert_eq!(storage.get_score(&at(offset), AlignState::Match), INF);
        }
    }

    #[test]
    fn update_score_if_lower_only_updates_on_improvement() {
        let g = create_test_graph1();
        let mut storage: BlockStorage<_, O, 3> = BlockStorage::new(&g);

        let node = AlignmentGraphNode::new(g.rank_to_node(2), 3u32);

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
    fn reset_clears_touched_rows_and_bits() {
        let g = create_test_graph1();
        let mut storage: BlockStorage<_, O, 3> = BlockStorage::new(&g);
        let a = AlignmentGraphNode::new(g.rank_to_node(1), 4u32);
        let b = AlignmentGraphNode::new(g.rank_to_node(9), 300u32);

        for node in [&a, &b] {
            storage.set_score(node, AlignState::Match, 3);
            storage.mark_reached(node);
        }
        storage.reset(&g);

        for node in [&a, &b] {
            assert_eq!(storage.get_score(node, AlignState::Match), INF);
            assert!(storage.reached(node.node()).is_empty());
            assert_eq!(storage.reached(node.node()).reached_in(0, u32::MAX).count(), 0);
        }
        assert!(storage.touched_rows.is_empty() && storage.touched_reached.is_empty());
    }

    #[test]
    fn reached_bits_survive_growth() {
        let mut bits = ReachedBits::default();
        for offset in [130, 700, 5] {
            bits.mark(offset);
        }
        assert_eq!(bits.reached_in(0, u32::MAX).collect::<Vec<_>>(), vec![5, 130, 700]);
        assert_eq!(bits.n_reached, 3);

        bits.mark(130);
        assert_eq!(bits.n_reached, 3);
    }

    /// Guards the pruning check's equivalence with the sorted offset vector
    /// it replaced: every scan must agree with a `BTreeSet`.
    #[test]
    fn bit_scans_match_btreeset() {
        let mut x = 42u64;
        let mut next = || {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (x >> 33) as u32 % 2000
        };

        let mut bits = ReachedBits::default();
        let mut set = BTreeSet::new();
        for _ in 0..300 {
            let offset = next();
            bits.mark(offset);
            set.insert(offset);
        }
        assert_eq!(bits.n_reached as usize, set.len());

        for x in 0..2000u32 {
            assert_eq!(bits.last_reached_below(x), set.range(..x).next_back().copied(), "below {x}");
            assert_eq!(bits.first_reached_above(x), set.range(x + 1..).next().copied(), "above {x}");
            assert_eq!(
                bits.reached_in(x, x + 37).collect::<Vec<_>>(),
                set.range(x..=x + 37).copied().collect::<Vec<_>>(),
                "in {x}..={}",
                x + 37
            );
        }
    }
}
