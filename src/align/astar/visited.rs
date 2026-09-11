//! Visited-state storage for the A* engine, ported from `stable`'s
//! `aligner::astar::AstarVisited` + `aligner::scoring::gap_affine`'s
//! `BlockedVisitedStorageAffine`/`AffineAstarData`.

use std::fmt::Write;
use std::marker::PhantomData;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::align::cost_models::AlignmentCostModel;
use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::traits::AlignedPair;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::offsets::OffsetType;
use crate::align::score::INF;
use crate::align::traits::AlignableGraph;
use crate::align::AlignState;
use crate::graph::bubbles::index::BubbleIndex;
use crate::align::astar::pruning::ReachedBubbleExitsMatch;
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

    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), std::fmt::Error>;
}

/// Node-rank × query-offset space tiled into `B×B` blocks, allocated on
/// demand. A cell is `[u32; S]` — the same layout `DPKernel` uses for a band
/// column (slot `i` = `AlignState as u8`), so `S` is `3` for affine, `5` for
/// two-piece, `1` for linear.
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

    pub fn write_tsv<W: Write>(&self, writer: &mut W, state_names: &[&str]) -> Result<(), std::fmt::Error> {
        writeln!(writer, "node_id\toffset\tmatrix\tscore")?;

        let mut rank_to_node_ix: Vec<usize> = (0..self.node_ranks.len()).collect();
        rank_to_node_ix.sort_unstable_by_key(|v| self.node_ranks[*v]);

        for (node_block_num, node_block) in self.node_blocks.iter().enumerate() {
            let node_rank_base = node_block_num * B;
            for (qry_block_num, block_data) in node_block.iter() {
                let qry_pos_base = O::new(qry_block_num.as_usize() * B);
                for (row_num, row) in block_data.iter().enumerate() {
                    let node_rank = node_rank_base + row_num;
                    if node_rank >= rank_to_node_ix.len() {
                        break;
                    }
                    let node_ix = rank_to_node_ix[node_rank];

                    for (col, cell) in row.iter().enumerate() {
                        let qry_pos = qry_pos_base + O::new(col);
                        for (state_ix, &score) in cell.iter().enumerate() {
                            if score < INF {
                                let name = state_names.get(state_ix).copied().unwrap_or("?");
                                writeln!(writer, "{node_ix}\t{qry_pos:?}\t{name}\t{score}")?;
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

/// `AstarVisited` impl for any three-state (Match/Deletion/Insertion) cost
/// model, composing the block storage with bubble-index-based pruning.
pub struct ThreeStateAstarData<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel,
{
    costs: C,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    visited: BlockedVisitedStorage<N, O, 3>,
    bubbles_reached_m: Vec<Vec<O>>,
}

/// Gap-affine visited-state storage (three-state data specialised to [`Affine`]).
pub type AffineAstarData<N, O> = ThreeStateAstarData<N, O, Affine>;

/// Linear-gap visited-state storage (three-state data specialised to [`Linear`]).
pub type LinearAstarData<N, O> = ThreeStateAstarData<N, O, Linear>;

impl<N, O, C> ThreeStateAstarData<N, O, C>
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

    /// Walk one backtrace step: given the current alignment-graph node and
    /// state, find the predecessor `(node, state)` pair whose score is
    /// consistent with the recorded transition, or `None` at the start of
    /// the alignment.
    fn get_backtrace<G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> Option<(AlignmentGraphNode<N, O>, AlignState)>
    where
        G: AlignableGraph<Node = N>,
    {
        let curr_score = self.visited.get_score(aln_node, aln_state);

        match aln_state {
            AlignState::Match => {
                if aln_node.offset() > O::zero() {
                    let is_match_or_end = ref_graph.is_symbol_equal(aln_node.node(), seq[aln_node.offset().as_usize() - 1])
                        || aln_node.node() == ref_graph.end_node();
                    let pred_offset = if aln_node.node() == ref_graph.end_node() {
                        aln_node.offset()
                    } else {
                        aln_node.offset() - O::one()
                    };

                    for p in ref_graph.predecessors(aln_node.node()) {
                        let pred = AlignmentGraphNode::new(p, pred_offset);
                        let pred_score = self.visited.get_score(&pred, AlignState::Match);

                        if (is_match_or_end && pred_score == curr_score)
                            || (!is_match_or_end && pred_score == curr_score - self.costs.mismatch() as u32)
                        {
                            return Some((pred, AlignState::Match));
                        }
                    }
                }

                let pred_score = self.visited.get_score(aln_node, AlignState::Deletion);
                if pred_score == curr_score {
                    return Some((*aln_node, AlignState::Deletion));
                }

                let pred_score = self.visited.get_score(aln_node, AlignState::Insertion);
                if pred_score == curr_score {
                    return Some((*aln_node, AlignState::Insertion));
                }
            }
            AlignState::Deletion => {
                let go_ge = (self.costs.gap_open() + self.costs.gap_extend()) as u32;
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    let pred_score = self.visited.get_score(&pred, AlignState::Match);
                    if pred_score == curr_score - go_ge {
                        return Some((pred, AlignState::Match));
                    }
                }

                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    let pred_score = self.visited.get_score(&pred, AlignState::Deletion);
                    if pred_score == curr_score - self.costs.gap_extend() as u32 {
                        return Some((pred, AlignState::Deletion));
                    }
                }
            }
            AlignState::Insertion => {
                if aln_node.offset() > O::zero() {
                    let go_ge = (self.costs.gap_open() + self.costs.gap_extend()) as u32;
                    let pred = AlignmentGraphNode::new(aln_node.node(), aln_node.offset() - O::one());
                    let pred_score = self.visited.get_score(&pred, AlignState::Match);
                    if pred_score == curr_score - go_ge {
                        return Some((pred, AlignState::Match));
                    }

                    let pred_score = self.visited.get_score(&pred, AlignState::Insertion);
                    if pred_score == curr_score - self.costs.gap_extend() as u32 {
                        return Some((pred, AlignState::Insertion));
                    }
                }
            }
            AlignState::Insertion2 | AlignState::Deletion2 => {
                panic!("Invalid gap-affine state {aln_state:?}!")
            }
        }

        None
    }
}

impl<N, O, C> AstarVisited<N, O> for ThreeStateAstarData<N, O, C>
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
        let Some((mut curr, mut curr_state)) = self.get_backtrace(ref_graph, seq, aln_node, AlignState::Match) else {
            panic!("No backtrace for alignment end state?");
        };

        let mut alignment = Vec::new();

        // `aln_node` (the alignment end state) is itself never pushed below --
        // the loop only records `curr` once its *predecessor* has been found.
        // For `Global` mode that's fine: the end state is always the virtual
        // end sentinel, which has no base to record. But `EndsFree`/semi-global
        // alignments can terminate on a real graph node (see `ends_free_is_end`),
        // whose own base would otherwise be silently dropped. Record it here,
        // unless it's just an alias for a just-closed indel at the same
        // (node, offset) -- that case is handled once the loop reaches it below.
        let is_indel_alias = curr.node() == aln_node.node() && curr.offset() == aln_node.offset();
        if aln_node.node() != ref_graph.end_node() && !is_indel_alias {
            alignment.push(AlignedPair::new(Some(aln_node.node()), Some(aln_node.offset().as_usize() - 1)));
        }

        while let Some((bt_node, bt_state)) = self.get_backtrace(ref_graph, seq, &curr, curr_state) {
            // If BT points towards indel, update the backtrace again to prevent double
            // using (node, query) pairs, since closing of indels is a zero cost edge.
            if curr_state == AlignState::Match && (bt_state == AlignState::Insertion || bt_state == AlignState::Deletion) {
                curr = bt_node;
                curr_state = bt_state;
                continue;
            }

            match curr_state {
                AlignState::Match => {
                    alignment.push(AlignedPair::new(Some(curr.node()), Some(curr.offset().as_usize() - 1)));
                }
                AlignState::Insertion => {
                    alignment.push(AlignedPair::new(None, Some(curr.offset().as_usize() - 1)));
                }
                AlignState::Deletion => {
                    alignment.push(AlignedPair::new(Some(curr.node()), None));
                }
                AlignState::Insertion2 | AlignState::Deletion2 => {
                    panic!("Unexpected align state {curr_state:?} in backtrace!")
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

    #[inline]
    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), std::fmt::Error> {
        self.visited.write_tsv(writer, &["M", "D", "I"])
    }
}

/// `AstarVisited` impl for the two-piece gap-affine cost model: same shape as
/// [`ThreeStateAstarData`] but with 5 states (`Match`/`Deletion`/`Insertion`/
/// `Deletion2`/`Insertion2`) and backtrace logic that checks both pieces.
pub struct TwoPieceAstarData<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    costs: TwoPieceAffine,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    visited: BlockedVisitedStorage<N, O, 5>,
    bubbles_reached_m: Vec<Vec<O>>,
}

impl<N, O> TwoPieceAstarData<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    pub fn new<G>(costs: TwoPieceAffine, ref_graph: &G, seq: &[u8], bubble_index: Arc<BubbleIndex<N>>) -> Self
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

    fn get_backtrace<G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> Option<(AlignmentGraphNode<N, O>, AlignState)>
    where
        G: AlignableGraph<Node = N>,
    {
        let curr_score = self.visited.get_score(aln_node, aln_state);

        match aln_state {
            AlignState::Match => {
                if aln_node.offset() > O::zero() {
                    let is_match_or_end = ref_graph.is_symbol_equal(aln_node.node(), seq[aln_node.offset().as_usize() - 1])
                        || aln_node.node() == ref_graph.end_node();
                    let pred_offset = if aln_node.node() == ref_graph.end_node() {
                        aln_node.offset()
                    } else {
                        aln_node.offset() - O::one()
                    };

                    for p in ref_graph.predecessors(aln_node.node()) {
                        let pred = AlignmentGraphNode::new(p, pred_offset);
                        let pred_score = self.visited.get_score(&pred, AlignState::Match);

                        if (is_match_or_end && pred_score == curr_score)
                            || (!is_match_or_end && pred_score == curr_score - self.costs.mismatch() as u32)
                        {
                            return Some((pred, AlignState::Match));
                        }
                    }
                }

                // Close deletion (either piece)
                if self.visited.get_score(aln_node, AlignState::Deletion) == curr_score {
                    return Some((*aln_node, AlignState::Deletion));
                }
                if self.visited.get_score(aln_node, AlignState::Deletion2) == curr_score {
                    return Some((*aln_node, AlignState::Deletion2));
                }

                // Close insertion (either piece)
                if self.visited.get_score(aln_node, AlignState::Insertion) == curr_score {
                    return Some((*aln_node, AlignState::Insertion));
                }
                if self.visited.get_score(aln_node, AlignState::Insertion2) == curr_score {
                    return Some((*aln_node, AlignState::Insertion2));
                }
            }
            AlignState::Deletion => {
                let go1_ge1 = (self.costs.gap_open() + self.costs.gap_extend()) as u32;
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    if self.visited.get_score(&pred, AlignState::Match) == curr_score - go1_ge1 {
                        return Some((pred, AlignState::Match));
                    }
                }
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    if self.visited.get_score(&pred, AlignState::Deletion) == curr_score - self.costs.gap_extend() as u32 {
                        return Some((pred, AlignState::Deletion));
                    }
                }
            }
            AlignState::Deletion2 => {
                let go2_ge2 = (self.costs.gap_open2() + self.costs.gap_extend2()) as u32;
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    if self.visited.get_score(&pred, AlignState::Match) == curr_score - go2_ge2 {
                        return Some((pred, AlignState::Match));
                    }
                }
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    if self.visited.get_score(&pred, AlignState::Deletion2) == curr_score - self.costs.gap_extend2() as u32 {
                        return Some((pred, AlignState::Deletion2));
                    }
                }
            }
            AlignState::Insertion => {
                if aln_node.offset() > O::zero() {
                    let go1_ge1 = (self.costs.gap_open() + self.costs.gap_extend()) as u32;
                    let pred = AlignmentGraphNode::new(aln_node.node(), aln_node.offset() - O::one());
                    if self.visited.get_score(&pred, AlignState::Match) == curr_score - go1_ge1 {
                        return Some((pred, AlignState::Match));
                    }
                    if self.visited.get_score(&pred, AlignState::Insertion) == curr_score - self.costs.gap_extend() as u32 {
                        return Some((pred, AlignState::Insertion));
                    }
                }
            }
            AlignState::Insertion2 => {
                if aln_node.offset() > O::zero() {
                    let go2_ge2 = (self.costs.gap_open2() + self.costs.gap_extend2()) as u32;
                    let pred = AlignmentGraphNode::new(aln_node.node(), aln_node.offset() - O::one());
                    if self.visited.get_score(&pred, AlignState::Match) == curr_score - go2_ge2 {
                        return Some((pred, AlignState::Match));
                    }
                    if self.visited.get_score(&pred, AlignState::Insertion2) == curr_score - self.costs.gap_extend2() as u32 {
                        return Some((pred, AlignState::Insertion2));
                    }
                }
            }
        }

        None
    }
}

impl<N, O> AstarVisited<N, O> for TwoPieceAstarData<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
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
        let Some((mut curr, mut curr_state)) = self.get_backtrace(ref_graph, seq, aln_node, AlignState::Match) else {
            panic!("No backtrace for alignment end state?");
        };

        let mut alignment = Vec::new();

        // See the analogous comment in `ThreeStateAstarData::backtrace`: the
        // end state itself is never pushed by the loop below, which is only
        // correct when it's the virtual end sentinel (`Global` mode). For
        // `EndsFree`/semi-global alignments ending on a real node, record it
        // here unless it's an alias for a just-closed indel at the same
        // (node, offset).
        let is_indel_alias = curr.node() == aln_node.node() && curr.offset() == aln_node.offset();
        if aln_node.node() != ref_graph.end_node() && !is_indel_alias {
            alignment.push(AlignedPair::new(Some(aln_node.node()), Some(aln_node.offset().as_usize() - 1)));
        }

        while let Some((bt_node, bt_state)) = self.get_backtrace(ref_graph, seq, &curr, curr_state) {
            if curr_state == AlignState::Match
                && matches!(
                    bt_state,
                    AlignState::Insertion | AlignState::Deletion | AlignState::Insertion2 | AlignState::Deletion2
                )
            {
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

    #[inline]
    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), std::fmt::Error> {
        self.visited.write_tsv(writer, &["M", "D1", "I1", "D2", "I2"])
    }
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
