use std::fmt;
use tracing::{debug, debug_span, span, trace, Level};

use super::AlignmentCostModel;
use crate::aligner::{
    astar::{
        queue::{LayeredQueue, QueueLayer},
        AlignState, AstarState,
    },
    extension::extend,
    fr_points::{to_node_pos, Diag, DiagType, NodeFrPoints, PosType, Score},
    utils::AlignedPair,
    AlignmentMode,
};
use crate::aligner::traits::{AlignableGraph, AlignableGraphNodePos};
use crate::graph::traits::GraphNodeId;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Affine {
    mismatch: u8,
    gap_open: u8,
    gap_extend: u8,
}

impl Affine {
    pub fn new(mismatch: u8, gap_open: u8, gap_extend: u8) -> Self {
        Self {
            mismatch,
            gap_open,
            gap_extend,
        }
    }
}

impl AlignmentCostModel for Affine {
    type AstarStateType<G, D, O> = AffineAstarState<G, D, O>
        where
            G: AlignableGraph,
            D: DiagType,
            O: PosType;

    fn initialize<G, D, O, F>(
        &self,
        graph: &G,
        seq: &[u8],
        alignment_mode: AlignmentMode,
        heuristic: F,
    ) -> Self::AstarStateType<G, D, O>
    where
        G: AlignableGraph,
        D: DiagType,
        O: PosType,
        F: Fn(&<Self::AstarStateType<G, D, O> as AstarState<G>>::AstarItem) -> usize,
    {
        let mut state = AffineAstarState::new(*self, graph, seq, alignment_mode);
        let initial_state = AffineAstarItem::new(
            Score::default(),
            graph.start_node(),
            Diag::default(),
            AlignState::Match,
        );
        
        state.update_if_further(&initial_state, 0);
        let h = heuristic(&initial_state);
        state.queue_item(initial_state, h);

        state
    }

    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.mismatch
    }

    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.gap_open
    }

    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.gap_extend
    }

    #[inline(always)]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        0
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AffineAstarItem<N, D> {
    pub score: Score,
    pub node: N,
    pub diag: Diag<D>,
    pub state: AlignState,
}

impl<N, D> AffineAstarItem<N, D> {
    pub fn new(score: Score, node: N, diag: Diag<D>, state: AlignState) -> Self {
        Self {
            score,
            node,
            diag,
            state,
        }
    }
}

impl<N, D> From<(Score, N, Diag<D>)> for AffineAstarItem<N, D> {
    fn from((score, node, diag): (Score, N, Diag<D>)) -> Self {
        Self::new(score, node, diag, AlignState::Match)
    }
}

impl<N, D> From<(Score, N, Diag<D>, AlignState)> for AffineAstarItem<N, D> {
    fn from((score, node, diag, state): (Score, N, Diag<D>, AlignState)) -> Self {
        Self::new(score, node, diag, state)
    }
}

#[derive(Clone, Debug, Default)]
struct AffineNodeDiagonals<D, O> {
    fr_points_m: NodeFrPoints<D, O>,
    fr_points_i: NodeFrPoints<D, O>,
    fr_points_d: NodeFrPoints<D, O>,
}

impl<D, O> AffineNodeDiagonals<D, O>
where
    D: DiagType,
    O: PosType,
{
    fn get_furthest<N>(&self, item: &AffineAstarItem<N, D>) -> Option<O> {
        match item.state {
            AlignState::Match => self.fr_points_m.get_furthest(item.score, item.diag),
            AlignState::Deletion => self.fr_points_d.get_furthest(item.score, item.diag),
            AlignState::Insertion => self.fr_points_i.get_furthest(item.score, item.diag),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }
    }

    fn is_further<N>(&self, item: &AffineAstarItem<N, D>, offset: O) -> bool {
        match item.state {
            AlignState::Match => self.fr_points_m.is_further(item.score, item.diag, offset),
            AlignState::Deletion => self.fr_points_d.is_further(item.score, item.diag, offset),
            AlignState::Insertion => self.fr_points_i.is_further(item.score, item.diag, offset),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }
    }

    fn update_if_further<N>(&mut self, item: &AffineAstarItem<N, D>, offset: O) -> bool {
        match item.state {
            AlignState::Match => self
                .fr_points_m
                .update_if_further(item.score, item.diag, offset),
            AlignState::Deletion => self
                .fr_points_d
                .update_if_further(item.score, item.diag, offset),
            AlignState::Insertion => self
                .fr_points_i
                .update_if_further(item.score, item.diag, offset),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }
    }

    fn is_visited<N>(&self, item: &AffineAstarItem<N, D>) -> bool {
        match item.state {
            AlignState::Match => self.fr_points_m.is_visited(item.score, item.diag),
            AlignState::Deletion => self.fr_points_d.is_visited(item.score, item.diag),
            AlignState::Insertion => self.fr_points_i.is_visited(item.score, item.diag),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }
    }

    fn set_visited<N>(&mut self, item: &AffineAstarItem<N, D>) {
        match item.state {
            AlignState::Match => self.fr_points_m.set_visited(item.score, item.diag, true),
            AlignState::Deletion => self.fr_points_d.set_visited(item.score, item.diag, true),
            AlignState::Insertion => self.fr_points_i.set_visited(item.score, item.diag, true),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }
    }
}

/// Holds all A* state information (reached points, queued points, ...) for
/// the gap-affine alignment model
pub struct AffineAstarState<G, D, O>
where
    G: AlignableGraph,
{
    /// Alignment costs
    costs: Affine,

    /// The length of the sequence we are aligning
    seq_length: usize,

    /// The kind of alignment we are performing
    alignment_mode: AlignmentMode,

    /// For each node in the POA graph, we store the maximum reached query offsets per diagonal
    fr_points: Vec<AffineNodeDiagonals<D, O>>,

    /// A* queue
    queue: LayeredQueue<AffineQueueLayer<G::NodeType, D>>,
}

impl<G, D, O> AffineAstarState<G, D, O>
where
    G: AlignableGraph,
    D: DiagType,
    O: PosType,
{
    fn new(costs: Affine, graph: &G, seq: &[u8], alignment_mode: AlignmentMode) -> Self {
        let node_count = graph.node_count();
        Self {
            costs,
            seq_length: seq.len(),
            alignment_mode,
            fr_points: vec![AffineNodeDiagonals::default(); node_count],
            queue: LayeredQueue::new(),
        }
    }

    fn relax_match<F>(
        &mut self,
        graph: &G,
        seq: &[u8],
        item: &AffineAstarItem<G::NodeType, D>,
        heuristic: F,
    ) where
        F: Fn(&AffineAstarItem<G::NodeType, D>) -> usize,
    {
        if item.node == graph.end_node() {
            return;
        }
        
        let fr_point = self.fr_points[item.node.index()]
            .get_furthest(item)
            .unwrap_or(O::zero());
        
        let node_len = graph.node_len(item.node);
        let node_pos = to_node_pos(item.diag, fr_point.as_usize());
        
        debug!("Relaxing match, {:?}, length: {node_len}, pos: {node_pos}", item.node);

        if node_pos == node_len - 1 {
            self.relax_match_at_succ(graph, seq, item, heuristic, fr_point);
        } else {
            // Not at the node end yet, so queue mismatch, deletion, and insertion states within the node
            if fr_point.as_usize() < self.seq_length {
                let new_item_mis = AffineAstarItem::new(
                    item.score + self.costs.mismatch(),
                    item.node,
                    item.diag,
                    AlignState::Match,
                );
                let extended_offset = extend(graph, seq, item.node, item.diag, fr_point.increase_one());
                
                trace!(
                    target: "poasta::aligner::cost_models::affine::extend",
                    score=new_item_mis.score.as_usize(),
                    node=new_item_mis.node.index(),
                    diag=new_item_mis.diag.as_isize(),
                    offset=fr_point.increase_one().as_usize(),
                    extended_offset=extended_offset.as_usize(),
                    state=tracing::field::debug(&item.state),
                    "extend"
                ); 
                
                if self.update_if_further(&new_item_mis, extended_offset.as_usize()) {
                    let h = heuristic(&new_item_mis);
                    self.queue_item(new_item_mis, h);
                }
                
                let new_item_ins = AffineAstarItem::new(
                    item.score + self.costs.gap_open() + self.costs.gap_extend(),
                    item.node,
                    item.diag + 1isize,
                    AlignState::Insertion,
                );
                
                if self.update_if_further(&new_item_ins, fr_point.increase_one().as_usize()) {
                    let h = heuristic(&new_item_ins);
                    self.queue_item(new_item_ins, h);
                }
            }
            
            let new_item_del = AffineAstarItem::new(
                item.score + self.costs.gap_open() + self.costs.gap_extend(),
                item.node,
                item.diag - 1isize,
                AlignState::Deletion,
            );
            
            if self.update_if_further(&new_item_del, fr_point.as_usize()) {
                let h = heuristic(&new_item_del);
                self.queue_item(new_item_del, h);
            }
        }
    }
    
    fn relax_match_at_succ(
        &mut self,
        graph: &G,
        seq: &[u8],
        item: &AffineAstarItem<G::NodeType, D>,
        heuristic: impl Fn(&AffineAstarItem<G::NodeType, D>) -> usize,
        curr_offset: O,
    ) {
        // At the end of a node, we need to check successors to extend the (mis)match
        let node_len = graph.node_len(item.node);
        let mut any_mismatch = false;
        
        debug!("At node end, checking successors.");
        
        for succ in graph.successors(item.node) {
            let new_diag = item.diag + node_len;
            let succ_offset = curr_offset.increase_one();
            let succ_node_pos = to_node_pos(new_diag, succ_offset.as_usize());
            let succ_seq = graph.node_seq(succ);
            
            debug!("- Successor: {succ:?}, new diagonal: {new_diag:?}, new offset: {succ_offset:?}");
            
            if succ == graph.end_node() {
                debug!("Succ is end node");

                // We reached the end node, immediately queue with zero cost
                let new_item =
                    AffineAstarItem::new(item.score, succ, new_diag, AlignState::Match);

                if self
                    .update_if_further(&new_item, curr_offset.increase_one().as_usize())
                {
                    let h = heuristic(&new_item);
                    self.queue_item(new_item, h);
                }
                
                any_mismatch = true;

                continue;
            }
            
            if succ_offset.as_usize() > self.seq_length {
                debug!("Moving beyond sequence, opening deletion only.");
                
                // Beyond the sequence length, but ensure we open deletion on the successor node
                let new_item_del = AffineAstarItem::new(
                    item.score + self.costs.gap_open() + self.costs.gap_extend(),
                    succ,
                    new_diag - 1isize,
                    AlignState::Deletion,
                );
                
                if self.update_if_further(&new_item_del, curr_offset.as_usize()) {
                    let h = heuristic(&new_item_del);
                    self.queue_item(new_item_del, h);
                }
            } else if succ_seq[succ_node_pos] != seq[succ_offset.as_usize()-1] {
                debug!("Mismatch at successor node {succ:?}, pos: {succ_node_pos}, qry offset: {:?}", 
                    succ_offset);
                
                let new_item_mis = AffineAstarItem::new(
                    item.score + self.costs.mismatch(),
                    succ,
                    new_diag,
                    AlignState::Match,
                );
                
                let extended_offset = extend(graph, seq, succ, new_diag, succ_offset);
                trace!(
                    target: "poasta::aligner::cost_models::affine::extend",
                    score=new_item_mis.score.as_usize(),
                    node=new_item_mis.node.index(),
                    diag=new_item_mis.diag.as_isize(),
                    offset=succ_offset.as_usize(),
                    extended_offset=extended_offset.as_usize(),
                    state=tracing::field::debug(&new_item_mis.state),
                    "extend"
                ); 
                
                if self.update_if_further(&new_item_mis, extended_offset.as_usize()) {
                    let h = heuristic(&new_item_mis);
                    self.queue_item(new_item_mis, h);
                }
                
                let new_item_del = AffineAstarItem::new(
                    item.score + self.costs.gap_open() + self.costs.gap_extend(),
                    succ,
                    new_diag - 1isize,
                    AlignState::Deletion,
                );
                
                if self.update_if_further(&new_item_del, curr_offset.as_usize()) {
                    let h = heuristic(&new_item_del);
                    self.queue_item(new_item_del, h);
                }
                
                any_mismatch = true;
            } else {
                debug!("Match at successor node {succ:?}, pos: {succ_node_pos}, qry offset: {:?}", 
                    succ_offset);
                let extended_offset = extend(graph, seq, succ, new_diag, succ_offset);
                
                let new_item_match = AffineAstarItem::new(
                    item.score,
                    succ,
                    new_diag,
                    AlignState::Match,
                );
                
                trace!(
                    target: "poasta::aligner::cost_models::affine::extend",
                    score=new_item_match.score.as_usize(),
                    node=new_item_match.node.index(),
                    diag=new_item_match.diag.as_isize(),
                    offset=succ_offset.as_usize(),
                    extended_offset=extended_offset.as_usize(),
                    state=tracing::field::debug(&new_item_match.state),
                    "extend"
                ); 

                if self.update_if_further(&new_item_match, extended_offset.as_usize()) {
                    let h = heuristic(&new_item_match);
                    self.queue_item(new_item_match, h);
                }
            }
        }
        
        if any_mismatch {
            // If we had a mismatch, we should also open an insertion on the current node, one diagonal up
            if curr_offset.as_usize() < self.seq_length {
                let new_item_ins = AffineAstarItem::new(
                    item.score + self.costs.gap_open() + self.costs.gap_extend(),
                    item.node,
                    item.diag + 1isize,
                    AlignState::Insertion,
                );

                if self.update_if_further(&new_item_ins, curr_offset.increase_one().as_usize()) {
                    let h = heuristic(&new_item_ins);
                    self.queue_item(new_item_ins, h);
                }
            }
        }
    }
    
    fn relax_deletion<F>(
        &mut self,
        graph: &G,
        seq: &[u8],
        item: &AffineAstarItem<G::NodeType, D>,
        heuristic: F,
    ) where
        F: Fn(&AffineAstarItem<G::NodeType, D>) -> usize,
    {
        let fr_point = self.fr_points[item.node.index()]
            .get_furthest(item)
            .unwrap_or(O::zero());

        let node_len = graph.node_len(item.node);
        let node_pos = to_node_pos(item.diag, fr_point.as_usize());

        if node_pos == node_len - 1 {
            // At the end of a node, we need to check successors to extend the deletion
            for succ in graph.successors(item.node) {
                if succ == graph.end_node() {
                    continue;
                }
                
                let new_item_del = AffineAstarItem::new(
                    item.score + self.costs.gap_extend(),
                    succ,
                    item.diag + node_len - 1isize,
                    AlignState::Deletion,
                );

                if self.update_if_further(&new_item_del, fr_point.as_usize()) {
                    let h = heuristic(&new_item_del);
                    self.queue_item(new_item_del, h)
                }
            }
        } else {
            let new_item_del = AffineAstarItem::new(
                item.score + self.costs.gap_extend(),
                item.node,
                item.diag - 1isize,
                AlignState::Deletion,
            );

            if self.update_if_further(&new_item_del, fr_point.as_usize()) {
                let h = heuristic(&new_item_del);
                self.queue_item(new_item_del, h)
            }
        }

        // Close deletion if we can
        let new_item_match =
            AffineAstarItem::new(item.score, item.node, item.diag, AlignState::Match);
        let extended_offset = extend(graph, seq, item.node, item.diag, fr_point);
        
        trace!(
            target: "poasta::aligner::cost_models::affine::extend",
            score=new_item_match.score.as_usize(),
            node=new_item_match.node.index(),
            diag=new_item_match.diag.as_isize(),
            offset=fr_point.as_usize(),
            extended_offset=extended_offset.as_usize(),
            state=tracing::field::debug(&new_item_match.state),
            "extend"
        ); 
        
        if self.update_if_further(&new_item_match, extended_offset.as_usize()) {
            let h = heuristic(&new_item_match);
            self.queue_item(new_item_match, h)
        }
    }

    fn relax_insertion<F>(
        &mut self,
        graph: &G,
        seq: &[u8],
        item: &AffineAstarItem<G::NodeType, D>,
        heuristic: F,
    ) where
        F: Fn(&AffineAstarItem<G::NodeType, D>) -> usize,
    {
        let fr_point = self.fr_points[item.node.index()]
            .get_furthest(item)
            .unwrap_or(O::zero());

        if fr_point.as_usize() < self.seq_length {
            let new_item_ins = AffineAstarItem::new(
                item.score + self.costs.gap_extend(),
                item.node,
                item.diag + 1isize,
                AlignState::Insertion,
            );

            if self.update_if_further(&new_item_ins, fr_point.increase_one().as_usize()) {
                let h = heuristic(&new_item_ins);
                self.queue_item(new_item_ins, h)
            }
        }

        // Close insertion if we can
        let new_item_match =
            AffineAstarItem::new(item.score, item.node, item.diag, AlignState::Match);
        let extended_offset = extend(graph, seq, item.node, item.diag, fr_point);
        
        trace!(
            target: "poasta::aligner::cost_models::affine::extend",
            score=new_item_match.score.as_usize(),
            node=new_item_match.node.index(),
            diag=new_item_match.diag.as_isize(),
            offset=fr_point.as_usize(),
            extended_offset=extended_offset.as_usize(),
            state=tracing::field::debug(&new_item_match.state),
            "extend"
        ); 

        if self.update_if_further(&new_item_match, extended_offset.as_usize()) {
            let h = heuristic(&new_item_match);
            self.queue_item(new_item_match, h)
        }
    }

    fn get_prev(
        &self,
        graph: &G,
        item: &AffineAstarItem<G::NodeType, D>,
        curr_offset: O,
    ) -> Option<(AffineAstarItem<G::NodeType, D>, O)> {
        let span = debug_span!("get_prev");
        let _enter = span.enter();
        let mut sources = Vec::new();
        let s_gap_open = Score::new(self.costs.gap_open() as u32 + self.costs.gap_extend() as u32);
        let s_gap_extend = Score::new(self.costs.gap_extend() as u32);
        let s_mism = Score::new(self.costs.mismatch() as u32);

        match item.state {
            AlignState::Match => {
                // Check for closed insertions
                sources.push(AffineAstarItem::new(
                    item.score,
                    item.node,
                    item.diag,
                    AlignState::Insertion,
                ));

                // Check for closed insertions
                sources.push(AffineAstarItem::new(
                    item.score,
                    item.node,
                    item.diag,
                    AlignState::Deletion,
                ));

                // Check for deletions that could have been the source, both on the current node as well as predecessors
                let mut possible_pred = vec![item.node];
                possible_pred.extend(graph.predecessors(item.node));

                // Check for mismatches, same diagonal, but potentially on a predecessor node
                if item.score >= s_mism {
                    for pred in &possible_pred {
                        let pred_len = if *pred == item.node {
                            0
                        } else {
                            graph.node_len(*pred)
                        };
                        let pred_diag = item.diag - pred_len as isize;

                        sources.push(AffineAstarItem::new(
                            item.score - s_mism,
                            *pred,
                            pred_diag,
                            AlignState::Match,
                        ));
                    }
                }

                // Check for match on a predecessor node, same score
                for pred in graph.predecessors(item.node) {
                    let pred_len = graph.node_len(pred);
                    let pred_diag = item.diag - pred_len as isize;

                    sources.push(AffineAstarItem::new(
                        item.score,
                        pred,
                        pred_diag,
                        AlignState::Match,
                    ));
                }
            }
            AlignState::Deletion => {
                let mut possible_pred = vec![item.node];
                possible_pred.extend(graph.predecessors(item.node));

                // Check for extended deletions
                if item.score >= s_gap_extend {
                    for pred in &possible_pred {
                        let pred_len = if *pred == item.node {
                            0
                        } else {
                            graph.node_len(*pred)
                        };
                        let pred_diag = item.diag - pred_len as isize + 1isize;

                        sources.push(AffineAstarItem::new(
                            item.score - s_gap_extend,
                            *pred,
                            pred_diag,
                            AlignState::Deletion,
                        ));
                    }
                }

                // Check for opened deletions
                if item.score >= s_gap_open {
                    for pred in &possible_pred {
                        let pred_len = if *pred == item.node {
                            0
                        } else {
                            graph.node_len(*pred)
                        };
                        let pred_diag = item.diag - pred_len as isize + 1isize;

                        sources.push(AffineAstarItem::new(
                            item.score - s_gap_open,
                            *pred,
                            pred_diag,
                            AlignState::Match,
                        ));
                    }
                }
            }
            AlignState::Insertion => {
                if curr_offset.value() > O::zero() {
                    let pred_diag = item.diag - 1isize;
                    if item.score >= s_gap_extend {
                        sources.push(AffineAstarItem::new(
                            item.score - s_gap_extend,
                            item.node,
                            pred_diag,
                            AlignState::Insertion,
                        ));
                    }

                    if item.score >= s_gap_open {
                        sources.push(AffineAstarItem::new(
                            item.score - s_gap_open,
                            item.node,
                            pred_diag,
                            AlignState::Match,
                        ));
                    }
                }
            }
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }

        // Iterator::max_by_key will return the last maximum element if multiple have the same offset,
        // so we ordered the above code to prioritize matches > mismatches > deletions > insertions.
        sources
            .into_iter()
            .filter(|s| self.fr_points[s.node.index()].is_visited(s))
            .filter_map(|s| {
                self.fr_points[s.node.index()]
                    .get_furthest(&s)
                    .and_then(|v| if v.value() <= curr_offset.value() { Some((s, v)) } else { None })
            })
            .inspect(|(item, offset)| {
                debug!(
                    "Checking source: {:?} [{:?}], offset: {:?}, diag: {:?}, node: {:?}",
                    item.state, item.score, offset.as_usize(), item.diag, item.node
                )
            })
            .max_by_key(|(_, offset)| offset.as_usize())
    }
}

impl<G, D, O> AstarState<G> for AffineAstarState<G, D, O>
where
    G: AlignableGraph,
    D: DiagType,
    O: PosType,
{
    type AstarItem = AffineAstarItem<G::NodeType, D>;

    fn pop_front(&mut self) -> Option<Self::AstarItem> {
        self.queue.pop()
    }

    fn is_further(&self, item: &Self::AstarItem, offset: usize) -> bool {
        self.fr_points[item.node.index()].is_further(item, O::new(offset))
    }

    fn is_end(&self, graph: &G, item: &Self::AstarItem) -> bool {
        if item.node != graph.end_node() {
            return false;
        }

        let end_diag = Diag::new(self.seq_length as isize + 1);
        debug!("End diag: {end_diag:?}, item.diag: {:?}", item.diag);
        if item.diag != end_diag {
            return false;
        }

        debug!(
            "Offset {:?}, required seq length: {:?}",
            self.get_offset(item),
            self.seq_length + 1
        );
        self.get_offset(item) == self.seq_length + 1
    }

    fn get_score(&self, item: &Self::AstarItem) -> Score {
        item.score
    }

    fn get_offset(&self, item: &Self::AstarItem) -> usize {
        self.fr_points[item.node.index()]
            .get_furthest(item)
            .map(|v| v.as_usize())
            .unwrap_or(0)
    }

    fn is_visited(&self, item: &Self::AstarItem) -> bool {
        self.fr_points[item.node.index()].is_visited(item)
    }

    fn set_visited(&mut self, item: &Self::AstarItem) {
        let offset = self.get_offset(item);
        trace!(
            target: "poasta::aligner::cost_models::affine::set_visited",
            score=item.score.as_usize(),
            node=item.node.index(),
            diag=item.diag.as_isize(),
            offset=offset,
            state=tracing::field::debug(&item.state)
        );
        self.fr_points[item.node.index()].set_visited(item)
    }

    fn update_if_further(&mut self, item: &Self::AstarItem, offset: usize) -> bool {
        self.fr_points[item.node.index()].update_if_further(item, O::new(offset))
    }

    fn queue_item(&mut self, item: Self::AstarItem, heuristic: usize) {
        let priority = item.score.as_usize() + heuristic;
        let offset = self.get_offset(&item);

        trace!(
            target: "poasta::aligner::cost_models::affine::queue_item",
            score=item.score.as_usize(),
            node=item.node.index(),
            diag=item.diag.as_isize(),
            offset=offset,
            priority=priority,
            state=tracing::field::debug(&item.state),
            "queue item",
        );

        self.queue.queue(item, priority);
    }

    fn relax<F>(&mut self, graph: &G, seq: &[u8], item: &Self::AstarItem, heuristic: F)
    where
        F: Fn(&Self::AstarItem) -> usize,
    {
        match item.state {
            AlignState::Match => self.relax_match(graph, seq, item, heuristic),
            AlignState::Deletion => self.relax_deletion(graph, seq, item, heuristic),
            AlignState::Insertion => self.relax_insertion(graph, seq, item, heuristic),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Unexpected state: {:?}", item.state)
            }
        }
    }

    fn backtrace(&self, graph: &G, end: &Self::AstarItem) -> Vec<AlignedPair<G::NodePosType>> {
        let span = span!(Level::INFO, "backtrace");
        let _enter = span.enter();

        let mut curr = end.clone();
        let mut curr_offset = O::new(self.seq_length + 1);
        curr_offset.set_visited(true);
        let mut alignment = Vec::new();

        while let Some((prev, prev_offset)) = self.get_prev(graph, &curr, curr_offset) {
            debug!(curr = ?curr, curr_offset = ?curr_offset);
            debug!(prev = ?prev, prev_offset = ?prev_offset);

            if curr.node == graph.end_node() {
                curr = prev;
                curr_offset = prev_offset;
                continue;
            }

            match curr.state {
                AlignState::Match => {
                    for query_offset in (prev_offset.as_usize() + 1..=curr_offset.as_usize()).rev()
                    {
                        let node_pos = to_node_pos(curr.diag, query_offset);

                        alignment.push(AlignedPair::new(
                            Some(<G::NodePosType as AlignableGraphNodePos>::new(
                                curr.node, node_pos,
                            )),
                            Some(query_offset - 1),
                        ));
                    }
                }
                AlignState::Deletion => {
                    let node_pos = to_node_pos(curr.diag, curr_offset.as_usize());
                    alignment.push(AlignedPair::new(
                        Some(<G::NodePosType as AlignableGraphNodePos>::new(
                            curr.node, node_pos,
                        )),
                        None,
                    ));
                }
                AlignState::Insertion => {
                    alignment.push(AlignedPair::new(None, Some(curr_offset.as_usize() - 1)));
                }
                AlignState::Deletion2 | AlignState::Insertion2 => {
                    panic!("Unexpected state: {:?}", curr.state)
                }
            }

            curr = prev;
            curr_offset = prev_offset;
        }
        alignment.reverse();
        alignment
    }
}

impl<'a, G, T> AstarState<G> for &'a mut T
where
    G: AlignableGraph,
    T: AstarState<G>,
{
    type AstarItem = T::AstarItem;

    fn pop_front(&mut self) -> Option<Self::AstarItem> {
        (**self).pop_front()
    }

    fn is_further(&self, item: &Self::AstarItem, offset: usize) -> bool {
        (**self).is_further(item, offset)
    }

    fn is_end(&self, graph: &G, item: &Self::AstarItem) -> bool {
        (**self).is_end(graph, item)
    }

    fn get_score(&self, item: &Self::AstarItem) -> Score {
        (**self).get_score(item)
    }

    fn get_offset(&self, item: &Self::AstarItem) -> usize {
        (**self).get_offset(item)
    }

    fn is_visited(&self, item: &Self::AstarItem) -> bool {
        (**self).is_visited(item)
    }

    fn set_visited(&mut self, item: &Self::AstarItem) {
        (**self).set_visited(item)
    }

    fn update_if_further(&mut self, item: &Self::AstarItem, offset: usize) -> bool {
        (**self).update_if_further(item, offset)
    }

    fn queue_item(&mut self, item: Self::AstarItem, heuristic: usize) {
        (**self).queue_item(item, heuristic);
    }

    fn relax<F>(&mut self, graph: &G, seq: &[u8], item: &Self::AstarItem, heuristic: F)
    where
        F: Fn(&Self::AstarItem) -> usize,
    {
        (**self).relax(graph, seq, item, heuristic);
    }

    fn backtrace(&self, graph: &G, end: &Self::AstarItem) -> Vec<AlignedPair<G::NodePosType>> {
        (**self).backtrace(graph, end)
    }
}

/// A queue layer (bucket) for the gap-affine scoring model
///
/// Keep queued alignment graph nodes in different alignment states
/// in different vectors. When popping a state from the queue, we prioritize
/// matches, then deletions, and then insertions. This reduces branch prediction misses
/// in the main A* loop as compared to a single queue with all states mixed.
#[derive(Clone)]
pub struct AffineQueueLayer<N, D> {
    queued_states_m: Vec<(Score, N, Diag<D>)>,
    queued_states_i: Vec<(Score, N, Diag<D>)>,
    queued_states_d: Vec<(Score, N, Diag<D>)>,
}

impl<N, D> QueueLayer for AffineQueueLayer<N, D>
where
    N: Clone + fmt::Debug,
    D: DiagType,
{
    type QueueItem = AffineAstarItem<N, D>;

    fn queue(&mut self, item: Self::QueueItem) {
        // Since we have separate queues for each alignment state, we don't store the enum value to save memory
        match item.state {
            AlignState::Match => self
                .queued_states_m
                .push((item.score, item.node, item.diag)),
            AlignState::Deletion => self
                .queued_states_d
                .push((item.score, item.node, item.diag)),
            AlignState::Insertion => self
                .queued_states_i
                .push((item.score, item.node, item.diag)),
            AlignState::Insertion2 | AlignState::Deletion2 => {
                panic!("Invalid gap-affine state {:?}", item.state)
            }
        }
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        // Prioritize deletions, then insertions, then matches. This makes sure closing indels also update FR points, preventing
        // match states lower on the diagonals.
        // Construct the `AffineAstarItem` from on the fly to include the alignment state again.
        self.queued_states_d
            .pop()
            .map(|(score, node, diag)| {
                AffineAstarItem::new(score, node, diag, AlignState::Deletion)
            })
            .or_else(|| {
                self.queued_states_i
                    .pop()
                    .map(|(score, node, diag)| {
                        AffineAstarItem::new(score, node, diag, AlignState::Insertion)
                    })
                    .or_else(|| {
                        self.queued_states_m.pop().map(|(score, node, diag)| {
                            AffineAstarItem::new(score, node, diag, AlignState::Match)
                        })
                    })
            })
    }

    fn is_empty(&self) -> bool {
        self.queued_states_m.is_empty()
            && self.queued_states_d.is_empty()
            && self.queued_states_i.is_empty()
    }

    fn capacity(&self) -> usize {
        self.queued_states_m.capacity()
            + self.queued_states_d.capacity()
            + self.queued_states_i.capacity()
    }
}

impl<N, D> Default for AffineQueueLayer<N, D>
where
    D: DiagType,
{
    fn default() -> Self {
        Self {
            queued_states_m: Vec::with_capacity(16),
            queued_states_i: Vec::with_capacity(4),
            queued_states_d: Vec::with_capacity(4),
        }
    }
}
