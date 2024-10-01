use std::fmt;

use super::AlignmentCostModel;
use crate::aligner::{
    astar::{
        queue::{LayeredQueue, QueueLayer}, AlignState, AlignableGraphNodeId, AlignableGraphNodePos, AlignableGraphRef, AstarState
    }, extension::{extend, ExtendResult}, fr_points::{to_node_pos, Diag, DiagType, NodeFrPoints, OffsetType, Score}, utils::AlignedPair, AlignmentMode
};

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
            G: AlignableGraphRef,
            D: DiagType,
            O: OffsetType;
        
    fn initialize<G, D, O>(&self, graph: G, seq: &[u8], alignment_mode: AlignmentMode) -> Self::AstarStateType<G, D, O>
    where
        G: AlignableGraphRef,
        D: DiagType,
        O: OffsetType 
    {
        AffineAstarState::new(*self, graph, seq, alignment_mode)
    }

    fn initial_states<G, D, O>(
        &self,
        graph: G,
        _: &[u8],
        _: AlignmentMode,
    ) -> impl Iterator<Item = <Self::AstarStateType<G, D, O> as AstarState<G>>::AstarItem>
    where
        G: AlignableGraphRef,
        D: DiagType,
        O: OffsetType,
    {
        let start = AffineAstarItem::new(
            Score::default(),
            graph.start_node(),
            Diag::default(),
            AlignState::Match,
        );

        vec![start].into_iter()
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
    O: OffsetType,
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
    G: AlignableGraphRef,
{
    /// Alignment costs
    costs: Affine, 
    
    /// Graph we are aligning to
    graph: G,

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
    G: AlignableGraphRef,
    D: DiagType,
    O: OffsetType,
{
    fn new(costs: Affine, graph: G, seq: &[u8], alignment_mode: AlignmentMode) -> Self {
        Self {
            costs,
            graph,
            seq_length: seq.len(),
            alignment_mode,
            fr_points: vec![AffineNodeDiagonals::default(); graph.node_bound()],
            queue: LayeredQueue::new(),
        }
    }
    
    fn extend_and_relax_match<F>(&mut self, seq: &[u8], item: &AffineAstarItem<G::NodeType, D>, heuristic: F)
    where
        F: Fn(&AffineAstarItem<G::NodeType, D>) -> usize,
    {
        let fr_point = self.fr_points[item.node.index()]
            .get_furthest(item)
            .unwrap_or(O::zero());

        eprintln!("Offset before extend: {fr_point:?}, diagonal: {:?}", item.diag);
        let extend_result = extend(self.graph, seq, item.node, item.diag, fr_point);
        
        match extend_result {
            ExtendResult::NodeEnd(extended_offset) => {
                eprintln!("(NodeEnd) Extended offset: {extended_offset:?} {:?}", item.node);
                self.update_if_further(item, extended_offset.as_usize());
                
                let parent_node_len = self.graph.node_len(item.node);
                let mut any_mismatch = false;
                for succ in self.graph.successors(item.node) {
                    let succ_diag = item.diag + parent_node_len;
                    let succ_seq = self.graph.node_seq(succ);
                    let succ_len = succ_seq.len();
                    eprintln!("Succ: {:?} (diag: {succ_diag:?})", succ);
                    
                    if succ == self.graph.end_node() {
                        eprintln!("Reached end node");
                        
                        // We reached the end node, immediately queue with zero cost
                        let new_item = AffineAstarItem::new(
                            item.score, 
                            succ, 
                            succ_diag, 
                            AlignState::Match
                        );
                        
                        if self.update_if_further(&new_item, extended_offset.increase_one().as_usize()) {
                            let h = heuristic(&new_item);
                            self.queue_item(new_item, h);
                        }
                        
                        continue;
                    }
                    
                    let mut mismatch = false;
                    if extended_offset.as_usize() < self.seq_length {
                        // Check if we have a match or mismatch
                        let succ_qry_offset = extended_offset.increase_one();
                        let succ_pos = to_node_pos(succ_diag, succ_qry_offset.as_usize());
                        eprintln!("Succ pos: {succ_pos} (len: {succ_len}), succ offset: {succ_qry_offset:?} (seq len: {})", self.seq_length);
                        let succ_item = if succ_seq[succ_pos] == seq[succ_qry_offset.as_usize() - 1] {
                            // We have a successor with a match
                            eprintln!("match {:?} == {:?}", char::from(succ_seq[succ_pos]), char::from(seq[succ_qry_offset.as_usize()-1]));
                            AffineAstarItem::new(
                                item.score,
                                succ,
                                succ_diag,
                                AlignState::Match,
                            )
                        } else {
                            // We have a successor with a mismatch
                            eprintln!("mismatch {:?} != {:?}", char::from(succ_seq[succ_pos]), char::from(seq[succ_qry_offset.as_usize()-1]));
                            mismatch = true;
                            AffineAstarItem::new(
                                item.score + self.costs.mismatch(),
                                succ,
                                succ_diag,
                                AlignState::Match,
                            )
                        };
                        
                        if self.fr_points[succ.index()].update_if_further(&succ_item, succ_qry_offset) {
                            let h = heuristic(&succ_item);
                            self.queue_item(succ_item, h)
                        }
                    }
                    
                    // Check if we should open a deletion
                    if extended_offset.as_usize() == self.seq_length || mismatch {
                        let new_item_del = AffineAstarItem::new(
                            item.score + self.costs.gap_open() + self.costs.gap_extend(),
                            succ,
                            succ_diag - 1isize,
                            AlignState::Deletion,
                        );
                        
                        if self.fr_points[succ.index()].update_if_further(&new_item_del, extended_offset) {
                            let h = heuristic(&new_item_del);
                            self.queue_item(new_item_del, h)
                        }
                    }
                    
                    any_mismatch |= mismatch;
                }
                
                // Besides moving to successor nodes, see if we can open an insertion
                if any_mismatch && extended_offset.as_usize() < self.seq_length {
                    let new_item_ins = AffineAstarItem::new(
                        item.score + self.costs.gap_open() + self.costs.gap_extend(),
                        item.node,
                        item.diag + 1isize,
                        AlignState::Insertion,
                    );
                    
                    if self.fr_points[item.node.index()].update_if_further(&new_item_ins, extended_offset.increase_one()) {
                        let h = heuristic(&new_item_ins);
                        self.queue_item(new_item_ins, h)
                    }
                }
            }
            ExtendResult::OtherEnd(extended_offset) => {
                eprintln!("Extended offset: {extended_offset:?}");
                
                self.update_if_further(item, extended_offset.as_usize());
                
                // Extension ended within a node, queue the mismatch and open indels, expanding to other diagonals 
                if extended_offset.as_usize() < self.seq_length {
                    let new_item_mis = AffineAstarItem::new(
                        item.score + self.costs.mismatch(),
                        item.node,
                        item.diag,
                        AlignState::Match,
                    );
                    
                    if self.fr_points[item.node.index()].update_if_further(&new_item_mis, extended_offset.increase_one()) {
                        let h = heuristic(&new_item_mis);
                        self.queue_item(new_item_mis, h)
                    }
                    
                    let new_item_ins = AffineAstarItem::new(
                        item.score + self.costs.gap_open() + self.costs.gap_extend(),
                        item.node,
                        item.diag + 1isize,
                        AlignState::Insertion,
                    );
                    
                    if self.fr_points[item.node.index()].update_if_further(&new_item_ins, extended_offset.increase_one()) {
                        let h = heuristic(&new_item_ins);
                        self.queue_item(new_item_ins, h)
                    }
                }
                
                let new_item_del = AffineAstarItem::new(
                    item.score + self.costs.gap_open() + self.costs.gap_extend(),
                    item.node,
                    item.diag - 1isize,
                    AlignState::Deletion,
                );
                
                if self.fr_points[item.node.index()].update_if_further(&new_item_del, extended_offset) {
                    let h = heuristic(&new_item_del);
                    self.queue_item(new_item_del, h)
                }
            }
        }
    }
    
    fn relax_deletion<F>(&mut self, _: &[u8], item: &AffineAstarItem<G::NodeType, D>, heuristic: F)
    where
        F: Fn(&AffineAstarItem<G::NodeType, D>) -> usize,
    {
        let fr_point = self.fr_points[item.node.index()]
            .get_furthest(item)
            .unwrap_or(O::zero());
        
        let node_len = self.graph.node_len(item.node);
        let node_pos = to_node_pos(item.diag, fr_point.as_usize());
        
        if node_pos == node_len - 1 {
            // At the end of a node, we need to check successors to extend the deletion
            for succ in self.graph.successors(item.node) {
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
        let new_item_match = AffineAstarItem::new(
            item.score,
            item.node,
            item.diag,
            AlignState::Match,
        );
        
        if self.update_if_further(&new_item_match, fr_point.as_usize()) {
            let h = heuristic(&new_item_match);
            self.queue_item(new_item_match, h)
        }
    }
    
    fn relax_insertion<F>(&mut self, _: &[u8], item: &AffineAstarItem<G::NodeType, D>, heuristic: F)
    where
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
        let new_item_match = AffineAstarItem::new(
            item.score,
            item.node,
            item.diag,
            AlignState::Match,
        );
        
        if self.update_if_further(&new_item_match, fr_point.as_usize()) {
            let h = heuristic(&new_item_match);
            self.queue_item(new_item_match, h)
        }
    }
    
    fn get_prev(&self, item: &AffineAstarItem<G::NodeType, D>, curr_offset: O) -> Option<(AffineAstarItem<G::NodeType, D>, O)> 
    where
        <G as AlignableGraphRef>::Predecessors: IntoIterator<Item = G::NodeType>,
    {
        let mut sources = Vec::new();
        
        match item.state {
            AlignState::Match => {
                // We could have entered this diagonal from a match at a node predecessor, mismatch on the same node or a predecessor,
                // a deletion on the same node or a predecessor, or an insertion on the same node.
                
                // check for match at a predecessor, which should have been reached with the same score as the current item.
                for pred in self.graph.predecessors(item.node) {
                    let pred_len = self.graph.node_len(pred);
                    let pred_match_diag = item.diag - pred_len;
                    
                    let pred_item = AffineAstarItem::new(
                        item.score, 
                        pred, 
                        pred_match_diag, 
                        AlignState::Match
                    );
                    eprintln!("- Checkin for match {:?}", pred_item);
                    if let Some(pred_offset) = self.fr_points[pred.index()]
                        .get_furthest(&pred_item) 
                    {
                        eprintln!("- Curr offset: {curr_offset:?}, pred offset: {pred_offset:?}, isvisited: {}", self.is_visited(&pred_item));
                        if self.is_visited(&pred_item) && pred_offset < curr_offset {
                            sources.push((pred_item, pred_offset));
                        }
                    }
                }
                
                // Check for a mismatch, first on the same node, then on a predecessor
                // Predecessor should have a score of current score - mismatch cost
                if item.score >= Score::new(self.costs.mismatch() as u32) {
                    let mut possible_pred = vec![item.node];
                    possible_pred.extend(self.graph.predecessors(item.node));
                    
                    for pred in possible_pred {
                        let pred_len = if pred != item.node { self.graph.node_len(pred) } else { 0 };
                        let pred_mismatch_diag = item.diag - pred_len;
                        let pred_item = AffineAstarItem::new(
                            item.score - self.costs.mismatch(), 
                            pred, 
                            pred_mismatch_diag, 
                            AlignState::Match
                        );
                        
                        eprintln!("- Checkin for mismatch {:?}", pred_item);
                        if let Some(pred_offset) = self.fr_points[pred.index()]
                            .get_furthest(&pred_item)
                        {
                            eprintln!("- Curr offset: {curr_offset:?}, pred offset: {pred_offset:?}, isvisited: {}", self.is_visited(&pred_item));
                            if self.is_visited(&pred_item) && pred_offset < curr_offset {
                                sources.push((pred_item, pred_offset));
                            }
                        }
                    }
                }
                
                // Check if we closed an insertion, will be the same diagonal, score and offset
                let pred_ins_diag = item.diag;
                let pred_ins_item = AffineAstarItem::new(
                    item.score,
                    item.node,
                    pred_ins_diag,
                    AlignState::Insertion,
                );
                eprintln!("- Checkin for insertion {:?}", pred_ins_item);
                if let Some(pred_offset) = self.fr_points[item.node.index()]
                    .get_furthest(&pred_ins_item)
                {
                    eprintln!("- Curr offset: {curr_offset:?}, pred offset: {pred_offset:?}, isvisited: {}", self.is_visited(&pred_ins_item));
                    if self.is_visited(&pred_ins_item) && pred_offset < curr_offset {
                        sources.push((pred_ins_item, pred_offset));
                    }
                }
                
                // Check if we closed a deletion, will be the same diagonal, score and offset
                let pred_del_diag = item.diag;
                let pred_del_item = AffineAstarItem::new(
                    item.score,
                    item.node,
                    pred_del_diag,
                    AlignState::Deletion,
                );
                eprintln!("- Checkin for deletion {:?}", pred_del_item);
                if let Some(pred_offset) = self.fr_points[item.node.index()]
                    .get_furthest(&pred_del_item)
                {
                    eprintln!("- Curr offset: {curr_offset:?}, pred offset: {pred_offset:?}, isvisited: {}", self.is_visited(&pred_del_item));
                    if self.is_visited(&pred_del_item) && pred_offset <= curr_offset {
                        sources.push((pred_del_item, pred_offset));
                    }
                }
            },
            AlignState::Deletion => {
                let mut possible_pred = vec![item.node];
                possible_pred.extend(self.graph.predecessors(item.node));
                
                // Check first for opened deletions
                let gap_open_cost = self.costs.gap_open() as u32 + self.costs.gap_extend() as u32;
                if item.score >= Score::new(gap_open_cost) {
                    for pred in possible_pred.iter().copied() {
                        let pred_len = self.graph.node_len(pred);
                        let pred_diag = if item.node != pred { item.diag - pred_len + 1isize } else { item.diag + 1isize };
                        let pred_open = AffineAstarItem::new(
                            item.score - self.costs.gap_open() - self.costs.gap_extend(),
                            pred,
                            pred_diag,
                            AlignState::Match,
                        );
                        
                        if let Some(pred_offset) = self.fr_points[pred.index()]
                            .get_furthest(&pred_open)
                        {
                            if self.is_visited(&pred_open) && pred_offset == curr_offset {
                                sources.push((pred_open, pred_offset));
                            }
                        }
                    }
                }
                
                // Check for extended deletions
                if item.score >= Score::new(self.costs.gap_extend() as u32) {
                    for pred in possible_pred.iter().copied() {
                        let pred_len = self.graph.node_len(pred);
                        let pred_diag = if item.node != pred { item.diag - pred_len + 1isize } else { item.diag + 1isize };
                        let pred_ext = AffineAstarItem::new(
                            item.score - self.costs.gap_extend(),
                            pred,
                            pred_diag,
                            AlignState::Deletion,
                        );
                        
                        if let Some(pred_offset) = self.fr_points[pred.index()]
                            .get_furthest(&pred_ext)
                        {
                            if self.is_visited(&pred_ext) && pred_offset == curr_offset {
                                sources.push((pred_ext, pred_offset));
                            }
                        }
                    }
                }
            },
            AlignState::Insertion => {
                if curr_offset > O::zero() {
                    let gap_open_cost = self.costs.gap_open() as u32 + self.costs.gap_extend() as u32;
                    let pred_diag = item.diag - 1isize;
                    if item.score >= Score::new(gap_open_cost) {
                        let pred_open = AffineAstarItem::new(
                            item.score - self.costs.gap_open() - self.costs.gap_extend(),
                            item.node,
                            pred_diag,
                            AlignState::Match,
                        );
                        
                        eprintln!("- Checkin for ins open {:?}", pred_open);
                        if let Some(pred_offset) = self.fr_points[item.node.index()]
                            .get_furthest(&pred_open)
                        {
                            eprintln!("- Curr offset: {curr_offset:?}, pred offset: {pred_offset:?}, isvisited: {}", self.is_visited(&pred_open));
                            if self.is_visited(&pred_open) && pred_offset.increase_one() == curr_offset {
                                sources.push((pred_open, pred_offset));
                            }
                        }
                    }
                    
                    if item.score >= Score::new(self.costs.gap_extend() as u32) {
                        let pred_ext = AffineAstarItem::new(
                            item.score - self.costs.gap_extend(),
                            item.node,
                            pred_diag,
                            AlignState::Insertion,
                        );
                        
                        eprintln!("- Checkin for ins ext {:?}", pred_ext);
                        if let Some(pred_offset) = self.fr_points[item.node.index()]
                            .get_furthest(&pred_ext)
                        {
                            eprintln!("- Curr offset: {curr_offset:?}, pred offset: {pred_offset:?}, isvisited: {}", self.is_visited(&pred_ext));
                            if self.is_visited(&pred_ext) && pred_offset.increase_one() == curr_offset {
                                sources.push((pred_ext, pred_offset));
                            }
                        }
                    }
                }
            },
            AlignState::Deletion2 | AlignState::Insertion2 => panic!("Invalid gap-affine state {:?}", item.state)
        }
        
        sources.into_iter()
            .max_by_key(|(_, offset)| *offset)
    }
}

impl<G, D, O> AstarState<G> for AffineAstarState<G, D, O>
where
    G: AlignableGraphRef,
    D: DiagType,
    O: OffsetType,
{
    type AstarItem = AffineAstarItem<G::NodeType, D>;

    fn pop_front(&mut self) -> Option<Self::AstarItem> {
        self.queue.pop()
    }

    fn is_further(&self, item: &Self::AstarItem, offset: usize) -> bool {
        self.fr_points[item.node.index()].is_further(item, O::new(offset))
    }
    
    fn is_end(&self, item: &Self::AstarItem) -> bool {
        if item.node != self.graph.end_node() {
            return false;
        }
        
        let end_diag = Diag::new(self.seq_length as isize + 1);
        eprintln!("End diag: {end_diag:?}, item.diag: {:?}", item.diag);
        if item.diag != end_diag {
            return false;
        }
        
        eprintln!("Offset {:?}, required seq length: {:?}", self.get_offset(item), self.seq_length + 1);
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
        self.fr_points[item.node.index()].set_visited(item)
    }

    fn update_if_further(&mut self, item: &Self::AstarItem, offset: usize) -> bool {
        self.fr_points[item.node.index()].update_if_further(item, O::new(offset))
    }

    fn queue_item(&mut self, item: Self::AstarItem, heuristic: usize) {
        let priority = item.score.as_usize() + heuristic;
        
        let offset = self.get_offset(&item);
        eprintln!("Queuing: {item:?}, offset: {offset}, priority: {priority}");

        self.queue.queue(item, priority);
    }

    fn relax<F>(&mut self, seq: &[u8], item: &Self::AstarItem, heuristic: F)
    where
        F: Fn(&Self::AstarItem) -> usize,
    {
        match item.state {
            AlignState::Match => self.extend_and_relax_match(seq, item, heuristic),
            AlignState::Deletion => self.relax_deletion(seq, item, heuristic),
            AlignState::Insertion => self.relax_insertion(seq, item, heuristic),
            AlignState::Deletion2 | AlignState::Insertion2 => {
                panic!("Unexpected state: {:?}", item.state)
            }
        }
    }
    
    fn backtrace(&self, end: &Self::AstarItem) -> Vec<AlignedPair<G::NodePosType>> {
        let mut curr = end.clone();
        let mut curr_offset = O::new(self.seq_length + 1);
        let mut alignment = Vec::new();
        
        while let Some((prev, prev_offset)) = self.get_prev(&curr, curr_offset) {
            eprintln!("Current: {curr:?}, offset: {curr_offset:?}");
            eprintln!("Prev: {prev:?}, offset: {prev_offset:?}");
            eprintln!();
            
            if curr.node == self.graph.end_node() {
                curr = prev;
                curr_offset = prev_offset;
                continue;
            }
            
            match curr.state {
                AlignState::Match => {
                    for query_offset in (prev_offset.as_usize()+1..=curr_offset.as_usize()).rev() {
                        let node_pos = to_node_pos(curr.diag, query_offset);
                        
                        alignment.push(AlignedPair::new(
                            Some(<G::NodePosType as AlignableGraphNodePos>::new(curr.node, node_pos)),
                            Some(query_offset-1)
                        ));
                    }
                },
                AlignState::Deletion => {
                    let node_pos = to_node_pos(curr.diag, curr_offset.as_usize());
                    alignment.push(AlignedPair::new(
                        Some(<G::NodePosType as AlignableGraphNodePos>::new(curr.node, node_pos)),
                        None
                    ));
                },
                AlignState::Insertion => {
                    alignment.push(AlignedPair::new(
                        None,
                        Some(curr_offset.as_usize()-1)
                    ));
                },
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
    G: AlignableGraphRef,
    T: AstarState<G>,
{
    type AstarItem = T::AstarItem;

    fn pop_front(&mut self) -> Option<Self::AstarItem> {
        (**self).pop_front()
    }

    fn is_further(&self, item: &Self::AstarItem, offset: usize) -> bool {
        (**self).is_further(item, offset)
    }
    
    fn is_end(&self, item: &Self::AstarItem) -> bool {
        (**self).is_end(item)
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

    fn relax<F>(&mut self, seq: &[u8], item: &Self::AstarItem, heuristic: F)
    where
        F: Fn(&Self::AstarItem) -> usize,
    {
        (**self).relax(seq, item, heuristic);
    }
    
    fn backtrace(&self, end: &Self::AstarItem) -> Vec<AlignedPair<G::NodePosType>> {
        (**self).backtrace(end)
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
            .map(|(score, node, diag)| AffineAstarItem::new(score, node, diag, AlignState::Deletion))
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

// struct AffineDiagonals<O, const B: usize> {
//     fr_points_m: [O; B],
//     fr_points_i: [O; B],
//     fr_points_d: [O; B],
// }

// struct BlockedAffineDiagonals<O, const B: usize=16> {
//     /// For each node in the POA graph, we store the query offsets per diagonal. For cache-efficiency,
//     /// we use a blocked storage scheme where we store the diagonals in blocks of size B.
//     fr_point_blocks: Vec<FxHashMap<usize, AffineDiagonals<O, B>>>,
// }

// impl<O, const B: usize> BlockedAffineDiagonals<O, B> {
//     fn new() -> Self {
//         if B & (B-1) != 0 {
//             panic!("BlockedAffineDiagonals block size must be a power of 2");
//         }

//         Self {
//             fr_point_blocks: vec![FxHashMap::default()]
//         }
//     }
// }
