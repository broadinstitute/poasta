use std::cmp::Ordering;
use rustc_hash::FxHashMap;
use crate::aligner::{AlignedPair, Alignment};

use crate::graphs::{AlignableRefGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, Score};
use crate::aligner::aln_graph::{AlignmentGraph, AlignmentGraphNode, AlignState};
use crate::aligner::astar::{AstarQueue, AstarQueuedItem, AstarVisited};
use crate::aligner::queue::{LayeredQueue, QueueLayer};
use crate::bubbles::index::BubbleIndex;
use crate::bubbles::reached::ReachedBubbleExits;

#[derive(Clone, Copy, Debug)]
pub struct GapAffine {
    cost_mismatch: u8,
    cost_gap_extend: u8,
    cost_gap_open: u8
}

impl GapAffine {
    pub fn new(cost_mismatch: u8, cost_gap_extend: u8, cost_gap_open: u8) -> Self {
        Self { cost_mismatch, cost_gap_extend, cost_gap_open }
    }
}

impl AlignmentCosts for GapAffine {
    type AlignmentGraphType = AffineAlignmentGraph;
    type QueueType<N, O> = AffineLayeredQueue<N, O>
        where N: NodeIndexType,
              O: OffsetType;

    fn new_alignment_graph(&self, aln_type: AlignmentType) -> Self::AlignmentGraphType {
        AffineAlignmentGraph::new(self, aln_type)
    }

    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.cost_mismatch
    }

    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.cost_gap_open
    }

    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.cost_gap_extend
    }

    #[inline(always)]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        0
    }

    #[inline]
    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize {
        let gap_open = match current_state {
            AlignState::Insertion | AlignState::Deletion => 0,
            AlignState::Match => self.cost_gap_open,
            _ => panic!("Invalid current state {:?} for gap affine scoring model!", current_state)
        };

        gap_open as usize + (length * self.cost_gap_extend as usize)
    }
}

pub struct AffineAlignmentGraph {
    /// Gap affine costs
    costs: GapAffine,

    /// Alignment type, e.g., global alignment or ends-free alignment
    aln_type: AlignmentType,
}

impl AffineAlignmentGraph {
    fn new(costs: &GapAffine, aln_type: AlignmentType) -> Self {
        Self {
            costs: *costs,
            aln_type,
        }
    }
}

impl AlignmentGraph for AffineAlignmentGraph {
    type CostModel = GapAffine;

    fn get_costs(&self) -> &Self::CostModel {
        &self.costs
    }

    fn initial_states<G, O>(&self, ref_graph: &G) -> Vec<AlignmentGraphNode<G::NodeIndex, O>>
        where G: AlignableRefGraph,
              O: OffsetType,
    {
        match self.aln_type {
            AlignmentType::Global => vec![AlignmentGraphNode::new(ref_graph.start_node(), O::zero())],
            AlignmentType::EndsFree {
                qry_free_begin: free_qry, qry_free_end: _,
                graph_free_begin: free_grap, graph_free_end: _ }
            => {
                todo!();
            }
        }
    }

    fn is_end<G, O>(&self, ref_graph: &G, seq: &[u8], node: &AlignmentGraphNode<G::NodeIndex, O>, aln_state: AlignState) -> bool
        where G: AlignableRefGraph,
              O: OffsetType,
    {
        aln_state == AlignState::Match
            && node.node() == ref_graph.end_node()
            && node.offset().as_usize() == seq.len()
    }

    fn expand_all<V, G, O, F>(
        &self, visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        state: AlignState,
        mut f: F
    )
        where V: AstarVisited<G::NodeIndex, O>,
              G: AlignableRefGraph,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState)
    {
        match state {
            AlignState::Match => {
                let child_offset = node.offset().increase_one();
                for ref_succ in ref_graph.successors(node.node()) {
                    let new_node_mis = AlignmentGraphNode::new(ref_succ, child_offset);

                    // Move to next (mis)match state
                    let score_delta = if ref_graph.is_symbol_equal(ref_succ, seq[child_offset.as_usize()-1]) {
                        0u8
                    } else {
                        self.costs.cost_mismatch
                    };
                    let new_score_mis = score + score_delta;

                    if visited_data
                        .update_score_if_lower(new_score_mis, &new_node_mis, AlignState::Match, node, state)
                    {
                        f(score_delta, new_node_mis, AlignState::Match)
                    }

                    // Open deletion
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
                    if visited_data
                        .update_score_if_lower(new_score_del, &new_node_del, AlignState::Deletion, node, state)
                    {
                        f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Deletion);

                        // Also immediately traverse zero-cost D->M edge
                        if visited_data
                            .update_score_if_lower(new_score_del, &new_node_del, AlignState::Match, &new_node_del, AlignState::Deletion)
                        {
                            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Match);
                        }
                    }
                }

                // Open insertion
                let new_node_ins = AlignmentGraphNode::new(node.node(), child_offset);
                let new_score_ins = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
                if child_offset.as_usize() <= seq.len() && visited_data
                    .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Insertion, node, state)
                {
                    f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);

                    // Also immediately traverse zero-cost I->M edge
                    if visited_data
                        .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Match, &new_node_ins, AlignState::Insertion)
                    {
                        f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Match);
                    }
                }
            },
            AlignState::Insertion => {
                let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
                let new_score_ins = score + self.costs.cost_gap_extend;

                if node.offset().as_usize() < seq.len() && visited_data
                    .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Insertion, node, state)
                {
                    f(self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);

                    // Also immediately traverse zero-cost I->M edge
                    if visited_data
                        .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Match, &new_node_ins, AlignState::Insertion)
                    {
                        f(self.costs.cost_gap_extend, new_node_ins, AlignState::Match);
                    }
                }
            },
            AlignState::Deletion => {
                for ref_succ in ref_graph.successors(node.node()) {
                    // Extend deletion
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + self.costs.cost_gap_extend;
                    if visited_data.update_score_if_lower(new_score_del, &new_node_del, AlignState::Deletion, node, state) {
                        f(self.costs.cost_gap_extend, new_node_del, AlignState::Deletion);

                        // Also immediately traverse zero-cost D->M edge
                        if visited_data
                            .update_score_if_lower(new_score_del, &new_node_del, AlignState::Match, &new_node_del, AlignState::Deletion)
                        {
                            f(self.costs.cost_gap_extend, new_node_del, AlignState::Match);
                        }
                    }
                }
            },
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {state:?}!")
        }
    }

    fn expand_ref_graph_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
        if visited_data
            .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Insertion, parent, AlignState::Match)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);

            // Also immediately traverse zero-cost I->M edge
            if visited_data
                .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Match, &new_node_ins, AlignState::Insertion)
            {
                f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Match);
            }
        }
    }

    fn expand_query_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: N,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        let new_node_del = AlignmentGraphNode::new(child, parent.offset());
        let new_score_del = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
        if visited_data
            .update_score_if_lower(new_score_del, &new_node_del, AlignState::Deletion, parent, AlignState::Match)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Deletion);

            // Also immediately traverse zero-cost D->M edge
            if visited_data
                .update_score_if_lower(new_score_del, &new_node_del, AlignState::Match, &new_node_del, AlignState::Deletion)
            {
                f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Match);
            }
        }
    }

    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        f(self.costs.cost_mismatch, *child, AlignState::Match);

        // Also queue indel states from parent
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;

        if visited_data
            .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Insertion, parent, AlignState::Match)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);

            // Also immediately traverse zero-cost I->M edge
            if visited_data
                .update_score_if_lower(new_score_ins, &new_node_ins, AlignState::Match, &new_node_ins, AlignState::Insertion)
            {
                f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Match);
            }
        }

        let new_node_del = AlignmentGraphNode::new(child.node(), parent.offset());
        let new_score_del = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
        if visited_data
            .update_score_if_lower(new_score_del, &new_node_del, AlignState::Deletion, parent, AlignState::Match)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Deletion);

            // Also immediately traverse zero-cost D->M edge
            if visited_data
                .update_score_if_lower(new_score_del, &new_node_del, AlignState::Match, &new_node_del, AlignState::Deletion)
            {
                f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Match);
            }
        }
    }

}


type AffineAstarVisitedNodes<N, O> = Vec<FxHashMap<O, AffineVisitedNodeData<N, O>>>;

#[derive(Debug, Clone)]
pub struct AffineVisitedNodeData<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    score: Score,
    prev: Option<(AlignmentGraphNode<N, O>, AlignState)>
}

impl<N, O> AffineVisitedNodeData<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    fn new(score: Score, parent: &AlignmentGraphNode<N, O>, parent_state: AlignState) -> Self {
        Self {
            score,
            prev: Some((*parent, parent_state))
        }
    }
}

impl<N, O> Default for AffineVisitedNodeData<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    fn default() -> Self {
        Self {
            score: Score::Unvisited,
            prev: None
        }
    }
}



pub struct AffineAstarData<N, O>
where N: NodeIndexType,
      O: OffsetType,
{
    costs: GapAffine,
    bubble_index: BubbleIndex<N>,

    nodes_m: AffineAstarVisitedNodes<N, O>,
    nodes_i: AffineAstarVisitedNodes<N, O>,
    nodes_d: AffineAstarVisitedNodes<N, O>,

    bubbles_reached_m: ReachedBubbleExits<GapAffine, O>,
    bubbles_reached_i: ReachedBubbleExits<GapAffine, O>,
    bubbles_reached_d: ReachedBubbleExits<GapAffine, O>,
}

impl<N, O> AffineAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    pub fn new<G>(costs: GapAffine, ref_graph: &G, bubble_index: BubbleIndex<G::NodeIndex>) -> Self
        where G: AlignableRefGraph<NodeIndex=N>,
    {
        Self {
            costs,
            bubble_index,

            nodes_m: vec![FxHashMap::default(); ref_graph.node_count_with_start_and_end()],
            nodes_i: vec![FxHashMap::default(); ref_graph.node_count_with_start_and_end()],
            nodes_d: vec![FxHashMap::default(); ref_graph.node_count_with_start_and_end()],

            bubbles_reached_m: ReachedBubbleExits::new(costs, AlignState::Match, ref_graph),
            bubbles_reached_i: ReachedBubbleExits::new(costs, AlignState::Insertion, ref_graph),
            bubbles_reached_d: ReachedBubbleExits::new(costs, AlignState::Deletion, ref_graph),
        }
    }

    fn get_backtrace(
        &self,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState
    ) -> Option<(AlignmentGraphNode<N, O>, AlignState)> {
        match aln_state {
            AlignState::Match => self.nodes_m[aln_node.node().index()].get(&aln_node.offset())
                .and_then(|v| v.prev),
            AlignState::Deletion => self.nodes_d[aln_node.node().index()].get(&aln_node.offset())
                .and_then(|v| v.prev),
            AlignState::Insertion => self.nodes_i[aln_node.node().index()].get(&aln_node.offset())
                .and_then(|v| v.prev),
            AlignState::Deletion2 | AlignState::Insertion2 => panic!("Invalid gap-affine state {aln_state:?}")
        }
    }
}

impl<N, O> AstarVisited<N, O> for AffineAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    fn get_score(&self, aln_node: AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        let node_map = match aln_state {
            AlignState::Match => &self.nodes_m,
            AlignState::Insertion => &self.nodes_i,
            AlignState::Deletion => &self.nodes_d,
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state: {aln_state:?}")
        };

        node_map[aln_node.node().index()].get(&aln_node.offset())
            .map(|v| v.score)
            .unwrap_or(Score::Unvisited)
    }
    fn visit(&mut self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) {
        if self.bubble_index.is_exit(aln_node.node()) {
            let bubble_map = match aln_state {
                AlignState::Match => &mut self.bubbles_reached_m,
                AlignState::Insertion => &mut self.bubbles_reached_i,
                AlignState::Deletion => &mut self.bubbles_reached_d,
                AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine alignment state {aln_state:?}")
            };

            bubble_map.mark_reached(aln_node.node(), aln_node.offset(), score);
        }
    }

    fn prune(&self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool {
        if aln_state == AlignState::Insertion
            && !self.bubbles_reached_i.can_improve_alignment(&self.bubble_index, aln_node, score) {
            return true
        }

        if aln_state == AlignState::Deletion
            && !self.bubbles_reached_d.can_improve_alignment(&self.bubble_index, aln_node, score) {
            return true
        }


        !self.bubbles_reached_m.can_improve_alignment(&self.bubble_index, aln_node, score)
    }

    fn update_score_if_lower(
        &mut self, score: Score,
        aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState,
        parent: &AlignmentGraphNode<N, O>, parent_state: AlignState
    ) -> bool {
        let mut is_lower = false;

        let node_map = match aln_state {
            AlignState::Match => &mut self.nodes_m,
            AlignState::Insertion => &mut self.nodes_i,
            AlignState::Deletion => &mut self.nodes_d,
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid align state for gap-affine: {aln_state:?}")
        };

        node_map[aln_node.node().index()].entry(aln_node.offset())
            .and_modify(|data| {
                match score.cmp(&data.score) {
                    Ordering::Less => {
                        data.score = score;
                        data.prev = Some((*parent, parent_state));
                        is_lower = true;
                    },
                    Ordering::Equal => {
                        // Check if we need to update back trace
                        if let Some(ref mut prev_state) = &mut data.prev {
                            if parent_state < prev_state.1 {
                                *prev_state = (*parent, parent_state);
                            }
                        }
                    },
                    Ordering::Greater => (),
                }
            })
            .or_insert_with(|| {
                is_lower = true;
                AffineVisitedNodeData::new(score, parent, parent_state)
            });

        is_lower
    }

    fn backtrace<G>(&self, ref_graph: &G, aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
        where G: AlignableRefGraph<NodeIndex=N>,
    {
        let Some((mut curr, mut curr_state)) = self.get_backtrace(aln_node, AlignState::Match) else {
            panic!("No backtrace for alignment end state?");
        };

        let mut alignment = Alignment::new();

        while let Some((bt_node, bt_state)) = self.get_backtrace(&curr, curr_state) {
            if curr_state == AlignState::Match && (bt_state == AlignState::Insertion || bt_state == AlignState::Deletion) {
                curr = bt_node;
                curr_state = bt_state;
                continue;
            }

            match curr_state {
                AlignState::Match => {
                    alignment.push(AlignedPair { rpos: Some(curr.node()), qpos: Some(curr.offset().as_usize() - 1) });
                },
                AlignState::Insertion => {
                    alignment.push(AlignedPair { rpos: None, qpos: Some(curr.offset().as_usize() - 1) });
                },
                AlignState::Deletion => {
                    alignment.push(AlignedPair { rpos: Some(curr.node()), qpos: None });
                },
                AlignState::Insertion2 | AlignState::Deletion2 =>
                    panic!("Unexpected align state {curr_state:?} in backtrace!")
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
}


/// A queue layer (bucket) for the gap-affine scoring model
///
/// Keep queued alignment graph nodes in different alignment states
/// in different vectors to reduce branch prediction misses in the main A* loop.
#[derive(Default, Clone)]
pub struct AffineQueueLayer<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    queued_states_m: Vec<(Score, AlignmentGraphNode<N, O>)>,
    queued_states_i: Vec<(Score, AlignmentGraphNode<N, O>)>,
    queued_states_d: Vec<(Score, AlignmentGraphNode<N, O>)>,
}

impl<N, O> QueueLayer for AffineQueueLayer<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    type QueueItem = AstarQueuedItem<N, O>;

    fn queue(&mut self, item: Self::QueueItem) {
        match item.aln_state() {
            AlignState::Match => self.queued_states_m.push((item.score(), item.aln_node())),
            AlignState::Insertion => self.queued_states_i.push((item.score(), item.aln_node())),
            AlignState::Deletion => self.queued_states_d.push((item.score(), item.aln_node())),
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {:?}", item.aln_state())
        }
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        self.queued_states_m
            .pop()
            .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Match))
            .or_else(|| self.queued_states_d
                .pop()
                .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Deletion))
                .or_else(|| self.queued_states_i
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Insertion))
                )
            )
    }

    fn is_empty(&self) -> bool {
        self.queued_states_m.is_empty()
            && self.queued_states_d.is_empty()
            && self.queued_states_i.is_empty()
    }
}


type AffineLayeredQueue<N, O> = LayeredQueue<AffineQueueLayer<N, O>>;

impl<N, O> AstarQueue<N, O> for AffineLayeredQueue<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>> {
        self.pop()
    }

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score, h: usize) {
        let priority: usize = usize::from(score) + h;
        let item = AstarQueuedItem(score, node, aln_state);

        self.queue(item, priority)
    }
}


//
//     fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
//         writeln!(writer, "node_id\toffset\tmatrix\tscore")?;
//
//         for i in 0..self.nodes_m.len() {
//             for (offset, data) in self.nodes_m[i].iter() {
//                 writeln!(writer, "{}\t{:?}\tmatch\t{}", i, offset, data.score)?;
//             }
//             for (offset, data) in self.nodes_i[i].iter() {
//                 writeln!(writer, "{}\t{:?}\tinsertion\t{}", i, offset, data.score)?;
//             }
//             for (offset, data) in self.nodes_d[i].iter() {
//                 writeln!(writer, "{}\t{:?}\tdeletion\t{}", i, offset, data.score)?;
//             }
//         }
//
//         Ok(())
//     }
// }
