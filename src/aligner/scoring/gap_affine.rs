use std::cmp::Ordering;
use std::error::Error;
use std::io::Write;
use rustc_hash::FxHashMap;
use smallvec::{smallvec, SmallVec};
use crate::aligner::{AlignedPair, Alignment};

use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentCosts;
use crate::aligner::state::{AlignState, StateGraph, StateGraphNode, Score};

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
    type StateGraphType<N, O> = StateGraphAffine<N, O>
    where
        N: NodeIndexType,
        O: OffsetType;

    fn new_state_graph<G, O>(&self, graph: &G) -> Self::StateGraphType<G::NodeIndex, O>
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        StateGraphAffine::new(*self, graph)
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
            AlignState::Start | AlignState::Match | AlignState::Mismatch => self.cost_gap_open,
            _ => panic!("Invalid current state {:?} for gap affine scoring model!", current_state)
        };

        gap_open as usize + (length * self.cost_gap_extend as usize)
    }
}


#[derive(Debug, Clone)]
pub struct AffineStateNode<N, O> {
    node: N,
    offset: O,
    state: AlignState
}

impl<N, O> AffineStateNode<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    fn new(node: N, offset: O, state: AlignState) -> Self {
        Self {
            node,
            offset,
            state
        }
    }
}

impl<N, O> StateGraphNode<N, O> for AffineStateNode<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    fn new(node: N, offset: O, state: AlignState) -> Self {
        Self { node, offset, state }
    }

    fn node(&self) -> N {
        self.node
    }

    fn offset(&self) -> O {
        self.offset
    }

    fn state(&self) -> AlignState {
        self.state
    }
}

#[derive(Debug, Clone)]
pub struct AffineNodeData<N, O> {
    score: Score,
    prev: Option<AffineStateNode<N, O>>
}

impl<N, O> AffineNodeData<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    fn new(score: Score, prev: AffineStateNode<N, O>) -> Self {
        Self {
            score,
            prev: Some(prev)
        }
    }
}

impl<N, O> Default for AffineNodeData<N, O> {
    fn default() -> Self {
        Self {
            score: Score::Unvisited,
            prev: None
        }
    }
}

#[derive(Clone)]
pub struct StateGraphAffine<N, O> {
    costs: GapAffine,
    nodes_m: Vec<FxHashMap<O, AffineNodeData<N, O>>>,
    nodes_i: Vec<FxHashMap<O, AffineNodeData<N, O>>>,
    nodes_d: Vec<FxHashMap<O, AffineNodeData<N, O>>>,
}

impl<N, O> StateGraphAffine<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    fn new<G: AlignableGraph>(costs: GapAffine, seq_graph: &G) -> Self {
        Self {
            costs,
            nodes_m: vec![FxHashMap::default(); seq_graph.node_count_with_start()],
            nodes_i: vec![FxHashMap::default(); seq_graph.node_count_with_start()],
            nodes_d: vec![FxHashMap::default(); seq_graph.node_count_with_start()],
        }
    }
}

impl<N, O> StateGraph<N, O> for StateGraphAffine<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    type StateNode = AffineStateNode<N, O>;
    type NewStatesContainer = SmallVec<[(Self::StateNode, u8); 1]>;

    fn get_score(&self, state: &Self::StateNode) -> Score {
        match state.state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                self.nodes_m[state.node().index()]
                    .get(&state.offset())
                    .map(|data| data.score)
                    .unwrap_or(Score::Unvisited)
            },
            AlignState::Insertion => {
                self.nodes_i[state.node().index()]
                    .get(&state.offset())
                    .map(|data| data.score)
                    .unwrap_or(Score::Unvisited)
            }
            AlignState::Deletion => {
                self.nodes_d[state.node().index()]
                    .get(&state.offset())
                    .map(|data| data.score)
                    .unwrap_or(Score::Unvisited)
            },
            _ => panic!("Invalid state {:?} for GapAffine scoring!", state.state())
        }
    }

    fn update_score(&mut self, state: &Self::StateNode, score: Score, prev: &Self::StateNode) {
        match state.state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                self.nodes_m[state.node().index()].entry(state.offset())
                    .and_modify(|entry| {
                        entry.score = score;
                        entry.prev = Some(prev.clone());
                    })
                    .or_insert_with(|| AffineNodeData::new(score, prev.clone()));
            },
            AlignState::Insertion => {
                self.nodes_i[state.node().index()].entry(state.offset())
                    .and_modify(|entry| {
                        entry.score = score;
                        entry.prev = Some(prev.clone());
                    })
                    .or_insert_with(|| AffineNodeData::new(score, prev.clone()));
            },
            AlignState::Deletion => {
                self.nodes_d[state.node().index()].entry(state.offset())
                    .and_modify(|entry| {
                        entry.score = score;
                        entry.prev = Some(prev.clone());
                    })
                    .or_insert_with(|| AffineNodeData::new(score, prev.clone()));
            },
            _ => panic!("Invalid state {:?} for GapAffine scoring!", state.state())
        }
    }

    fn get_prev(&self, state: &Self::StateNode) -> Option<&Self::StateNode> {
        match state.state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                self.nodes_m[state.node().index()]
                    .get(&state.offset())
                    .and_then(|data| data.prev.as_ref())
            },
            AlignState::Insertion => {
                self.nodes_i[state.node().index()]
                    .get(&state.offset())
                    .and_then(|data| data.prev.as_ref())
            }
            AlignState::Deletion => {
                self.nodes_d[state.node().index()]
                    .get(&state.offset())
                    .and_then(|data| data.prev.as_ref())
            },
            _ => panic!("Invalid state {:?} for GapAffine scoring!", state.state())
        }
    }

    fn new_match_state(&mut self, parent: &Self::StateNode, child: N, current_score: Score) -> Option<Self::StateNode> {
        let child_offset = parent.offset().increase_one();
        let child_state = Self::StateNode::new(child, child_offset, AlignState::Match);

        match current_score.cmp(&self.get_score(&child_state)) {
            Ordering::Less => {
                self.update_score(&child_state, current_score, parent);
                return Some(child_state)
            },
            Ordering::Equal => {
                // Check if we need to update back trace
                if let Some(prev_state) = self.get_prev(&child_state) {
                    if parent.state() < prev_state.state() {
                        self.update_score(&child_state, current_score, parent);
                    }
                }
            },
            Ordering::Greater => (),
        }

        None
    }

    fn new_mismatch_state(&mut self, parent: &Self::StateNode, child: N, current_score: Score) -> Option<(Self::StateNode, u8)> {
        let child_offset = parent.offset().increase_one();
        let child_state = Self::StateNode::new(child, child_offset, AlignState::Mismatch);

        let mismatch_score = current_score + self.costs.mismatch();
        match mismatch_score.cmp(&self.get_score(&child_state)) {
            Ordering::Less => {
                self.update_score(&child_state, mismatch_score, parent);
                return Some((child_state, self.costs.mismatch()))
            },
            Ordering::Equal => {
                // Check if we need to update back trace
                if let Some(prev_state) = self.get_prev(&child_state) {
                    if parent.state() < prev_state.state() {
                        self.update_score(&child_state, mismatch_score, parent);
                    }
                }
            },
            Ordering::Greater => (),
        }

        None
    }

    fn open_or_extend_insertion(&mut self, parent: &Self::StateNode, score: Score) -> Self::NewStatesContainer {
        let mut valid_ins_states: SmallVec<[(Self::StateNode, u8); 1]> = smallvec![];

        match parent.state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                let gap_open = self.costs.gap_open() + self.costs.gap_extend();
                let gap_open_score = score + gap_open;
                let ins_state = StateGraphNode::new(
                    parent.node(),
                    parent.offset().increase_one(),
                    AlignState::Insertion
                );

                match gap_open_score.cmp(&self.get_score(&ins_state)) {
                    Ordering::Less => {
                        self.update_score(&ins_state, gap_open_score, parent);
                        valid_ins_states.push((ins_state, gap_open))
                    },
                    Ordering::Equal => {
                        // Check if we need to update back trace
                        if let Some(prev_state) = self.get_prev(&ins_state) {
                            if parent.state() < prev_state.state() {
                                self.update_score(&ins_state, gap_open_score, parent);
                            }
                        }
                    },
                    Ordering::Greater => (),
                }
            },
            AlignState::Insertion => {
                let gap_extend_score = score + self.costs.gap_extend();
                let ins_state = StateGraphNode::new(
                    parent.node(),
                    parent.offset().increase_one(),
                    AlignState::Insertion
                );

                match gap_extend_score.cmp(&self.get_score(&ins_state)) {
                    Ordering::Less => {
                        self.update_score(&ins_state, gap_extend_score, parent);
                        valid_ins_states.push((ins_state, self.costs.gap_extend()))
                    },
                    Ordering::Equal => {
                        // Check if we need to update back trace
                        if let Some(prev_state) = self.get_prev(&ins_state) {
                            if parent.state() < prev_state.state() {
                                self.update_score(&ins_state, gap_extend_score, parent);
                            }
                        }
                    },
                    Ordering::Greater => (),
                }
            },
            _ => ()
        }

        valid_ins_states
    }

    fn extend_insertion(&mut self, parent: &Self::StateNode, score: Score) -> Option<(Self::StateNode, u8)> {
        if parent.state() != AlignState::Insertion {
            return None
        }

        let gap_extend_score = score + self.costs.gap_extend();
        let ins_state = StateGraphNode::new(
            parent.node(),
            parent.offset().increase_one(),
            AlignState::Insertion
        );

        match gap_extend_score.cmp(&self.get_score(&ins_state)) {
            Ordering::Less => {
                self.update_score(&ins_state, gap_extend_score, parent);
                Some((ins_state, self.costs.gap_extend()))
            },
            Ordering::Equal => {
                // Check if we need to update back trace
                if let Some(prev_state) = self.get_prev(&ins_state) {
                    if parent.state() < prev_state.state() {
                        self.update_score(&ins_state, gap_extend_score, parent);
                    }
                }

                None
            },
            Ordering::Greater => None,
        }
    }

    fn open_or_extend_deletion(
        &mut self,
        parent: &Self::StateNode,
        child_node: N,
        score: Score
    ) -> Self::NewStatesContainer {
        let mut valid_del_states: SmallVec<[(Self::StateNode, u8); 1]> = smallvec![];

        match parent.state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                let gap_open = self.costs.gap_open() + self.costs.gap_extend();
                let gap_open_score = score + gap_open;
                let del_state = StateGraphNode::new(
                    child_node,
                    parent.offset(),
                    AlignState::Deletion
                );

                match gap_open_score.cmp(&self.get_score(&del_state)) {
                    Ordering::Less => {
                        self.update_score(&del_state, gap_open_score, parent);
                        valid_del_states.push((del_state, gap_open))
                    },
                    Ordering::Equal => {
                        // Check if we need to update back trace
                        if let Some(prev_state) = self.get_prev(&del_state) {
                            if parent.state() < prev_state.state() {
                                self.update_score(&del_state, gap_open_score, parent);
                            }
                        }
                    },
                    Ordering::Greater => (),
                }
            },
            AlignState::Deletion => {
                let gap_extend_score = score + self.costs.gap_extend();
                let del_state = StateGraphNode::new(
                    child_node,
                    parent.offset(),
                    AlignState::Deletion
                );

                match gap_extend_score.cmp(&self.get_score(&del_state)) {
                    Ordering::Less => {
                        self.update_score(&del_state, gap_extend_score, parent);
                        valid_del_states.push((del_state, self.costs.gap_extend()))
                    },
                    Ordering::Equal => {
                        // Check if we need to update back trace
                        if let Some(prev_state) = self.get_prev(&del_state) {
                            if parent.state() < prev_state.state() {
                                self.update_score(&del_state, gap_extend_score, parent);
                            }
                        }
                    },
                    Ordering::Greater => (),
                }
            },
            _ => ()
        }

        valid_del_states
    }

    fn extend_deletion(&mut self, parent: &Self::StateNode, child: N, score: Score) -> Option<(Self::StateNode, u8)> {
        if parent.state() != AlignState::Deletion {
            return None
        }

        let gap_extend_score = score + self.costs.gap_extend();
        let del_state = StateGraphNode::new(
            child,
            parent.offset(),
            AlignState::Deletion
        );

        match gap_extend_score.cmp(&self.get_score(&del_state)) {
            Ordering::Less => {
                self.update_score(&del_state, gap_extend_score, parent);
                Some((del_state, self.costs.gap_extend()))
            },
            Ordering::Equal => {
                // Check if we need to update back trace
                if let Some(prev_state) = self.get_prev(&del_state) {
                    if parent.state() < prev_state.state() {
                        self.update_score(&del_state, gap_extend_score, parent);
                    }
                }

                None
            },
            Ordering::Greater => None,
        }
    }

    fn backtrace(&self, end_state: &Self::StateNode) -> Alignment<N> {
        let mut curr = Some(end_state);
        let mut alignment = Alignment::new();

        while let Some(state) = curr {
            let Some(bt) = self.get_prev(&state) else {
                break;
            };

            match state.state() {
                AlignState::Match | AlignState::Mismatch => {
                    alignment.push(AlignedPair { rpos: Some(state.node()), qpos: Some(state.offset().as_usize() - 1) });
                },
                AlignState::Insertion => {
                    alignment.push(AlignedPair { rpos: None, qpos: Some(state.offset().as_usize() - 1) });
                },
                AlignState::Deletion => {
                    alignment.push(AlignedPair { rpos: Some(state.node()), qpos: None });
                },
                AlignState::Start | AlignState::Insertion2 | AlignState::Deletion2 =>
                    panic!("Unexpected align state in backtrace!")
            }

            curr = match bt.state() {
                AlignState::Start => None,
                _ => Some(bt)
            }
        }

        alignment.reverse();
        alignment
    }

    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        writeln!(writer, "node_id\toffset\tmatrix\tscore")?;

        for i in 0..self.nodes_m.len() {
            for (offset, data) in self.nodes_m[i].iter() {
                writeln!(writer, "{}\t{:?}\tmatch\t{}", i, offset, data.score)?;
            }
            for (offset, data) in self.nodes_i[i].iter() {
                writeln!(writer, "{}\t{:?}\tinsertion\t{}", i, offset, data.score)?;
            }
            for (offset, data) in self.nodes_d[i].iter() {
                writeln!(writer, "{}\t{:?}\tdeletion\t{}", i, offset, data.score)?;
            }
        }

        Ok(())
    }
}
