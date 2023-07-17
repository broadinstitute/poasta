use std::fmt::{Display, Formatter};
use std::collections::HashSet;

use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::io::state_tree::write_tree_gml;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentCostsAffine, AlignmentCostsEdit, AlignmentCostsLinear, AlignmentStateTree};
use crate::aligner::state::{AlignState, StateTree, StateTreeNode, Backtrace, TreeIndexType};

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
    type StateTreeType<N, O, Ix> = GapAffineStateTree<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType;

    fn to_new_state_tree<N, O, Ix>(&self) -> Self::StateTreeType<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType
    {
        GapAffineStateTree::new(*self)
    }
}

impl AlignmentCostsEdit for GapAffine {
    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.cost_mismatch
    }
}

impl AlignmentCostsLinear for GapAffine {
    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.cost_gap_extend
    }
}

impl AlignmentCostsAffine for GapAffine {
    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.cost_gap_open
    }
}


pub struct GapAffineStateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    costs: GapAffine,
    tree: StateTree<N, O, Ix>,
    visited_m: HashSet<(N, O)>,
    visited_d: HashSet<(N, O)>,
    visited_i: HashSet<(N, O)>,
}

impl<N, O, Ix> GapAffineStateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    fn new(costs: GapAffine) -> Self {
        Self {
            costs,
            tree: StateTree::new(),
            visited_m: HashSet::default(),
            visited_d: HashSet::default(),
            visited_i: HashSet::default()
        }
    }
}

impl<N, O, Ix> AlignmentStateTree<N, O, Ix> for GapAffineStateTree<N, O, Ix>
where
    N: NodeIndexType,
    O: OffsetType,
    Ix: TreeIndexType,
{
    fn add_node(&mut self, node: StateTreeNode<N, O, Ix>) -> Ix {
        self.tree.add_node(node)
    }

    #[inline(always)]
    fn get_node(&self, node_ix: Ix) -> &StateTreeNode<N, O, Ix> {
        self.tree.get_node(node_ix)
    }

    fn num_nodes(&self) -> usize {
        self.tree.num_nodes()
    }

    fn visited(&self, node: N, offset: O, state: AlignState) -> bool {
        let key = (node, offset);
        match state {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => self.visited_m.contains(&key),
            AlignState::Deletion => self.visited_d.contains(&key),
            AlignState::Insertion => self.visited_i.contains(&key),
            _ => panic!("Invalid alignment state for GapAffine!")
        }
    }

    fn mark_visited(&mut self, node: N, offset: O, state: AlignState) {
        let key = (node, offset);

        match state {
            AlignState::Start | AlignState::Match | AlignState::Mismatch => self.visited_m.insert(key),
            AlignState::Deletion => self.visited_d.insert(key),
            AlignState::Insertion => self.visited_i.insert(key),
            _ => panic!("Invalid alignment state for GapAffine!")
        };
    }

    fn close_indels_for(&mut self, node_indices: &[Ix]) -> Vec<Ix> {
        node_indices.iter().filter_map(|v| {
            let node = self.tree.get_node(*v);

            match node.state() {
                AlignState::Deletion | AlignState::Insertion => {
                    let key = (node.node(), node.offset());
                    if !self.visited_m.contains(&key) {
                        self.visited_m.insert(key);

                        let new_state = StateTreeNode::new(
                            node.node(), node.offset(), AlignState::Match, Backtrace::ClosedIndel(*v));
                        let new_ix = self.tree.add_node(new_state);

                        Some(new_ix)
                    } else {
                        None
                    }
                },
                _ => None
            }
        }).collect()
    }

    fn generate_next<G>(
        &mut self,
        graph: &G,
        seq_len: usize,
        curr_ix: Ix
    ) -> Vec<(u8, Ix)>
    where
        G: AlignableGraph<NodeIndex=N>,
        Ix: TreeIndexType,
    {
        let mut new_states = Vec::default();
        match self.get_node(curr_ix).state() {
            AlignState::Start | AlignState::Match | AlignState::Mismatch  => {
                // For each successor we can either enter a mismatch state or open a deletion
                for succ in graph.successors(self.get_node(curr_ix).node()) {
                    // Mismatch, only if there's still query sequence to match
                    if self.get_node(curr_ix).offset().as_usize() < seq_len {
                        let key_mismatch = (succ, self.get_node(curr_ix).offset().increase_one());
                        if !self.visited_m.contains(&key_mismatch) {
                            self.visited_m.insert(key_mismatch);

                            let new_state = StateTreeNode::new(
                                succ, self.get_node(curr_ix).offset().increase_one(),
                                AlignState::Mismatch, Backtrace::SingleStep(curr_ix));
                            let new_ix = self.tree.add_node(new_state);
                            new_states.push((self.costs.mismatch(), new_ix));
                        }
                    }

                    // Open deletion
                    let key_del = (succ, self.get_node(curr_ix).offset());
                    if !self.visited_d.contains(&key_del) {
                        self.visited_d.insert(key_del);

                        let new_state = StateTreeNode::new(
                            succ, self.get_node(curr_ix).offset(), AlignState::Deletion, Backtrace::SingleStep(curr_ix));
                        let new_ix = self.tree.add_node(new_state);
                        new_states.push((self.costs.gap_open() + self.costs.gap_extend(), new_ix));
                    }
                }

                if self.get_node(curr_ix).offset().as_usize() < seq_len {
                    // Open insertion
                    let key_ins = (self.get_node(curr_ix).node(), self.get_node(curr_ix).offset().increase_one());
                    if !self.visited_i.contains(&key_ins) {
                        self.visited_i.insert(key_ins);

                        let new_state = StateTreeNode::new(
                            self.get_node(curr_ix).node(), self.get_node(curr_ix).offset().increase_one(),
                            AlignState::Insertion, Backtrace::SingleStep(curr_ix));
                        let new_ix = self.tree.add_node(new_state);
                        new_states.push((self.costs.gap_open() + self.costs.gap_extend(), new_ix));
                    }
                }
            },
            AlignState::Deletion => {
                // Extend deletion for each successor
                for succ in graph.successors(self.get_node(curr_ix).node()) {
                    let key_del = (succ, self.get_node(curr_ix).offset());
                    if !self.visited_d.contains(&key_del) {
                        self.visited_d.insert(key_del);

                        let new_state = StateTreeNode::new(
                            succ, self.get_node(curr_ix).offset(),
                            AlignState::Deletion, Backtrace::SingleStep(curr_ix));
                        let new_ix = self.tree.add_node(new_state);
                        new_states.push((self.costs.cost_gap_extend, new_ix));
                    }
                }
            },
            AlignState::Insertion => {
                if self.get_node(curr_ix).offset().as_usize() < seq_len {
                    let key_ins = (self.get_node(curr_ix).node(), self.get_node(curr_ix).offset().increase_one());
                    if !self.visited_i.contains(&key_ins) {
                        self.visited_i.insert(key_ins);

                        let new_state = StateTreeNode::new(
                            self.get_node(curr_ix).node(), self.get_node(curr_ix).offset().increase_one(),
                            AlignState::Insertion, Backtrace::SingleStep(curr_ix));
                        let new_ix = self.tree.add_node(new_state);
                        new_states.push((self.costs.cost_gap_extend, new_ix));
                    }
                }
            }
            _ => panic!("Invalid alignment state for gap affine!")
        };

        new_states
    }
}

impl<N, O, Ix> Display for GapAffineStateTree<N, O, Ix>
    where
        N: NodeIndexType,
        O: OffsetType,
        Ix: TreeIndexType,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write_tree_gml(f, self)?;

        Ok(())
    }
}
