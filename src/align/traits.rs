use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};
use std::error::Error;

pub trait AlignableGraph: GraphWithNodeOrdering<NodeType = Self::Node> {
    type Node: GraphNodeId; // Mostly here to constrain subtrait associated types

    fn get_node_symbol(&self, p: Self::Node) -> u8;
}

/// Per-alignment counters emitted by the alignment engines.
///
/// `cells_computed` counts each DP state as one cell. The canonical DP
/// touches `(m + 1) * n_real * n_states` cells, so
/// `fraction_of_full_matrix = cells_computed / ((m+1) * n_real * n_states)`
/// (1.0 ⇒ equivalent to the full DP).
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignmentStats {
    pub max_bandwidth: usize,
    pub cells_computed: usize,
    pub fraction_of_full_matrix: f64,
}

/// Result produced by an alignment engine after a successful alignment.
pub trait AlignResult {
    fn alignment_score(&self) -> u32;
    fn alignment_stats(&self) -> AlignmentStats;
}

pub trait AlignmentEngine<ToAlign> {
    type Graph: AlignableGraph;
    type Success: AlignResult;
    type Error: Error;

    fn align(&self, graph: &Self::Graph, to_align: ToAlign) -> Result<Self::Success, Self::Error>;
}
