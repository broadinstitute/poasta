//! The interface between graphs and alignment engines.
//!
//! [`AlignableGraph`] is what an engine needs from a graph, [`AlignmentEngine`]
//! is what the aligner needs from an engine, and [`AlignOutput`] carries the
//! score, the aligned pairs and the search statistics back.

use std::error::Error;

use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};

pub trait AlignableGraph: GraphWithNodeOrdering<NodeType = Self::Node> {
    type Node: GraphNodeId; // Mostly here to constrain subtrait associated types

    fn get_node_symbol(&self, p: Self::Node) -> u8;

    /// Whether the graph node `p` carries the query byte `symbol`.
    ///
    /// A* compares query bytes against graph symbols in its hottest loop
    /// (the greedy match-extension DFA); this provided default is trivial,
    /// but gives implementations a chance to specialise (e.g. IUPAC
    /// ambiguity codes) without changing the trait's callers.
    #[inline]
    fn is_symbol_equal(&self, p: Self::Node, symbol: u8) -> bool {
        self.get_node_symbol(p) == symbol
    }
}

/// Per-alignment counters emitted by the alignment engine.
///
/// `cells_computed` counts every alignment state the search expanded.
/// `fraction_of_full_matrix` divides that by the size of the full dynamic
/// programming matrix, `(m + 1) * real_nodes * states`, so `1.0` means the
/// search did as much work as a full-matrix DP would have.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignmentStats {
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
/// An aligned pair of residues. The first element represent
/// the position with a node of the graph, and the second element
/// represents the query sequence position.
///
/// In case of on insertion or deletion, set one of the elements to `None`.
#[derive(Debug, Clone, Copy)]
pub struct AlignedPair<N>(pub(crate) Option<N>, pub(crate) Option<usize>);

impl<N> AlignedPair<N>
where
    N: GraphNodeId,
{
    pub fn new(node_pos: Option<N>, query_pos: Option<usize>) -> Self {
        AlignedPair(node_pos, query_pos)
    }

    #[inline(always)]
    pub fn node(&self) -> Option<N> {
        self.0
    }

    #[inline(always)]
    pub fn query_pos(&self) -> Option<usize> {
        self.1
    }

    pub fn is_aligned(&self) -> bool {
        self.0.is_some() && self.1.is_some()
    }

    /// `true` when this pair represents a query insertion relative to the
    /// graph: a query position with no corresponding graph node.
    #[inline(always)]
    pub fn is_insertion(&self) -> bool {
        self.0.is_none() && self.1.is_some()
    }

    /// `true` when this pair represents a deletion in the query relative to
    /// the graph: a graph node with no corresponding query position.
    #[inline(always)]
    pub fn is_deletion(&self) -> bool {
        self.0.is_some() && self.1.is_none()
    }
}

pub struct AlignOutput<G: AlignableGraph> {
    pub score: u32,
    pub alignment: Vec<AlignedPair<G::Node>>,
    pub stats: AlignmentStats,
}

impl<G: AlignableGraph> AlignResult for AlignOutput<G> {
    fn alignment_score(&self) -> u32 {
        self.score
    }

    fn alignment_stats(&self) -> AlignmentStats {
        self.stats
    }
}
