use crate::align::traits::{AlignResult, AlignableGraph};
use crate::graph::traits::GraphNodeId;

pub mod band_doubling;
pub mod dp;

// Re-export so existing code in engine sub-modules can keep using `AlignmentStats`
// from this module path.
pub use crate::align::traits::AlignmentStats;

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
