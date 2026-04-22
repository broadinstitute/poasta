use crate::align::traits::AlignableGraph;
use crate::graph::traits::GraphNodeId;

pub mod band_doubling;
pub mod dp;

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

pub struct AlignResult<G: AlignableGraph> {
    pub score: u32,
    pub alignment: Vec<AlignedPair<G::Node>>,
    pub stats: AlignmentStats,
}
