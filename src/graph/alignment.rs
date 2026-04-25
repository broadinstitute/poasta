//! Contains types for representing sequence-to-graph alignments.
use std::error::Error;

use petgraph::graph::IndexType;

use super::poa::POANodeIndex;

pub trait AddAlignment<T> {
    type Error: Error;

    fn add_alignment(
        &mut self,
        sequence_name: &str,
        sequence: &[u8],
        alignment: Option<&T>,
        weights: &[usize],
    ) -> Result<(), Self::Error>;
}

/// Represents an alignment pairing between a node in the graph and a symbol
/// in the aligned sequence.
#[derive(Clone, Copy, Debug)]
pub struct SeqAlignedPair<Ix>
where
    Ix: IndexType,
{
    pub rpos: Option<POANodeIndex<Ix>>,
    pub qpos: Option<usize>,
}

/// The alignment of a single sequence to a partial order graph
pub type SeqAlignment<Ix> = Vec<SeqAlignedPair<Ix>>;
