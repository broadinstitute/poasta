//! Contains types for representing sequence-to-graph alignments.
use std::error::Error;

use petgraph::graph::IndexType;

use super::poa::POANodeIndex;


pub trait AddAlignment<T> {
    type Error: Error;
    
    fn add_alignment(&mut self, 
        sequence_name: &str, 
        sequence: &[u8], 
        alignment: Option<&T>, 
        weights: &[usize]
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





/// Refers to a specific position with a POA graph node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct POANodePos<Ix>(pub POANodeIndex<Ix>, pub usize)
where
    Ix: IndexType;

impl<Ix> POANodePos<Ix>
where
    Ix: IndexType,
{
    #[inline]
    pub fn node_equal(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<Ix> AlignableGraphNodePos for POANodePos<Ix>
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;

    #[inline]
    fn new(node: Self::NodeType, pos: usize) -> Self {
        POANodePos(node, pos)
    }

    #[inline]
    fn node(&self) -> Self::NodeType {
        self.0
    }

    #[inline]
    fn pos(&self) -> usize {
        self.1
    }
}
