pub mod heuristic;

/// Enum representing the alignment state of a particular cell in the alignment matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignState {
    Match,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
}


pub trait AstarAlignableGraph {
    type NodeType;
    
    fn node_count(&self) -> usize;
}


pub trait AstarState {
    type QueueItem;
}