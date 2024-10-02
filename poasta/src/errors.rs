use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::io;

use petgraph::graph::IndexType;
use petgraph::algo::Cycle;

use crate::graph::poa::POANodeIndex;

// use crate::debug::messages::DebugOutputMessage;

#[derive(Debug)]
pub enum GraphError<Ix>
where 
    Ix: IndexType,
{
    /// Error variant when the graph contains a cycle
    CycleError(Cycle<POANodeIndex<Ix>>),
    
    /// Error when a required alignment is not present
    EmptyAlignment,
}

impl<Ix> From<Cycle<POANodeIndex<Ix>>> for GraphError<Ix> 
where
    Ix: IndexType,
{
    fn from(value: Cycle<POANodeIndex<Ix>>) -> Self {
        Self::CycleError(value)
    }
}

impl<Ix> Error for GraphError<Ix>
where 
    Ix: IndexType 
{ }

impl<Ix> Display for GraphError<Ix> 
where 
    Ix: IndexType,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CycleError(cycle) => write!(f, "The graph contains a cycle: {:?} visited twice.", cycle.node_id()),
            Self::EmptyAlignment => write!(f, "Alignment required when adding a sequence to a non-empty graph."),
        }
    }
}

#[derive(Debug)]
pub enum PoastaIOError {
    FileReadError { source: io::Error },
    FileWriteError { source: io::Error },
    OtherError { source: io::Error }
}

impl Error for PoastaIOError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match *self {
            Self::FileReadError { ref source } => Some(source),
            Self::FileWriteError { ref source } => Some(source),
            Self::OtherError { ref source } => Some(source),
        }
    }
}

impl Display for PoastaIOError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::FileReadError { source } => write!(f, "Error reading file: {}", source),
            Self::FileWriteError { source } => write!(f, "Error writing file: {}", source),
            Self::OtherError { source } => write!(f, "Other IO error: {}", source),
        }
    }
}

impl From<io::Error> for PoastaIOError {
    fn from(value: io::Error) -> Self {
        Self::OtherError { source: value }
    }
}



#[derive(Debug)]
pub enum PoastaError {
    /// The size of the weights vector is not equal to the length of the sequence
    WeightsUnequalSize(usize, usize),

    /// The alignment did not include any valid positions
    InvalidAlignment,

    /// Error variant when something went wrong in the alignment
    AlignmentError,

    /// Other formatting errors
    FormatError(std::fmt::Error),

    /// Other miscellaneous poasta errors
    Other,
}

impl Error for PoastaError { }


impl From<std::fmt::Error> for PoastaError {
    fn from(value: std::fmt::Error) -> Self {
        Self::FormatError(value)
    }
}

impl Display for PoastaError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::WeightsUnequalSize(seq_len, weights_len) =>
                write!(f, "The length of the weights vector ({weights_len}) does not match the length of the sequence ({seq_len})!"),
            Self::InvalidAlignment =>
                write!(f, "The specified alignment did not include any valid sequence positions!"),
            Self::AlignmentError =>
                write!(f, "Something went wrong with the alignment!"),
            Self::FormatError(ref err) =>
                std::fmt::Display::fmt(err, f),
            // Self::DebugError { source: _ } =>
            //     write!(f, "Could not log debug data!"),
            Self::Other =>
                write!(f, "Poasta error!")
        }
    }
}
