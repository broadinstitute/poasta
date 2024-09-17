use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::io;
use std::sync::mpsc::SendError;
use petgraph::algo::Cycle;

use crate::graph::poa::POANodeIndex;

// use crate::debug::messages::DebugOutputMessage;

#[derive(Debug)]
pub enum GraphError<Ix> {
    /// Error variant when the graph contains a cycle
    CycleError(Cycle<POANodeIndex<Ix>>),
    
    /// Error when a required alignment is not present
    EmptyAlignment,
}

impl<Ix> From<Cycle<POANodeIndex<Ix>>> for GraphError<Ix> {
    fn from(value: Cycle<POANodeIndex<Ix>>) -> Self {
        Self::CycleError(value)
    }
}

impl<Ix> Error for GraphError<Ix>
where Ix: Debug 
{ }

impl<Ix> Display for GraphError<Ix> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CycleError(cycle) => write!(f, "The graph contains a cycle: {:?} visited twice.", cycle.node_id()),
            Self::EmptyAlignment => write!(f, "Alignment required when adding a sequence to a non-empty graph."),
        }
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

    /// Error variant when we couldn't read from a file
    FileReadError { source: io::Error },

    /// Error variant when we could not serialize the graph to binary representation
    SerializationError { source: bincode::Error },

    /// Other IO errors
    IOError(io::Error),

    /// Other formatting errors
    FormatError(std::fmt::Error),

    /// Debug output error
    // DebugError { source: SendError<DebugOutputMessage> },

    /// Other miscellaneous poasta errors
    Other,
}

impl Error for PoastaError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match *self {
            Self::FileReadError { ref source } => Some(source),
            Self::SerializationError { ref source } => Some(source),
            Self::IOError(ref source) => Some(source),
            // Self::DebugError { ref source} => Some(source),
            _ => None
        }
    }
}


impl From<io::Error> for PoastaError {
    fn from(value: io::Error) -> Self {
        Self::IOError(value)
    }
}

impl From<bincode::Error> for PoastaError {
    fn from(value: bincode::Error) -> Self {
        Self::SerializationError {
            source: value
        }
    }
}

// impl From<SendError<DebugOutputMessage>> for PoastaError {
//     fn from(value: SendError<DebugOutputMessage>) -> Self {
//         Self::DebugError { source: value }
//     }
// }

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
            Self::FileReadError { source: _ } =>
                write!(f, "Could not read from file!"),
            Self::SerializationError { source: _ } =>
                write!(f, "Error loading/saving the graph from a POASTA graph file!"),
            Self::IOError(ref err) =>
                std::fmt::Display::fmt(err, f),
            Self::FormatError(ref err) =>
                std::fmt::Display::fmt(err, f),
            // Self::DebugError { source: _ } =>
            //     write!(f, "Could not log debug data!"),
            Self::Other =>
                write!(f, "Poasta error!")
        }
    }
}
