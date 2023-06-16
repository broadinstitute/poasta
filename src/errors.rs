use std::error::Error;
use std::fmt::{Display, Formatter};
use std::io;
use petgraph::algo::Cycle;

#[derive(Debug)]
pub enum PoastaError {
    /// The size of the weights vector is not equal to the length of the sequence
    WeightsUnequalSize(usize, usize),

    /// The alignment did not include any valid positions
    InvalidAlignment,

    /// Error variant when something went wrong in the alignment
    AlignmentError,

    /// Error indicating the graph is in an invalid state
    GraphError,

    /// Error variant when we couldn't read from a file
    FileReadError { source: io::Error },

    /// Error variant when we could not serialize the graph to binary representation
    SerializationError { source: bincode::Error },

    /// Other IO errors
    IOError(io::Error),

    /// Other miscellaneous poasta errors
    Other,
}

impl Error for PoastaError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match *self {
            Self::FileReadError { ref source } => Some(source),
            Self::SerializationError { ref source } => Some(source),
            Self::IOError(ref source) => Some(source),
            _ => None
        }
    }
}


impl<N> From<Cycle<N>> for PoastaError {
    fn from(_: Cycle<N>) -> Self {
        Self::GraphError
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


impl Display for PoastaError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::WeightsUnequalSize(seq_len, weights_len) =>
                write!(f, "The length of the weights vector ({weights_len}) does not match the length of the sequence ({seq_len})!"),
            Self::InvalidAlignment =>
                write!(f, "The specified alignment did not include any valid sequence positions!"),
            Self::GraphError =>
                write!(f, "The graph is in an invalid state (possibly a cycle?)."),
            Self::AlignmentError =>
                write!(f, "Something went wrong with the alignment!"),
            Self::FileReadError { source: _ } =>
                write!(f, "Could not read from file!"),
            Self::SerializationError { source: _ } =>
                write!(f, "Could not serialize the graph to file!"),
            Self::IOError(ref err) =>
                err.fmt(f),
            Self::Other =>
                write!(f, "Poasta error!")
        }
    }
}
