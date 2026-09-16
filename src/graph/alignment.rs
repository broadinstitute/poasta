//! The hook through which an alignment result is folded into a graph.
use std::error::Error;

pub trait AddAlignment<T> {
    type Error: Error;

    fn add_alignment(
        &mut self,
        sequence_name: &str,
        sequence: &[u8],
        alignment: Option<&T>,
        weights: &[u32],
    ) -> Result<(), Self::Error>;
}
