//! The partial order alignment graph and everything that works on it:
//! [`poa::POAGraph`] itself, the [`traits`] engines see it through,
//! [`alignment::AddAlignment`] for folding alignments in, [`consensus`],
//! superbubble detection in [`bubbles`], and file formats in [`io`].
//!
//! Nothing here depends on the aligner, except that `poa.rs` implements
//! `AddAlignment` for the engine's `AlignOutput`.

pub mod alignment;
pub mod bubbles;
pub mod consensus;
pub mod io;
pub mod poa;
pub mod traits;
pub mod utils;

#[cfg(test)]
pub(crate) mod mock;
