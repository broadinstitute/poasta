//! POASTA: partial order alignment with an A* search.
//!
//! - [`graph`]: the POA graph, its traits, consensus, superbubbles and I/O.
//! - [`align`]: aligning a query to a graph, with the engine in [`align::astar`].
//! - [`cli`] and [`debug`] (feature `cli`): argument structs and the
//!   `--debug-output-dir` tracing layer used by the binaries.
//!
//! See `docs/architecture.md` for the full walkthrough.

pub mod align;
#[cfg(feature = "cli")]
pub mod cli;
#[cfg(feature = "cli")]
pub mod debug;
pub mod errors;
pub mod graph;
