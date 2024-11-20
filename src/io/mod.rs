pub mod fasta;
pub mod graph;

pub use graph::load_graph_from_fasta_msa;
#[cfg(feature = "serialize")]
pub use graph::{load_graph, save_graph};
