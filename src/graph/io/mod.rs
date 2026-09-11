//! Graph and sequence I/O: FASTA and FASTQ input, FASTA MSA and GFA graphs
//! in both directions, GAF for read alignments, DOT for debugging.

pub mod dot;
pub mod fasta;
pub mod gaf;
pub mod gfa;
pub mod seq;

pub const CONSENSUS_NAME: &str = "consensus";
