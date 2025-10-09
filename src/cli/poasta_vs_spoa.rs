use clap::Parser;

use std::path::PathBuf;


#[derive(Parser)]
pub struct PoastaVsSpoaArgs {
    /// Sequences to align in FASTA format
    #[clap(help_heading = "Input/Output")]
    pub sequences: PathBuf,

    /// Working directory and output directory for files
    #[arg(short, long)]
    #[clap(help_heading = "Input/Output")]
    pub output_dir: PathBuf,

    /// Number of sequences to use for graph construction. The remaining
    /// sequences will be used for testing.
    #[arg(short, long, default_value = "25")]
    #[clap(help_heading = "Input")]
    pub num_graph: usize,

    /// Penalty for mismatching bases
    #[arg(short = 'n', default_value = "4")]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_mismatch: Option<u8>,

    /// Penalty for opening a new gap
    #[arg(short = 'g', default_value = "6")]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_gap_open: Option<u8>,

    /// Penalty for extending a gap
    #[arg(short = 'e', default_value = "2")]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_gap_extend: Option<u8>,
}
