use clap::Parser;

use std::path::PathBuf;

use crate::cli::poasta::CostModelKind;

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

    /// Which cost model to use for POASTA alignment.
    /// Note: SPOA always uses affine gap costs for graph construction and scoring.
    #[arg(long, value_enum, default_value = "affine")]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_model: CostModelKind,

    /// Penalty for mismatching bases.
    #[arg(short = 'n', long, default_value_t = 4)]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_mismatch: u8,

    /// Penalty for opening a new gap. Ignored when --cost-model=linear.
    #[arg(short = 'g', long, default_value_t = 6)]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_gap_open: u8,

    /// Penalty for extending a gap.
    #[arg(short = 'e', long, default_value_t = 2)]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_gap_extend: u8,

    /// Penalty for opening a second-piece gap. Only used with --cost-model=two-piece.
    #[arg(short = 'G', long, default_value_t = 24)]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_gap_open2: u8,

    /// Penalty for extending a second-piece gap. Only used with --cost-model=two-piece.
    #[arg(short = 'E', long, default_value_t = 1)]
    #[clap(help_heading = "Alignment configuration")]
    pub cost_gap_extend2: u8,
}
