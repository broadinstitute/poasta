//! Lasagna CLI surface: flat parser, no subcommand.

use std::path::PathBuf;

use clap::Parser;

use crate::cli::{AlignmentSpan, CostModelKind, LogLevel};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Align reads (FASTA/FASTQ) to an existing GFA partial-order graph, emit GAF."
)]
pub struct LasagnaCli {
    /// Log verbosity level. Overridden by the RUST_LOG env var when set.
    #[arg(long, value_enum, default_value = "info")]
    pub log_level: LogLevel,

    /// Input partial order graph in GFA format (optionally gzipped). The graph
    /// is loaded read-only. Lasagna does not modify it.
    #[clap(help_heading = "Inputs")]
    pub graph: PathBuf,

    /// Sequences to align in FASTA or FASTQ format (optionally gzipped).
    #[clap(help_heading = "Inputs")]
    pub sequences: PathBuf,

    /// Number of aligner worker threads.
    #[arg(short = 'j', long, default_value_t = 1)]
    #[clap(help_heading = "Processing")]
    pub num_threads: usize,

    /// Output filename. If not given, defaults to stdout. Output is always
    /// GAF. The filename is not used to choose the format.
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    pub output: Option<PathBuf>,

    /// Alignment span. Only `global` is currently supported on lasagna.
    #[arg(short = 'm', long, value_enum, default_value = "global")]
    #[clap(help_heading = "Alignment configuration")]
    pub alignment_span: AlignmentSpan,

    /// Which cost model to use.
    #[arg(long, value_enum, default_value = "affine")]
    #[clap(help_heading = "Cost model")]
    pub cost_model: CostModelKind,

    /// Penalty for mismatching bases.
    #[arg(short = 'n', long, default_value_t = 4)]
    #[clap(help_heading = "Cost model")]
    pub cost_mismatch: u8,

    /// Penalty for opening a new gap. Ignored when --cost-model=linear.
    #[arg(short = 'g', long, default_value_t = 6)]
    #[clap(help_heading = "Cost model")]
    pub cost_gap_open: u8,

    /// Penalty for extending a gap.
    #[arg(short = 'e', long, default_value_t = 2)]
    #[clap(help_heading = "Cost model")]
    pub cost_gap_extend: u8,

    /// Penalty for opening a second-piece gap. Only used with --cost-model=two-piece.
    #[arg(short = 'G', long, default_value_t = 24)]
    #[clap(help_heading = "Cost model")]
    pub cost_gap_open2: u8,

    /// Penalty for extending a second-piece gap. Only used with --cost-model=two-piece.
    #[arg(short = 'E', long, default_value_t = 1)]
    #[clap(help_heading = "Cost model")]
    pub cost_gap_extend2: u8,

    /// Write per-query debug files to this directory: `graph_for_<query>.dot`,
    /// `<query>_visited.tsv` and `<query>_queued.tsv`. Filenames are derived
    /// from query names.
    #[arg(long)]
    #[clap(help_heading = "Debug")]
    pub debug_output_dir: Option<PathBuf>,
}

impl LasagnaCli {
    /// Sanity-check the combination of CLI flags. Mirrors the pattern in
    /// `poasta::cli::poasta::AlignArgs::validate`.
    pub fn validate(&self) -> Result<(), String> {
        if !matches!(self.alignment_span, AlignmentSpan::Global) {
            return Err(format!(
                "alignment span '{:?}' is not yet supported on lasagna; only 'global' is available",
                self.alignment_span
            ));
        }
        Ok(())
    }
}
