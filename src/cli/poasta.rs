use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};

use crate::cli::{AlignmentSpan, CostModelKind, LogLevel};

/// The various output formats supported by Poasta
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum OutputType {
    /// Output a tabular MSA in FASTA file format
    Fasta,

    /// Output the graph as GFA
    Gfa,
}

/// Which alignment engine the CLI should dispatch to.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum EngineKind {
    /// Priority-queue-driven A* aligner (default). The only engine that supports
    /// `--alignment-span` ends-free and semi-global alignment.
    Astar,

    /// Scalar Ukkonen-style band-doubling Gotoh aligner (default).
    BandDoubling,

    /// Canonical O(N*m) Gotoh DP (reference / oracle — slow on large inputs).
    CanonicalDp,
}

/// Which A* heuristic to use. Ignored by non-`astar` engines.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum HeuristicKind {
    /// Admissible min-gap-cost lower bound (default): tighter bound, visits
    /// fewer states than Dijkstra.
    MinGapCost,

    /// Plain Dijkstra (`h ≡ 0`): exact, but explores every state within the
    /// optimal-score radius.
    Dijkstra,
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct CliArgs {
    /// Log verbosity level. Overridden by the RUST_LOG env var when set.
    #[arg(long, global = true, value_enum, default_value = "info")]
    pub log_level: LogLevel,

    #[command(subcommand)]
    pub command: Option<CliSubcommand>,
}

#[derive(Subcommand, Debug)]
pub enum CliSubcommand {
    /// Perform multiple sequence alignment and create or update POA graphs
    Align(AlignArgs),

    /// Convert POASTA POA graphs to various output formats
    View(ViewArgs),

    /// Print graph statistics
    Stats(StatsArgs),
}

#[derive(Args, Debug)]
pub struct AlignArgs {
    /// Sequences to align in FASTA or FASTQ format, optionally gzipped.
    #[clap(help_heading = "Inputs")]
    pub sequences: PathBuf,

    /// Input partial order graph (as a FASTA MSA, optionally gzipped) to align sequences to.
    /// If not specified, a new graph is created from the input sequences.
    #[arg(short = 'I', long)]
    #[clap(help_heading = "Inputs")]
    pub graph: Option<PathBuf>,

    /// Output filename. If not given, defaults to stdout.
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    pub output: Option<PathBuf>,

    /// Output file type. When absent, inferred from --output extension (.gfa → gfa, else fasta).
    #[arg(value_enum, short = 'O', long)]
    #[clap(help_heading = "Outputs")]
    pub output_type: Option<OutputType>,

    /// Include a consensus sequence (Lee's heaviest-bundling) alongside the normal output.
    /// For FASTA output, appended as the last record named "consensus".
    /// For GFA output, emitted as an additional P line named "consensus".
    #[arg(long, conflicts_with = "consensus_only")]
    #[clap(help_heading = "Outputs")]
    pub include_consensus: bool,

    /// Output only the consensus sequence as a single FASTA record. Incompatible with GFA output.
    #[arg(long, conflicts_with = "include_consensus")]
    #[clap(help_heading = "Outputs")]
    pub consensus_only: bool,

    /// Alignment span. Only 'global' is currently supported.
    #[arg(short = 'm', long, value_enum, default_value = "global")]
    #[clap(help_heading = "Alignment configuration")]
    pub alignment_span: AlignmentSpan,

    /// Which alignment engine to use.
    #[arg(long, value_enum, default_value = "astar")]
    #[clap(help_heading = "Alignment configuration")]
    pub engine: EngineKind,

    /// Starting bandwidth for band-doubling. Ignored by canonical-dp and astar.
    #[arg(long, default_value_t = 1)]
    #[clap(help_heading = "Alignment configuration")]
    pub initial_k: usize,

    /// Which heuristic the astar engine uses. Ignored by other engines.
    #[arg(long, value_enum, default_value = "min-gap-cost")]
    #[clap(help_heading = "Alignment configuration")]
    pub heuristic: HeuristicKind,

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

    /// Write debug output files (DOT graph, band TSV, cell TSV) to this directory.
    /// One set of files is produced per aligned sequence. Only supported with the
    /// band-doubling engine.
    #[arg(long)]
    #[clap(help_heading = "Debug")]
    pub debug_output_dir: Option<PathBuf>,
}

impl AlignArgs {
    /// Sanity-check the combination of CLI flags.
    pub fn validate(&self) -> Result<(), String> {
        if !matches!(self.alignment_span, AlignmentSpan::Global)
            && !matches!(self.engine, EngineKind::Astar)
        {
            return Err(format!(
                "alignment span '{:?}' requires --engine astar; \
                 band-doubling and canonical-dp only support global alignment",
                self.alignment_span
            ));
        }
        Ok(())
    }
}

#[derive(Args, Debug)]
pub struct StatsArgs {
    /// The POASTA graph or an existing MSA in FASTA format to analyze
    pub graph: PathBuf,
}

#[derive(Args, Debug)]
pub struct ViewArgs {
    /// Input POA graph
    pub graph: PathBuf,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Output file type
    #[arg(value_enum, short = 'O', long)]
    pub output_type: OutputType,
}
