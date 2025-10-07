use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};


/// The various output formats supported by Poasta
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum OutputType {
    /// Default Poasta graph file format
    Poasta,

    /// Output a tabular MSA in FASTA file format
    Fasta,

    /// Output the graph as GFA
    Gfa,

    /// Output the graph in DOT format for visualization
    Dot,
}

/// An enum indicating what kind of alignment to perform
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum AlignmentSpan {
    /// Perform global alignment
    Global,

    /// Perform semi-global alignment, i.e., globally align query but allow free gaps in the graph
    /// at the beginning and end
    SemiGlobal,

    /// Perform ends-free alignment, i.e., indels at the beginning or end on either the query or
    /// graph are free
    EndsFree,
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct CliArgs {
    /// Set verbosity level. Use multiple times to increase the verbosity level.
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,

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
    /// Sequences to align in FASTA format.
    #[clap(help_heading = "Inputs")]
    pub sequences: PathBuf,

    /// Input partial order graph to align sequences to. If not specified,
    /// will create a new graph from input sequences.
    #[arg(short = 'I', long)]
    #[clap(help_heading = "Inputs")]
    pub graph: Option<PathBuf>,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    pub output: Option<PathBuf>,

    /// Output file type.
    #[arg(value_enum, short = 'O', long)]
    #[clap(help_heading = "Outputs")]
    pub output_type: Option<OutputType>,

    /// Output debug information (intermediate graphs, aligner state)
    /// and write files to the given directory
    #[cfg(feature = "debug_output")]
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    pub debug_output: Option<PathBuf>,

    /// Alignment span, either global alignment, semi-global alignment or ends-free alignment.
    #[arg(short = 'm', long, default_value = "global")]
    #[clap(help_heading = "Alignment configuration")]
    pub alignment_span: AlignmentSpan,

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

#[derive(Args, Debug)]
pub struct StatsArgs {
    /// The POASTA graph or an existing MSA in FASTA format to analyze
    graph: PathBuf,
}

#[derive(Args, Debug)]
pub struct ViewArgs {
    /// Input POA graph
    graph: PathBuf,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output file type
    #[arg(value_enum, short = 'O', long)]
    output_type: OutputType,
}
