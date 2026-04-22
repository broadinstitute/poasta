//! POASTA command-line driver.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

use clap::Parser;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

use poasta::align::cost_models::affine::Affine;
use poasta::align::cost_models::linear::Linear;
use poasta::align::cost_models::two_piece::TwoPieceAffine;
use poasta::align::engine::band_doubling::BandDoublingEngineScalar;
use poasta::align::engine::dp::CanonicalDP;
use poasta::align::engine::AlignResult;
use poasta::align::traits::AlignmentEngine;
use poasta::align::PoastaAligner;
use poasta::cli::poasta::{
    AlignArgs, CliArgs, CliSubcommand, CostModelKind, EngineKind, OutputType,
};
use poasta::errors::PoastaIOError;
use poasta::graph::io::fasta::{
    load_graph_from_fasta_msa, poa_graph_to_fasta, FastaOutputOptions,
};
use poasta::graph::io::gfa::{poa_graph_to_gfa, GfaOutputOptions};
use poasta::graph::io::seq::open_sequences;
use poasta::graph::poa::POAGraph;
use poasta::graph::traits::GraphBase;

fn main() -> ExitCode {
    let args = CliArgs::parse();

    init_tracing(args.log_level.as_tracing_str());

    match args.command {
        Some(CliSubcommand::Align(align_args)) => match run_align(align_args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::FAILURE
            }
        },
        Some(CliSubcommand::View(_)) => {
            eprintln!("error: `view` subcommand is not yet implemented on this branch");
            ExitCode::FAILURE
        }
        Some(CliSubcommand::Stats(_)) => {
            eprintln!("error: `stats` subcommand is not yet implemented on this branch");
            ExitCode::FAILURE
        }
        None => {
            eprintln!("error: no subcommand given (try `poasta align --help`)");
            ExitCode::FAILURE
        }
    }
}

fn init_tracing(default_level: &str) {
    let filter =
        EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(default_level));
    tracing_subscriber::registry()
        .with(filter)
        .with(fmt::layer().with_writer(io::stderr))
        .init();
}

#[derive(Debug)]
enum CliError {
    Io(PoastaIOError),
    Config(String),
    Align(String),
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Config(m) => write!(f, "{m}"),
            Self::Align(m) => write!(f, "alignment failed: {m}"),
        }
    }
}

impl From<PoastaIOError> for CliError {
    fn from(e: PoastaIOError) -> Self {
        Self::Io(e)
    }
}

impl From<io::Error> for CliError {
    fn from(e: io::Error) -> Self {
        Self::Io(PoastaIOError::OtherError { source: e })
    }
}

fn run_align(args: AlignArgs) -> Result<(), CliError> {
    args.validate().map_err(CliError::Config)?;

    tracing::info!(
        engine = ?args.engine,
        cost_model = ?args.cost_model,
        initial_k = args.initial_k,
        input = %args.sequences.display(),
        seed_graph = ?args.graph.as_ref().map(|p| p.display().to_string()),
        "configuring aligner"
    );

    let graph = match &args.graph {
        Some(path) => {
            let reader = BufReader::new(open_seed_graph(path)?);
            let g = load_graph_from_fasta_msa::<u32, _>(reader)?;
            tracing::info!(nodes = g.node_count(), sequences = g.sequences.len(), "loaded seed graph");
            g
        }
        None => POAGraph::<u32>::new(),
    };

    let sequences: Vec<(String, Vec<u8>)> = open_sequences(&args.sequences)?
        .collect::<Result<Vec<_>, _>>()?;
    tracing::info!(count = sequences.len(), "read input sequences");

    let graph = dispatch(&args, graph, &sequences)?;

    write_output(&args, &graph)?;

    tracing::info!(
        nodes = graph.node_count(),
        sequences = graph.sequences.len(),
        "done"
    );
    Ok(())
}

fn open_seed_graph(path: &Path) -> Result<Box<dyn io::Read>, PoastaIOError> {
    let file = File::open(path).map_err(|source| PoastaIOError::FileReadError { source })?;
    let gzipped = path
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);
    if gzipped {
        Ok(Box::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn dispatch(
    args: &AlignArgs,
    graph: POAGraph<u32>,
    sequences: &[(String, Vec<u8>)],
) -> Result<POAGraph<u32>, CliError> {
    match args.cost_model {
        CostModelKind::Affine => {
            let costs = Affine::new(
                args.cost_match,
                args.cost_mismatch,
                args.cost_gap_open,
                args.cost_gap_extend,
            );
            match args.engine {
                EngineKind::BandDoubling => {
                    let engine = BandDoublingEngineScalar::<Affine, u32>::new(costs)
                        .with_initial_k(args.initial_k);
                    run_with(engine, graph, sequences)
                }
                EngineKind::CanonicalDp => {
                    let engine = CanonicalDP::<Affine, POAGraph<u32>>::new(costs);
                    run_with(engine, graph, sequences)
                }
            }
        }
        CostModelKind::Linear => {
            let costs = Linear::new(args.cost_match, args.cost_mismatch, args.cost_gap_extend);
            match args.engine {
                EngineKind::BandDoubling => {
                    let engine = BandDoublingEngineScalar::<Linear, u32>::new(costs)
                        .with_initial_k(args.initial_k);
                    run_with(engine, graph, sequences)
                }
                EngineKind::CanonicalDp => {
                    let engine = CanonicalDP::<Linear, POAGraph<u32>>::new(costs);
                    run_with(engine, graph, sequences)
                }
            }
        }
        CostModelKind::TwoPiece => {
            let costs = TwoPieceAffine::new(
                args.cost_match,
                args.cost_mismatch,
                args.cost_gap_open,
                args.cost_gap_extend,
                args.cost_gap_open2,
                args.cost_gap_extend2,
            );
            match args.engine {
                EngineKind::BandDoubling => {
                    let engine = BandDoublingEngineScalar::<TwoPieceAffine, u32>::new(costs)
                        .with_initial_k(args.initial_k);
                    run_with(engine, graph, sequences)
                }
                EngineKind::CanonicalDp => {
                    let engine = CanonicalDP::<TwoPieceAffine, POAGraph<u32>>::new(costs);
                    run_with(engine, graph, sequences)
                }
            }
        }
    }
}

fn run_with<E>(
    engine: E,
    graph: POAGraph<u32>,
    sequences: &[(String, Vec<u8>)],
) -> Result<POAGraph<u32>, CliError>
where
    E: for<'s> AlignmentEngine<&'s [u8], Graph = POAGraph<u32>, Success = AlignResult<POAGraph<u32>>>,
    for<'s> <E as AlignmentEngine<&'s [u8]>>::Error: std::fmt::Display,
{
    let mut aligner = PoastaAligner::new(engine, graph);
    let iter = sequences
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_slice()));
    aligner
        .align_all_named(iter)
        .map_err(|e| CliError::Align(format!("{e}")))?;
    let run_stats = aligner.run_stats();
    tracing::info!(
        n_alignments = run_stats.n_alignments,
        avg_max_bandwidth = run_stats.avg_max_bandwidth(),
        max_bandwidth_overall = run_stats.max_bandwidth_overall,
        avg_cells_computed = run_stats.avg_cells_computed(),
        avg_fraction_of_full_matrix = run_stats.avg_fraction(),
        "run alignment stats (averaged across sequences)",
    );
    Ok(aligner.into_graph())
}

fn write_output(args: &AlignArgs, graph: &POAGraph<u32>) -> Result<(), CliError> {
    let output_type = resolve_output_type(args);

    if args.consensus_only && matches!(output_type, OutputType::Gfa) {
        return Err(CliError::Config(
            "--consensus-only is incompatible with GFA output; \
             rerun without -O gfa or without --consensus-only"
                .into(),
        ));
    }

    let mut writer: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).map_err(|source| PoastaIOError::FileWriteError { source })?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    tracing::info!(
        ?output_type,
        output = ?args.output.as_ref().map(|p| p.display().to_string()),
        include_consensus = args.include_consensus,
        consensus_only = args.consensus_only,
        "writing output"
    );

    match output_type {
        OutputType::Fasta => {
            let opts = FastaOutputOptions {
                include_consensus: args.include_consensus,
                consensus_only: args.consensus_only,
            };
            poa_graph_to_fasta(graph, &mut writer, opts)?;
        }
        OutputType::Gfa => {
            let opts = GfaOutputOptions {
                include_consensus: args.include_consensus,
            };
            poa_graph_to_gfa(graph, &mut writer, opts)?;
        }
    }

    writer.flush().map_err(|source| PoastaIOError::FileWriteError { source })?;
    Ok(())
}

fn resolve_output_type(args: &AlignArgs) -> OutputType {
    if let Some(t) = args.output_type {
        return t;
    }
    if let Some(path) = &args.output {
        let name = path
            .file_name()
            .map(|v| v.to_string_lossy().to_lowercase())
            .unwrap_or_default();
        let stripped = name.strip_suffix(".gz").unwrap_or(&name);
        if stripped.ends_with(".gfa") {
            return OutputType::Gfa;
        }
        if stripped.ends_with(".fa")
            || stripped.ends_with(".fasta")
            || stripped.ends_with(".fna")
            || stripped.ends_with(".msa")
        {
            return OutputType::Fasta;
        }
    }
    OutputType::Fasta
}
