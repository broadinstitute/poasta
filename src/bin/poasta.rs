//! POASTA command-line driver.

use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Write};
use std::ops::Bound;
use std::path::Path;
use std::process::ExitCode;

use clap::Parser;
use tracing_subscriber::{EnvFilter, fmt, prelude::*};

use poasta::align::AlignmentMode;
use poasta::align::PoastaAligner;
use poasta::align::cost_models::AlignmentCostModel;
use poasta::align::cost_models::affine::Affine;
use poasta::align::cost_models::linear::Linear;
use poasta::align::cost_models::two_piece::TwoPieceAffine;
use poasta::align::traits::AlignOutput;
use poasta::align::astar::{AstarEngine, AstarEngineTwoPiece, DijkstraHeuristic};
use poasta::align::traits::AlignmentEngine;
use poasta::cli::CostModelKind;
use poasta::cli::poasta::{AlignArgs, CliArgs, CliSubcommand, HeuristicKind, OutputType};
use poasta::debug::subscriber::AlignStateLayer;
use poasta::debug::write_graph_dot;
use poasta::errors::PoastaIOError;
use poasta::graph::io::fasta::{FastaOutputOptions, load_graph_from_fasta_msa, poa_graph_to_fasta};
use poasta::graph::io::gfa::{GfaOutputOptions, poa_graph_to_gfa};
use poasta::graph::io::seq::open_sequences;
use poasta::graph::poa::POAGraph;
use poasta::graph::traits::GraphBase;

fn main() -> ExitCode {
    let args = CliArgs::parse();

    let debug_dir = match &args.command {
        Some(CliSubcommand::Align(a)) => a.debug_output_dir.clone(),
        _ => None,
    };
    if let Some(dir) = &debug_dir {
        if let Err(e) = fs::create_dir_all(dir) {
            eprintln!("error: could not create debug output directory: {e}");
            return ExitCode::FAILURE;
        }
    }
    init_tracing(args.log_level.as_tracing_str(), debug_dir.as_deref());

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

/// Install the log formatter on stderr, filtered by `RUST_LOG` or
/// `default_level`, and, when `debug_dir` is given, the A* state dumper.
/// The level filter is attached to the formatter only, so trace-level A*
/// events reach the dumper without being printed.
fn init_tracing(default_level: &str, debug_dir: Option<&Path>) {
    let env_filter =
        EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(default_level));
    let fmt_layer = fmt::layer().with_writer(io::stderr).with_filter(env_filter);
    let debug_layer = debug_dir.map(AlignStateLayer::new);

    if cfg!(feature = "disable_tracing") && debug_dir.is_some() {
        eprintln!("warning: built with `disable_tracing`; --debug-output-dir produces no state files");
    }

    tracing_subscriber::registry()
        .with(fmt_layer)
        .with(debug_layer)
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
    tracing::info!(
        cost_model = ?args.cost_model,
        heuristic = ?args.heuristic,
        input = %args.sequences.display(),
        seed_graph = ?args.graph.as_ref().map(|p| p.display().to_string()),
        "configuring aligner"
    );

    let graph = match &args.graph {
        Some(path) => {
            let reader = BufReader::new(open_seed_graph(path)?);
            let g = load_graph_from_fasta_msa::<u32, _>(reader)?;
            tracing::info!(
                nodes = g.node_count(),
                sequences = g.sequences.len(),
                "loaded seed graph"
            );
            g
        }
        None => POAGraph::<u32>::new(),
    };

    let sequences: Vec<(String, Vec<u8>)> =
        open_sequences(&args.sequences)?.collect::<Result<Vec<_>, _>>()?;
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

/// Translate the CLI's `--alignment-span` preset into the engine-level
/// `AlignmentMode`. `SemiGlobal` frees the graph's leading/trailing indels
/// only (query is aligned end-to-end).
fn alignment_mode(args: &AlignArgs) -> AlignmentMode {
    use poasta::cli::AlignmentSpan;

    match args.alignment_span {
        AlignmentSpan::Global => AlignmentMode::Global,
        AlignmentSpan::SemiGlobal => AlignmentMode::EndsFree {
            qry_free_begin: Bound::Excluded(0),
            qry_free_end: Bound::Excluded(0),
            graph_free_begin: Bound::Unbounded,
            graph_free_end: Bound::Unbounded,
        },
    }
}

fn dispatch(
    args: &AlignArgs,
    graph: POAGraph<u32>,
    sequences: &[(String, Vec<u8>)],
) -> Result<POAGraph<u32>, CliError> {
    let mode = alignment_mode(args);

    match args.cost_model {
        CostModelKind::Affine => {
            let costs = Affine::new(args.cost_mismatch, args.cost_gap_open, args.cost_gap_extend);
            run_astar(costs, mode, args, graph, sequences)
        }
        CostModelKind::Linear => {
            let costs = Linear::new(args.cost_mismatch, args.cost_gap_extend);
            run_astar(costs, mode, args, graph, sequences)
        }
        CostModelKind::TwoPiece => {
            let costs = TwoPieceAffine::new(
                args.cost_mismatch,
                args.cost_gap_open,
                args.cost_gap_extend,
                args.cost_gap_open2,
                args.cost_gap_extend2,
            );
            let engine = AstarEngineTwoPiece::<POAGraph<u32>>::new(costs).with_mode(mode);
            match args.heuristic {
                HeuristicKind::MinGapCost => run_with(engine, args, graph, sequences),
                HeuristicKind::Dijkstra => {
                    run_with(engine.with_heuristic(DijkstraHeuristic), args, graph, sequences)
                }
            }
        }
    }
}

fn run_astar<C>(
    costs: C,
    mode: AlignmentMode,
    args: &AlignArgs,
    graph: POAGraph<u32>,
    sequences: &[(String, Vec<u8>)],
) -> Result<POAGraph<u32>, CliError>
where
    C: AlignmentCostModel + Copy + 'static,
{
    let engine = AstarEngine::<POAGraph<u32>, C>::new(costs).with_mode(mode);
    match args.heuristic {
        HeuristicKind::MinGapCost => run_with(engine, args, graph, sequences),
        HeuristicKind::Dijkstra => {
            run_with(engine.with_heuristic(DijkstraHeuristic), args, graph, sequences)
        }
    }
}

fn run_with<E>(
    engine: E,
    args: &AlignArgs,
    graph: POAGraph<u32>,
    sequences: &[(String, Vec<u8>)],
) -> Result<POAGraph<u32>, CliError>
where
    E: for<'s> AlignmentEngine<
            &'s [u8],
            Graph = POAGraph<u32>,
            Success = AlignOutput<POAGraph<u32>>,
        >,
    for<'s> <E as AlignmentEngine<&'s [u8]>>::Error: std::fmt::Display,
{
    let mut aligner = PoastaAligner::new(engine, graph);

    for (name, seq) in sequences {
        if let Some(dir) = &args.debug_output_dir {
            write_graph_dot(dir, name, aligner.graph(), seq)?;
        }

        aligner
            .align_one_named(name, seq)
            .map_err(|e| CliError::Align(format!("{e}")))?;
    }

    let run_stats = aligner.run_stats();
    tracing::info!(
        n_alignments = run_stats.n_alignments,
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

    writer
        .flush()
        .map_err(|source| PoastaIOError::FileWriteError { source })?;
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
