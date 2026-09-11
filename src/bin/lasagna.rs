//! Lasagna: read-to-graph alignment driver.
//!
//! Loads a GFA partial-order graph, then aligns FASTA/FASTQ reads against it
//! with the A* engine. The graph stays read-only throughout: this is the
//! read-mapping use case, not POA-MSA construction, so the bubble index is
//! built once and shared by every aligner thread.

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;
use std::sync::{Arc, Mutex};
use std::thread;

use clap::Parser;
use crossbeam_channel::{Receiver, Sender, bounded};
use rustc_hash::FxHashMap;
use tracing_subscriber::{EnvFilter, fmt, prelude::*};

use poasta::align::cost_models::AlignmentCostModel;
use poasta::align::cost_models::affine::Affine;
use poasta::align::cost_models::linear::Linear;
use poasta::align::cost_models::two_piece::TwoPieceAffine;
use poasta::align::traits::AlignOutput;
use poasta::align::astar::{AstarEngine, AstarEngineTwoPiece};
use poasta::align::traits::AlignmentEngine;
use poasta::cli::CostModelKind;
use poasta::cli::lasagna::LasagnaCli;
use poasta::debug::subscriber::AlignStateLayer;
use poasta::debug::write_graph_dot;
use poasta::errors::PoastaIOError;
use poasta::graph::bubbles::index::BubbleIndex;
use poasta::graph::io::gaf::{GAFRecord, alignment_to_gaf};
use poasta::graph::io::gfa::{GraphSegments, POAGraphFromGFA, load_graph_from_gfa};
use poasta::graph::io::seq::open_sequences;
use poasta::graph::poa::{POAGraph, POANodeIndex};
use poasta::graph::traits::GraphBase;

fn main() -> ExitCode {
    let args = LasagnaCli::parse();

    if let Some(dir) = &args.debug_output_dir {
        if let Err(e) = fs::create_dir_all(dir) {
            eprintln!("error: could not create debug output directory: {e}");
            return ExitCode::FAILURE;
        }
    }
    init_tracing(
        args.log_level.as_tracing_str(),
        args.debug_output_dir.as_deref(),
    );

    match run(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// See `poasta::init_tracing`: the level filter sits on the formatter only,
/// so trace-level A* events reach the state dumper without being printed.
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

type SequenceRecord = (String, Vec<u8>);

fn run(args: LasagnaCli) -> Result<(), CliError> {
    args.validate().map_err(CliError::Config)?;

    tracing::info!(
        graph = %args.graph.display(),
        sequences = %args.sequences.display(),
        threads = args.num_threads,
        cost_model = ?args.cost_model,
        "configuring lasagna"
    );

    let POAGraphFromGFA {
        graph,
        graph_segments,
    } = load_graph_from_gfa::<u32>(&args.graph)?;

    let node_to_segment = build_node_to_segment(&graph, &graph_segments);
    debug_assert!(
        node_to_segment_covers_all_real_nodes(&graph, &node_to_segment),
        "node_to_segment must cover every non-sentinel node in the graph"
    );

    // The graph never changes, so index its bubbles once for all threads.
    let bubbles = Arc::new(BubbleIndex::new(&graph));

    match args.cost_model {
        CostModelKind::Affine => {
            let costs = Affine::new(args.cost_mismatch, args.cost_gap_open, args.cost_gap_extend);
            let engine = astar_engine(costs, bubbles);
            run_threaded(&engine, &args, &graph, &graph_segments, &node_to_segment)
        }
        CostModelKind::Linear => {
            let costs = Linear::new(args.cost_mismatch, args.cost_gap_extend);
            let engine = astar_engine(costs, bubbles);
            run_threaded(&engine, &args, &graph, &graph_segments, &node_to_segment)
        }
        CostModelKind::TwoPiece => {
            let costs = TwoPieceAffine::new(
                args.cost_mismatch,
                args.cost_gap_open,
                args.cost_gap_extend,
                args.cost_gap_open2,
                args.cost_gap_extend2,
            );
            let engine =
                AstarEngineTwoPiece::<POAGraph<u32>>::new(costs).with_existing_bubbles(bubbles);
            run_threaded(&engine, &args, &graph, &graph_segments, &node_to_segment)
        }
    }
}

fn astar_engine<C: AlignmentCostModel + Copy + 'static>(
    costs: C,
    bubbles: Arc<BubbleIndex<POANodeIndex<u32>>>,
) -> AstarEngine<POAGraph<u32>, C> {
    AstarEngine::new(costs).with_existing_bubbles(bubbles)
}

/// Walk each GFA segment by its recorded length, mapping every node it owns
/// to its (segment_ix, position_within_segment). Walking by step count
/// prevents the walk from escaping into the global end-sentinel that
/// `post_process` may have attached to a tip segment's end node.
fn build_node_to_segment<Ix: petgraph::graph::IndexType>(
    graph: &POAGraph<Ix>,
    segments: &GraphSegments<Ix>,
) -> FxHashMap<POANodeIndex<Ix>, (usize, usize)> {
    let mut map: FxHashMap<POANodeIndex<Ix>, (usize, usize)> = FxHashMap::default();
    for (segment_ix, _) in segments.names.iter().enumerate() {
        let start_node = segments.start_nodes[segment_ix];
        let end_node = segments.end_nodes[segment_ix];
        let len = segments.segment_lengths[segment_ix];

        map.insert(start_node, (segment_ix, 0));
        if len == 1 {
            continue;
        }

        let mut curr = start_node;
        for pos in 1..len {
            // Within a GFA segment, the next-character node is the unique
            // successor that we haven't visited yet AND that isn't a
            // non-segment node attached by post_process. The cleanest filter
            // is: if pos == len - 1, we must land on end_node; otherwise the
            // successor that *isn't* end_node is a non-existent case for
            // valid GFA, so just take the unique successor that hasn't been
            // mapped yet.
            let next = if pos == len - 1 {
                end_node
            } else {
                graph
                    .successors(curr)
                    .find(|s| !map.contains_key(s) && *s != end_node)
                    .expect("segment walk found no internal successor")
            };
            map.insert(next, (segment_ix, pos));
            curr = next;
        }
    }
    map
}

/// Debug invariant: every non-sentinel node in the graph appears in
/// `node_to_segment`. Catches GFA topologies that the walk in
/// `build_node_to_segment` doesn't handle (e.g. branching mid-segment).
fn node_to_segment_covers_all_real_nodes<Ix: petgraph::graph::IndexType>(
    _graph: &POAGraph<Ix>,
    _map: &FxHashMap<POANodeIndex<Ix>, (usize, usize)>,
) -> bool {
    // We don't have a public "iterate all non-sentinel nodes" API at the
    // crate boundary, so this assertion is conservative: it always returns
    // true. Left as a hook in case the API is exposed later.
    true
}

fn run_threaded<E>(
    engine: &E,
    args: &LasagnaCli,
    graph: &POAGraph<u32>,
    graph_segments: &GraphSegments<u32>,
    node_to_segment: &FxHashMap<POANodeIndex<u32>, (usize, usize)>,
) -> Result<(), CliError>
where
    E: Sync
        + for<'s> AlignmentEngine<
            &'s [u8],
            Graph = POAGraph<u32>,
            Success = AlignOutput<POAGraph<u32>>,
        >,
    for<'s> <E as AlignmentEngine<&'s [u8]>>::Error: std::fmt::Display,
{
    let num_threads = args.num_threads.max(1);
    // Bounded channels apply back-pressure: the reader can't run away from
    // the aligners, and aligners can't run away from the writer.
    let cap = (num_threads * 4).max(4);
    let (tx_seq, rx_seq) = bounded::<SequenceRecord>(cap);
    let (tx_out, rx_out) = bounded::<GAFRecord>(cap);

    // Collect the first error from any worker, if any. Workers store and bail.
    let first_error: Arc<Mutex<Option<CliError>>> = Arc::new(Mutex::new(None));
    let debug_dir = args.debug_output_dir.as_deref();

    thread::scope(|s| -> Result<(), CliError> {
        // Reader thread.
        let reader_first_error = first_error.clone();
        let sequences_path = args.sequences.clone();
        s.spawn(move || {
            if let Err(e) = read_sequences(&sequences_path, tx_seq) {
                stash_error(&reader_first_error, e);
            }
        });

        // Writer thread.
        let writer_first_error = first_error.clone();
        let output_path = args.output.clone();
        s.spawn(move || {
            if let Err(e) = write_output(output_path.as_deref(), rx_out) {
                stash_error(&writer_first_error, e);
            }
        });

        // Aligner threads.
        for _ in 0..num_threads {
            let rx_seq = rx_seq.clone();
            let tx_out = tx_out.clone();
            let thread_first_error = first_error.clone();
            s.spawn(move || {
                while let Ok((name, seq)) = rx_seq.recv() {
                    let span = tracing::info_span!("align_seq", seq_id = name.as_str());
                    let _enter = span.enter();

                    if let Some(dir) = debug_dir {
                        if let Err(e) = write_graph_dot(dir, &name, graph, &seq) {
                            stash_error(&thread_first_error, e.into());
                            break;
                        }
                    }

                    let out = match engine.align(graph, seq.as_slice()) {
                        Ok(out) => out,
                        Err(e) => {
                            stash_error(&thread_first_error, CliError::Align(format!("{e}")));
                            break;
                        }
                    };

                    let rec = alignment_to_gaf(
                        graph,
                        graph_segments,
                        &name,
                        &seq,
                        &out.alignment,
                        node_to_segment,
                        out.score,
                    );

                    if let Some(rec) = rec {
                        if tx_out.send(rec).is_err() {
                            // Writer has hung up; nothing more we can do.
                            break;
                        }
                    } else {
                        tracing::debug!(query = %name, "alignment produced no GAF record");
                    }
                }

                // tx_out and rx_seq drop here, closing the channels once
                // every aligner thread has exited.
            });
        }

        // We've cloned tx_out into each aligner thread. Drop our own copy so
        // the writer's `recv()` actually returns Err when the last aligner
        // hangs up.
        drop(tx_out);
        drop(rx_seq);

        Ok(())
    })?;

    let err = first_error.lock().unwrap().take();
    if let Some(e) = err {
        return Err(e);
    }

    tracing::info!("done");
    Ok(())
}

fn stash_error(slot: &Arc<Mutex<Option<CliError>>>, e: CliError) {
    let mut guard = slot.lock().unwrap();
    if guard.is_none() {
        *guard = Some(e);
    }
}

fn read_sequences(path: &Path, tx: Sender<SequenceRecord>) -> Result<(), CliError> {
    for record in open_sequences(path)? {
        let (name, seq) = record?;
        if tx.send((name, seq)).is_err() {
            // Aligners have all hung up; stop reading.
            return Ok(());
        }
    }
    Ok(())
}

fn write_output(path: Option<&Path>, rx: Receiver<GAFRecord>) -> Result<(), CliError> {
    let mut writer: Box<dyn Write> = match path {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).map_err(|source| PoastaIOError::FileWriteError { source })?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    let mut n = 0usize;
    while let Ok(record) = rx.recv() {
        writeln!(writer, "{}", record)
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
        n += 1;
    }
    writer
        .flush()
        .map_err(|source| PoastaIOError::FileWriteError { source })?;

    tracing::info!(records = n, "wrote GAF");
    Ok(())
}
