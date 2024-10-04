use std::error::Error;
use std::fs::File;
use std::io::{self, BufWriter, IsTerminal};
use std::io::BufReader;
use std::path::Path;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use noodles::fasta;

use poasta::aligner::utils::print_alignment;
use poasta::aligner::{AlignmentMode, GraphAligner};
use poasta::graph::io::graph_to_dot;
use poasta::graph::poa::IndexType;
use tracing::{span, info};
use tracing::{trace_span, Level};
use tracing::Subscriber;
use tracing_subscriber::registry::LookupSpan;
use tracing_subscriber::Registry;
use tracing_subscriber::prelude::*;

use poasta::aligner::astar::heuristic::Dijkstra;
use poasta::aligner::cost_models::affine::Affine;
use poasta::graph::poa::POASeqGraph;
use poasta::errors::PoastaError;

mod cli;
mod debug_tracing;

use debug_tracing::filter::{align_state_filter, only_align_states};
use debug_tracing::subscriber::AlignStateLayer;

/// Any object that supports writing and checking if it is a terminal.
trait PoastaWrite: io::Write + io::IsTerminal {}
impl<T> PoastaWrite for T where T: io::Write + io::IsTerminal {}

/// Build our base tracing subscriber with stderr logging.
fn build_base_subscriber() -> impl Subscriber + for<'span> LookupSpan<'span> {
    let stderr_log = tracing_subscriber::fmt::layer()
        .with_target(false)
        .with_file(false)
        .with_writer(io::stderr)
        .with_ansi(!io::stderr().is_terminal())
        .with_filter(align_state_filter());
    
    Registry::default()
        .with(stderr_log)
}

fn main() -> Result<(), Box<dyn Error + 'static>> {
    let args = cli::CliArgs::parse();
    
    match &args.command {
        Some(cli::CliSubcommand::Align(v)) => align_subcommand(v)?,
        Some(cli::CliSubcommand::View(v)) => (),
        Some(cli::CliSubcommand::Stats(v)) => (),
        None => {
            eprintln!("No subcommand given!");
            
            Err(PoastaError::Other)?
        }
    };

    Ok(())
}

#[cfg(feature = "debug_output")]
fn configure_debug_output(align_args: &cli::AlignArgs) -> Option<AlignStateLayer> {
    align_args.debug_output.as_deref()
        .map(|dir| {
            AlignStateLayer::new(dir)
        })
}

#[cfg(not(feature = "debug_output"))]
fn configure_debug_output(_: &cli::AlignArgs) -> Option<AlignStateLayer> {
    None
}

fn align_subcommand(align_args: &cli::AlignArgs) -> Result<(), Box<dyn Error + 'static>> {
    let span = span!(Level::INFO, "align_subcommand");
    let _enter = span.enter();
    
    build_base_subscriber()
        .with(configure_debug_output(align_args)
            .map(|l| l.with_filter(only_align_states()))
        )
        .init();
    
    // TODO: separate implementations for edit distance/linear gap penalties
    let scoring = Affine::new(
        align_args.cost_mismatch.unwrap_or(4),
        align_args.cost_gap_extend.unwrap_or(2),
        align_args.cost_gap_open.unwrap_or(6),
    );
    
    let mut graph = poasta::graph::poa::POASeqGraph::<u32>::new();
    let aligner = poasta::aligner::PoastaAligner::<Affine, Dijkstra, _, i32, u32>::new(scoring);
    
    perform_alignment(align_args, &mut graph, &aligner, &align_args.sequences)?;
    
    Ok(())
}

fn perform_alignment<Ix, A>(
    align_args: &cli::AlignArgs,
    graph: &mut POASeqGraph<Ix>,
    aligner: &A,
    sequences_fname: &Path,
) -> Result<(), Box<dyn Error + 'static>>
where
    Ix: IndexType,
    A: GraphAligner<POASeqGraph<Ix>>,
{
    // Let's read the sequences from the given FASTA
    let is_gzipped = sequences_fname
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
    let reader_inner: Box<dyn std::io::BufRead> = if is_gzipped {
        Box::new(
            File::open(sequences_fname)
                .map(MultiGzDecoder::new)
                .map(BufReader::new)?,
        )
    } else {
        Box::new(File::open(sequences_fname).map(BufReader::new)?)
    };
    let mut reader = fasta::Reader::new(reader_inner);

    let mut i = 1;
    for result in reader.records() {
        let record = result?;
        let weights: Vec<usize> = vec![1; record.sequence().len()];
        let seq_name = std::str::from_utf8(record.name()).unwrap();
        
        let span = trace_span!(target: "poasta::aligner", "align_seq", seq_id=seq_name);
        let _enter = span.enter();
        
        
        if graph.is_empty() {
            info!("Creating initial graph from {}...", seq_name);
            graph.add_aligned_sequence(seq_name, record.sequence(), &weights, None)?;
        } else {
            info!("Aligning #{i} {}... ", seq_name);
            
            #[cfg(feature = "debug_output")]
            if let Some(debug_dir) = &align_args.debug_output {
                let fname = debug_dir.join(format!("graph_for_{}.dot", seq_name));
                let mut file = File::create(&fname)
                    .map(BufWriter::new)?;
                graph_to_dot(&mut file, graph)?;
            }
            
            let result = aligner.align(graph, record.sequence(), AlignmentMode::Global)?;
            
            info!("Done. Alignment Score: {:?}", result.score);
            info!(
                "{}",
                print_alignment(graph, record.sequence().as_ref(), &result.alignment)
            );

            graph.add_aligned_sequence(
                seq_name,
                record.sequence(),
                &weights,
                Some(&result.alignment),
            )?;
        }

        i += 1;
    }

    Ok(())
}