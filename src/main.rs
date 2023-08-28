use std::fs;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::{stdout, IsTerminal, Write, BufReader};

use clap::{Parser, Subcommand, Args, ValueEnum};
use noodles::fasta;
use anyhow::{Result, Context};

use petgraph::graph::IndexType;
use serde::de::DeserializeOwned;
use flate2::read::GzDecoder;

use poasta::aligner::alignment::print_alignment;
use poasta::debug::DebugOutputWriter;
use poasta::debug::messages::DebugOutputMessage;

use poasta::graphs::poa::{POAGraph, POAGraphWithIx};
use poasta::aligner::PoastaAligner;
use poasta::aligner::scoring::{AlignmentCosts, GapAffine};
use poasta::errors::PoastaError;
use poasta::io::graph::load_graph_from_fasta_msa;
use poasta::io::load_graph;


trait Output: Write + IsTerminal { }
impl<T> Output for T where T: Write + IsTerminal { }


/// The various output formats supported by Poasta
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum OutputType {
    /// Default Poasta graph file format
    POASTA,

    /// Output a tabular MSA in FASTA file format
    FASTA,

    /// Output the graph as GFA
    GFA,

    /// Output the graph in DOT format for visualization
    DOT,
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct CliArgs {
    /// Set verbosity level. Use multiple times to increase the verbosity level.
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    #[command(subcommand)]
    command: Option<CliSubcommand>,
}


#[derive(Subcommand, Debug)]
enum CliSubcommand {
    /// Perform multiple sequence alignment and create or update POA graphs
    Align(AlignArgs)
}

#[derive(Args, Debug)]
struct AlignArgs {
    /// Sequences to align in FASTA format.
    sequences: PathBuf,

    /// Input partial order graph to align sequences to. If not specified, will create a new graph from input sequences.
    #[arg(short, long)]
    graph: Option<PathBuf>,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output file type.
    #[arg(value_enum, short='O', long)]
    output_type: Option<OutputType>,

    /// Output debug information (intermediate graphs, aligner state) and write files to the given directory
    #[arg(short, long)]
    debug_output: Option<PathBuf>,

}

fn perform_alignment<Ix, C>(
    graph: &mut POAGraph<Ix>,
    aligner: &mut PoastaAligner<C>,
    debug_writer: Option<&DebugOutputWriter>,
    sequences_fname: &Path
) -> Result<()>
where
    Ix: IndexType + DeserializeOwned,
    C: AlignmentCosts,
{
    // Let's read the sequences from the given FASTA
    let is_gzipped = sequences_fname.file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
    let reader_inner: Box<dyn std::io::BufRead> = if is_gzipped {
        Box::new(File::open(sequences_fname)
            .map(GzDecoder::new)
            .map(BufReader::new)?)
    } else {
        Box::new(File::open(sequences_fname)
            .map(BufReader::new)?)
    };
    let mut reader = fasta::Reader::new(reader_inner);

    let mut i = 1;
    for result in reader.records() {
        let record = result?;
        let weights: Vec<usize> = vec![1; record.sequence().len()];

        if let Some(debug) = debug_writer {
            debug.log(DebugOutputMessage::NewSequence {
                seq_name: record.name().to_string(),
                seq_length: record.sequence().len(),
                max_rank: graph.max_rank()
            });

            if !graph.is_empty() {
                debug.log(DebugOutputMessage::new_from_graph(graph));
            }
        }

        if graph.is_empty() {
            // eprintln!("Creating initial graph from {}...", record.name());
            graph.add_alignment_with_weights(record.name(), record.sequence(), None, &weights)?;
        } else {
            eprint!("Aligning #{i} {}... ", record.name());
            let (score, alignment) = aligner.align::<u32, usize, _, _, _>(graph, record.sequence());
            eprintln!("Done. Alignment Score: {:?}", score);
            // eprintln!();
            // eprintln!("{}", print_alignment(graph, record.sequence(), &alignment));
            // eprintln!();

            graph.add_alignment_with_weights(record.name(), record.sequence(), Some(&alignment), &weights)?;
        }

        i += 1;
    }

    Ok(())
}

fn align_subcommand(align_args: &AlignArgs) -> Result<()> {
    let debug_writer = align_args.debug_output.as_ref().map(|v| {
        DebugOutputWriter::new(v)
    });

    let mut graph = if let Some(path) = &align_args.graph {
        let fasta_extensions = vec![".fa", ".fa.gz", ".fna", ".fna.gz", ".fasta", ".fasta.gz"];
        let path_as_str = path.to_string_lossy();
        eprintln!("Ext: {:?}, {:?}", path, path_as_str.ends_with(".fa"));
        if fasta_extensions.into_iter().any(|ext| path_as_str.ends_with(ext)) {
            load_graph_from_fasta_msa(path)?
        } else {
            let file_in = File::open(path)?;
            load_graph(&file_in)?
        }
    } else {
        POAGraphWithIx::USIZE(POAGraph::new())
    };

    // TODO: make configurable through CLI
    let scoring = GapAffine::new(4, 2, 6);
    let mut aligner: PoastaAligner<GapAffine> = if let Some(ref debug) = debug_writer {
        PoastaAligner::new_with_debug_output(scoring, debug)
    } else {
        PoastaAligner::new(scoring)
    };

    match graph {
        POAGraphWithIx::U8(ref mut g) =>
            perform_alignment(g, &mut aligner, debug_writer.as_ref(), align_args.sequences.as_ref())?,
        POAGraphWithIx::U16(ref mut g) =>
            perform_alignment(g, &mut aligner, debug_writer.as_ref(), align_args.sequences.as_ref())?,
        POAGraphWithIx::U32(ref mut g) =>
            perform_alignment(g, &mut aligner, debug_writer.as_ref(), align_args.sequences.as_ref())?,
        POAGraphWithIx::USIZE(ref mut g) =>
            perform_alignment(g, &mut aligner, debug_writer.as_ref(), align_args.sequences.as_ref())?,
    }

    // Determine where to write the graph to
    let mut writer: Box<dyn Output> = if let Some(path) = &align_args.output {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?
        }

        let file = File::create(path)?;
        Box::new(file) as Box<dyn Output>
    } else {
        Box::new(stdout()) as Box<dyn Output>
    };

    let output_type = align_args.output_type.unwrap_or(OutputType::POASTA);
    match output_type {
        OutputType::POASTA => {
            if !writer.is_terminal() {
                poasta::io::save_graph(&graph, writer)?
            } else {
                eprintln!("WARNING: not writing binary graph data to terminal standard output!");
            }
        },
        OutputType::DOT => {
            write!(writer, "{}", &graph)?
        }
        _ => return Err(PoastaError::Other).with_context(|| "Other output formats not supported yet!".to_string())
    }

    if let Some(debug) = debug_writer {
        eprintln!("Waiting for debug writer thread to finish...");
        debug.log(DebugOutputMessage::Terminate);
        debug.join()?;
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = CliArgs::parse();

    match &args.command {
        Some(CliSubcommand::Align(v)) => {
            align_subcommand(v)?
        },
        None => {
            return Err(PoastaError::Other).with_context(|| "No subcommand given.".to_string())
        }
    };

    Ok(())
}
