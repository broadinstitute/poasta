use std::fs::File;
use std::path::PathBuf;
use std::io::{stdout, IsTerminal, Write};

use clap::{Parser, Subcommand, Args, ValueEnum};
use noodles::fasta;
use anyhow::{Result, Context};

use petgraph::dot::Dot;
use poasta::debug::DebugOutputWriter;
use poasta::debug::messages::DebugOutputMessage;

use poasta::graphs::poa::POAGraph;
use poasta::aligner::PoastaAligner;
use poasta::aligner::scoring::gap_affine::GapAffine;
use poasta::errors::PoastaError;
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

fn align(align_args: &AlignArgs) -> Result<()> {
    let debug_writer = align_args.debug_output.as_ref().map(|v| {
        DebugOutputWriter::new(v)
    });

    let mut graph = if let Some(path) = &align_args.graph {
        let file_in = File::open(path)?;
        load_graph(&file_in)?
    } else {
        POAGraph::new()
    };

    // TODO: make configurable through CLI
    let scoring = GapAffine::new(4, 2, 6);
    let mut aligner: PoastaAligner<GapAffine> = if let Some(ref debug) = debug_writer {
        PoastaAligner::new_with_debug_output(scoring, debug)
    } else {
        PoastaAligner::new(scoring)
    };

    // Let's read the sequences from the given FASTA
    let mut reader = fasta::reader::Builder::default().build_from_path(&align_args.sequences)
        .with_context(|| "Could not read FASTA file!".to_string())?;

    for result in reader.records() {
        let record = result?;
        let weights: Vec<usize> = vec![1; record.sequence().len()];
        eprintln!("Aligning {}", record.name());

        if let Some(ref debug) = debug_writer {
            debug.log(DebugOutputMessage::NewSequence {
                seq_name: record.name().to_string(),
                seq_length: record.sequence().len(),
                max_rank: graph.max_rank()
            });

            if !graph.is_empty() {
                debug.log(DebugOutputMessage::new_from_graph(&graph));
            }
        }

        if graph.is_empty() {
            graph.add_alignment_with_weights(record.name(), record.sequence(), None, &weights)?;
        } else {
            let alignment = aligner.align::<u32, u32, _, _, _>(&graph, record.sequence());
            graph.add_alignment_with_weights(record.name(), record.sequence(), Some(&alignment), &weights)?;
        }
    }

    // Determine where to write the graph to
    let mut writer: Box<dyn Output> = if let Some(path) = &align_args.output {
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
            let transformed = graph.graph.map(
                |ix, data|
                    format!("{:?} ({:?})", char::from(data.symbol), graph.get_node_rank(ix)),
                |_, data|
                    format!("{}, {:?}", data.weight, data.sequence_ids)
            );

            let dot = Dot::new(&transformed);
            write!(writer, "{:?}", dot)?
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
            align(v)?
        },
        None => {
            return Err(PoastaError::Other).with_context(|| "No subcommand given.".to_string())
        }
    };

    Ok(())
}
