use std::fs;
use std::fs::File;
use std::path::PathBuf;
use std::io::{stdout, IsTerminal, Write};

use clap::{Parser, Subcommand, Args, ValueEnum};
use noodles::fasta;
use anyhow::{Result, Context};

use petgraph::dot::Dot;

use poasta::graph::POAGraph;
use poasta::wavefront::aligner::WavefrontPOAligner;
use poasta::wavefront::compute::gap_affine::WFComputeGapAffine;
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
    if let Some(debug_output_dir) = &align_args.debug_output {
        fs::create_dir_all(debug_output_dir)?;
    }

    let mut graph = if let Some(path) = &align_args.graph {
        let file_in = File::open(path)?;
        load_graph(&file_in)?
    } else {
        POAGraph::new()
    };

    let mut aligner: WavefrontPOAligner<WFComputeGapAffine<u32>> = if let Some(debug_output_dir) = &align_args.debug_output {
        WavefrontPOAligner::new_with_debug_output(debug_output_dir)
    } else {
        WavefrontPOAligner::new()
    };

    // Let's read the sequences from the given FASTA
    let mut reader = fasta::reader::Builder::default().build_from_path(&align_args.sequences)
        .with_context(|| "Could not read FASTA file!".to_string())?;

    for (i, result) in reader.records().enumerate() {
        let record = result?;
        let weights: Vec<usize> = vec![1; record.sequence().len()];
        eprintln!("Aligning {}", record.name());

        if graph.is_empty() {
            graph.add_alignment_with_weights(record.name(), record.sequence(), None, &weights)?;
        } else {
            let alignment = aligner.align(&graph, record.sequence());
            graph.add_alignment_with_weights(record.name(), record.sequence(), Some(&alignment), &weights)?;
        }

        if let Some(debug_output_dir) = &align_args.debug_output {
            let transformed = graph.graph.map(
                |ix, data|
                    format!("{:?} ({:?})", char::from(data.symbol), graph.get_node_rank(ix)),
                |_, data|
                    format!("{}, {:?}", data.weight, data.sequence_ids)
            );

            let fname = debug_output_dir.join(format!("graph{i}.dot"));
            let dot = Dot::new(&transformed);
            let mut dot_file = File::create(fname)?;
            write!(dot_file, "{:?}", dot)?
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
