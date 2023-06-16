extern crate clap;
extern crate noodles;
extern crate anyhow;

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

}

fn align(align_args: &AlignArgs) -> Result<()> {
    let mut graph = if let Some(path) = &align_args.graph {
        let file_in = File::open(path)?;
        load_graph(&file_in)?
    } else {
        POAGraph::new()
    };

    let mut aligner: WavefrontPOAligner<WFComputeGapAffine<u32>> = WavefrontPOAligner::new("output");

    // Let's read the sequences from the given FASTA
    let mut reader = fasta::reader::Builder::default().build_from_path(&align_args.sequences)
        .with_context(|| "Could not read FASTA file!".to_string())?;

    for result in reader.records() {
        let record = result?;
        let weights: Vec<usize> = vec![1; record.sequence().len()];
        eprintln!("Aligning {}", record.name());

        if graph.is_empty() {
            graph.add_alignment_with_weights(record.name(), record.sequence(), None, &weights)?;
        } else {
            let alignment = aligner.align(&graph, record.sequence());
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


// fn main() -> Result<(), String> {
//     let mut poa_graph = POAGraph::new();
//
//     let seq1 = b"AATGGTTGTCACGTCAGT";
//     let weights1 = vec![1; seq1.len()];
//
//     let seq2 = b"ATTGTAAAGTCTCGTCGGT";
//     let weights2 = vec![1; seq2.len()];
//
//     // EMBOSS_001         1 AATGGT--TGTCACGTCAGT     18
//     //                       ||.||  .|||.||||.||
//     // EMBOSS_001         1 -ATTGTAAAGTCTCGTCGGT     19
//
//     let alignment = vec![
//         AlignedPair{ rpos: Some(NodeIndex::from(0)), qpos: None},
//         AlignedPair{ rpos: Some(NodeIndex::from(1)), qpos: Some(0)},
//         AlignedPair{ rpos: Some(NodeIndex::from(2)), qpos: Some(1)},
//         AlignedPair{ rpos: Some(NodeIndex::from(3)), qpos: Some(2)},
//         AlignedPair{ rpos: Some(NodeIndex::from(4)), qpos: Some(3)},
//         AlignedPair{ rpos: Some(NodeIndex::from(5)), qpos: Some(4)},
//         AlignedPair{ rpos: None, qpos: Some(5)},
//         AlignedPair{ rpos: None, qpos: Some(6)},
//         AlignedPair{ rpos: Some(NodeIndex::from(6)), qpos: Some(7)},
//         AlignedPair{ rpos: Some(NodeIndex::from(7)), qpos: Some(8)},
//         AlignedPair{ rpos: Some(NodeIndex::from(8)), qpos: Some(9)},
//         AlignedPair{ rpos: Some(NodeIndex::from(9)), qpos: Some(10)},
//         AlignedPair{ rpos: Some(NodeIndex::from(10)), qpos: Some(11)},
//         AlignedPair{ rpos: Some(NodeIndex::from(11)), qpos: Some(12)},
//         AlignedPair{ rpos: Some(NodeIndex::from(12)), qpos: Some(13)},
//         AlignedPair{ rpos: Some(NodeIndex::from(13)), qpos: Some(14)},
//         AlignedPair{ rpos: Some(NodeIndex::from(14)), qpos: Some(15)},
//         AlignedPair{ rpos: Some(NodeIndex::from(15)), qpos: Some(16)},
//         AlignedPair{ rpos: Some(NodeIndex::from(16)), qpos: Some(17)},
//         AlignedPair{ rpos: Some(NodeIndex::from(17)), qpos: Some(18)},
//         AlignedPair{ rpos: Some(NodeIndex::from(18)), qpos: None},
//         AlignedPair{ rpos: Some(NodeIndex::from(19)), qpos: None},
//     ];
//
//     poa_graph.add_alignment_with_weights(seq1, None, &weights1)?;
//     poa_graph.add_alignment_with_weights(seq2, Some(&alignment), &weights2)?;
//
//     let transformed = poa_graph.graph.map(
//         |ix, data|
//             format!("{:?} ({:?})", char::from(data.symbol), poa_graph.get_node_rank(ix)),
//         |_, data|
//             format!("{}, {:?}", data.weight, data.sequence_ids)
//     );
//
//     let dot = Dot::new(&transformed);
//     println!("{:?}", dot);
//
//     //        let seq1 = b"AATGGTTGTCACGT------CAGT";
//     let seq3 = b"TTGTCAACATCAGTA";
//
//     let mut aligner: WavefrontPOAligner<u32, WFComputeGapAffine<u32>> = WavefrontPOAligner::new(&poa_graph, "output");
//     let alignment = aligner.align(seq3);
//
//     for AlignedPair{ rpos, qpos} in alignment.iter() {
//         let rpos_str = rpos.map(|v| format!("{}", poa_graph.get_node_rank(v))).unwrap_or(String::new());
//         let qpos_str = qpos.map(|v| format!("{}", v)).unwrap_or(String::new());
//
//         let rnuc = rpos.map(|nix| char::from(poa_graph.graph[nix].symbol)).unwrap_or('|');
//         let qnuc = qpos.map(|p| char::from(seq3[p])).unwrap_or('|');
//
//         let aln_char = if rpos.is_some() && qpos.is_some() {
//             if rnuc == qnuc { '-' } else { '*' }
//         } else {
//             '-'
//         };
//
//         eprintln!("{rpos_str:>5} {rnuc} {aln_char} {qnuc} {qpos_str:<5}");
//     }
//
//     Ok(())
// }
