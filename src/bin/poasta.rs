use std::fs;
use std::fs::File;
use std::io::{stdout, BufReader, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::ops::Bound;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use noodles::fasta;

use flate2::read::MultiGzDecoder;
use petgraph::graph::IndexType;
use serde::de::DeserializeOwned;

use poasta::aligner::config::{AffineMinGapCost, Affine2PieceMinGapCost, AlignmentConfig};
use poasta::debug::messages::DebugOutputMessage;
use poasta::debug::DebugOutputWriter;

use poasta::aligner::scoring::{AlignmentType, GapAffine, GapAffine2Piece};
use poasta::aligner::PoastaAligner;
use poasta::errors::PoastaError;
use poasta::graphs::poa::{POAGraph, POAGraphWithIx};
use poasta::graphs::AlignableRefGraph;
use poasta::io::fasta::poa_graph_to_fasta;
use poasta::io::graph::{graph_to_dot, graph_to_gfa, load_graph_from_fasta_msa};
use poasta::io::load_graph;

trait Output: Write + IsTerminal {}
impl<T> Output for T where T: Write + IsTerminal {}

/// The various output formats supported by Poasta
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum OutputType {
    /// Default Poasta graph file format
    Poasta,

    /// Output a tabular MSA in FASTA file format
    Fasta,

    /// Output the graph as GFA
    Gfa,

    /// Output the graph in DOT format for visualization
    Dot,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
/// An enum indicating what kind of alignment to perform
enum AlignmentSpan {
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
struct CliArgs {
    /// Set verbosity level. Use multiple times to increase the verbosity level.
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    #[command(subcommand)]
    command: Option<PoastaSubcommand>,
}

#[derive(Subcommand, Debug)]
enum PoastaSubcommand {
    /// Perform multiple sequence alignment and create or update POA graphs
    Align(AlignArgs),

    /// Convert POASTA POA graphs to various output formats
    View(ViewArgs),

    /// Print graph statistics
    Stats(StatsArgs),
}

#[derive(Args, Debug)]
struct AlignArgs {
    /// Sequences to align in FASTA format.
    #[clap(help_heading = "Inputs")]
    sequences: PathBuf,

    /// Input partial order graph to align sequences to. If not specified,
    /// will create a new graph from input sequences.
    #[arg(short = 'I', long)]
    #[clap(help_heading = "Inputs")]
    graph: Option<PathBuf>,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    output: Option<PathBuf>,

    /// Output file type.
    #[arg(value_enum, short = 'O', long)]
    #[clap(help_heading = "Outputs")]
    output_type: Option<OutputType>,

    /// Output debug information (intermediate graphs, aligner state)
    /// and write files to the given directory
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    debug_output: Option<PathBuf>,

    /// Alignment span, either global alignment, semi-global alignment or ends-free alignment.
    #[arg(short = 'm', long, default_value = "global")]
    #[clap(help_heading = "Alignment configuration")]
    alignment_span: AlignmentSpan,

    /// Penalty for mismatching bases
    #[arg(short = 'n', default_value = "4")]
    #[clap(help_heading = "Alignment configuration")]
    cost_mismatch: Option<u8>,

    /// Penalty for opening a new gap. Use comma-separated values for two-piece affine (e.g., "6,24")
    #[arg(short = 'g', default_value = "6")]
    #[clap(help_heading = "Alignment configuration")]
    cost_gap_open: String,

    /// Penalty for extending a gap. Use comma-separated values for two-piece affine (e.g., "2,1")
    #[arg(short = 'e', default_value = "2")]
    #[clap(help_heading = "Alignment configuration")]
    cost_gap_extend: String,
}

#[derive(Args, Debug)]
struct StatsArgs {
    /// The POASTA graph or an existing MSA in FASTA format to analyze
    graph: PathBuf,
}

#[derive(Args, Debug)]
struct ViewArgs {
    /// Input POA graph
    graph: PathBuf,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output file type
    #[arg(value_enum, short = 'O', long)]
    output_type: OutputType,
}

fn perform_alignment<N, C>(
    graph: &mut POAGraph<N>,
    aligner: &mut PoastaAligner<C>,
    debug_writer: Option<&DebugOutputWriter>,
    sequences_fname: &Path,
) -> Result<()>
where
    N: IndexType + DeserializeOwned,
    C: AlignmentConfig,
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

    for result in reader.records() {
        let record = result?;
        let weights: Vec<usize> = vec![1; record.sequence().len()];
        
        let seq_name = std::str::from_utf8(record.name())?;

        if let Some(debug) = debug_writer {
            debug.log(DebugOutputMessage::NewSequence {
                seq_name: seq_name.to_string(),
                sequence: String::from_utf8_lossy(record.sequence().as_ref()).to_string(),
                max_rank: graph.node_count_with_start_and_end(),
            });

            if !graph.is_empty() {
                debug.log(DebugOutputMessage::new_from_graph(graph));
            }
        }

        if graph.is_empty() {
            // eprintln!("Creating initial graph from {}...", record.name());
            graph.add_alignment_with_weights(seq_name, record.sequence().as_ref(), None, &weights)?;
        } else {
            // eprint!("Aligning #{i} {}... ", seq_name);
            let result = aligner.align::<u32, _>(graph, record.sequence().as_ref());
            // eprintln!("Done. Alignment Score: {:?}", result.score);
            // eprintln!();
            // eprintln!(
            //     "{}",
            //     print_alignment(graph, record.sequence().as_ref(), &result.alignment)
            // );
            // eprintln!();
            // eprintln!();

            graph.add_alignment_with_weights(
                seq_name,
                record.sequence().as_ref(),
                Some(&result.alignment),
                &weights,
            )?;
        }
    }

    Ok(())
}

fn parse_gap_penalties(gap_str: &str) -> Result<Vec<u8>> {
    gap_str
        .split(',')
        .map(|s| s.trim().parse::<u8>().context("Invalid gap penalty value"))
        .collect()
}

fn align_subcommand(align_args: &AlignArgs) -> Result<()> {
    let debug_writer = align_args
        .debug_output
        .as_ref()
        .map(DebugOutputWriter::init);
        
    // Convert alignment span to AlignmentType
    let alignment_type = match align_args.alignment_span {
        AlignmentSpan::Global => AlignmentType::Global,
        AlignmentSpan::SemiGlobal => AlignmentType::EndsFree {
            qry_free_begin: Bound::Unbounded,
            qry_free_end: Bound::Unbounded,
            graph_free_begin: Bound::Unbounded,
            graph_free_end: Bound::Unbounded,
        },
        AlignmentSpan::EndsFree => AlignmentType::EndsFree {
            qry_free_begin: Bound::Unbounded,
            qry_free_end: Bound::Unbounded,
            graph_free_begin: Bound::Unbounded,
            graph_free_end: Bound::Unbounded,
        },
    };

    let mut graph = if let Some(path) = &align_args.graph {
        let fasta_extensions = vec![".fa", ".fa.gz", ".fna", ".fna.gz", ".fasta", ".fasta.gz"];
        let path_as_str = path.to_string_lossy();
        if fasta_extensions
            .into_iter()
            .any(|ext| path_as_str.ends_with(ext))
        {
            load_graph_from_fasta_msa(path)?
        } else {
            let file_in = File::open(path)?;
            load_graph(&file_in)?
        }
    } else {
        POAGraphWithIx::U32(POAGraph::new())
    };

    // Parse gap penalties
    let gap_open_values = parse_gap_penalties(&align_args.cost_gap_open)?;
    let gap_extend_values = parse_gap_penalties(&align_args.cost_gap_extend)?;
    
    // Check if we're using two-piece affine mode
    let use_two_piece = gap_open_values.len() == 2 && gap_extend_values.len() == 2;
    
    if use_two_piece {
        if gap_open_values.len() != 2 || gap_extend_values.len() != 2 {
            return Err(anyhow::anyhow!("Two-piece affine mode requires exactly 2 values for both gap-open and gap-extend (e.g., -g 8,24 -e 2,1)"));
        }
        
        let gap_open1 = gap_open_values[0];
        let gap_extend1 = gap_extend_values[0];
        let gap_open2 = gap_open_values[1];
        let gap_extend2 = gap_extend_values[1];
        
        if gap_extend1 <= gap_extend2 {
            eprintln!("Warning: gap_extend1 ({}) should be greater than gap_extend2 ({}) for two-piece model", gap_extend1, gap_extend2);
            eprintln!("Using standard affine gap model instead.");
            let scoring = GapAffine::new(
                align_args.cost_mismatch.unwrap_or(4),
                gap_extend1,
                gap_open1,
            );
            let mut aligner = if let Some(ref debug) = debug_writer {
                PoastaAligner::new_with_debug(AffineMinGapCost(scoring), alignment_type, debug)
            } else {
                PoastaAligner::new(AffineMinGapCost(scoring), alignment_type)
            };
            
            match graph {
                POAGraphWithIx::U8(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
                POAGraphWithIx::U16(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
                POAGraphWithIx::U32(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
                POAGraphWithIx::USIZE(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
            }
        } else {
            let scoring = GapAffine2Piece::new(
                align_args.cost_mismatch.unwrap_or(4),
                gap_extend1,
                gap_open1,
                gap_extend2,
                gap_open2,
            );
            let mut aligner = if let Some(ref debug) = debug_writer {
                PoastaAligner::new_with_debug(Affine2PieceMinGapCost(scoring), alignment_type, debug)
            } else {
                PoastaAligner::new(Affine2PieceMinGapCost(scoring), alignment_type)
            };
            
            match graph {
                POAGraphWithIx::U8(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
                POAGraphWithIx::U16(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
                POAGraphWithIx::U32(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
                POAGraphWithIx::USIZE(ref mut g) => perform_alignment(
                    g,
                    &mut aligner,
                    debug_writer.as_ref(),
                    align_args.sequences.as_ref(),
                )?,
            }
        }
    } else {
        // Use standard affine gap model
        if gap_open_values.len() != 1 || gap_extend_values.len() != 1 {
            return Err(anyhow::anyhow!("Standard affine mode requires exactly 1 value for both gap-open and gap-extend (e.g., -g 6 -e 2)"));
        }
        
        let scoring = GapAffine::new(
            align_args.cost_mismatch.unwrap_or(4),
            gap_extend_values[0],
            gap_open_values[0],
        );
        let mut aligner = if let Some(ref debug) = debug_writer {
            PoastaAligner::new_with_debug(AffineMinGapCost(scoring), alignment_type, debug)
        } else {
            PoastaAligner::new(AffineMinGapCost(scoring), alignment_type)
        };

        match graph {
            POAGraphWithIx::U8(ref mut g) => perform_alignment(
                g,
                &mut aligner,
                debug_writer.as_ref(),
                align_args.sequences.as_ref(),
            )?,
            POAGraphWithIx::U16(ref mut g) => perform_alignment(
                g,
                &mut aligner,
                debug_writer.as_ref(),
                align_args.sequences.as_ref(),
            )?,
            POAGraphWithIx::U32(ref mut g) => perform_alignment(
                g,
                &mut aligner,
                debug_writer.as_ref(),
                align_args.sequences.as_ref(),
            )?,
            POAGraphWithIx::USIZE(ref mut g) => perform_alignment(
                g,
                &mut aligner,
                debug_writer.as_ref(),
                align_args.sequences.as_ref(),
            )?,
        }
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

    let output_type = align_args.output_type.unwrap_or(OutputType::Poasta);
    match output_type {
        OutputType::Poasta => {
            if !writer.is_terminal() {
                poasta::io::save_graph(&graph, writer)?
            } else {
                eprintln!("WARNING: not writing binary graph data to terminal standard output!");
            }
        }
        OutputType::Dot => write!(writer, "{}", &graph)?,
        OutputType::Fasta => match graph {
            POAGraphWithIx::U8(ref g) => poa_graph_to_fasta(g, &mut writer),
            POAGraphWithIx::U16(ref g) => poa_graph_to_fasta(g, &mut writer),
            POAGraphWithIx::U32(ref g) => poa_graph_to_fasta(g, &mut writer),
            POAGraphWithIx::USIZE(ref g) => poa_graph_to_fasta(g, &mut writer),
        }?,
        OutputType::Gfa => match graph {
            POAGraphWithIx::U8(ref g) => graph_to_gfa(&mut writer, g),
            POAGraphWithIx::U16(ref g) => graph_to_gfa(&mut writer, g),
            POAGraphWithIx::U32(ref g) => graph_to_gfa(&mut writer, g),
            POAGraphWithIx::USIZE(ref g) => graph_to_gfa(&mut writer, g),
        }?,
    }

    if let Some(debug) = debug_writer {
        eprintln!("Waiting for debug writer thread to finish...");
        debug.log(DebugOutputMessage::Terminate);
        debug.join()?;
    }

    Ok(())
}

fn view_subcommand(view_args: &ViewArgs) -> Result<()> {
    let fasta_extensions = vec![".fa", ".fa.gz", ".fna", ".fna.gz", ".fasta", ".fasta.gz"];
    let path_as_str = view_args.graph.to_string_lossy();
    let graph = if fasta_extensions
        .into_iter()
        .any(|ext| path_as_str.ends_with(ext))
    {
        load_graph_from_fasta_msa(&view_args.graph)?
    } else {
        let file_in = File::open(&view_args.graph)?;
        load_graph(&file_in)?
    };

    // Determine where to write the graph to
    let mut writer: Box<dyn Output> = if let Some(path) = &view_args.output {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?
        }

        let file = File::create(path)?;
        Box::new(file) as Box<dyn Output>
    } else {
        Box::new(stdout()) as Box<dyn Output>
    };

    match view_args.output_type {
        OutputType::Poasta => {
            if !writer.is_terminal() {
                poasta::io::save_graph(&graph, writer)?
            } else {
                eprintln!("WARNING: not writing binary graph data to terminal standard output!");
            }
        }
        OutputType::Dot => match graph {
            POAGraphWithIx::U8(ref g) => graph_to_dot(&mut writer, g),
            POAGraphWithIx::U16(ref g) => graph_to_dot(&mut writer, g),
            POAGraphWithIx::U32(ref g) => graph_to_dot(&mut writer, g),
            POAGraphWithIx::USIZE(ref g) => graph_to_dot(&mut writer, g),
        }?,
        OutputType::Fasta => match graph {
            POAGraphWithIx::U8(ref g) => poa_graph_to_fasta(g, &mut writer),
            POAGraphWithIx::U16(ref g) => poa_graph_to_fasta(g, &mut writer),
            POAGraphWithIx::U32(ref g) => poa_graph_to_fasta(g, &mut writer),
            POAGraphWithIx::USIZE(ref g) => poa_graph_to_fasta(g, &mut writer),
        }?,
        OutputType::Gfa => match graph {
            POAGraphWithIx::U8(ref g) => graph_to_gfa(&mut writer, g),
            POAGraphWithIx::U16(ref g) => graph_to_gfa(&mut writer, g),
            POAGraphWithIx::U32(ref g) => graph_to_gfa(&mut writer, g),
            POAGraphWithIx::USIZE(ref g) => graph_to_gfa(&mut writer, g),
        }?,
    }

    Ok(())
}

fn stats_subcommand(stats_args: &StatsArgs) -> Result<()> {
    let fasta_extensions = vec![".fa", ".fa.gz", ".fna", ".fna.gz", ".fasta", ".fasta.gz"];
    let path_as_str = stats_args.graph.to_string_lossy();
    let graph = if fasta_extensions
        .into_iter()
        .any(|ext| path_as_str.ends_with(ext))
    {
        load_graph_from_fasta_msa(&stats_args.graph)?
    } else {
        let file_in = File::open(&stats_args.graph)?;
        load_graph(&file_in)?
    };

    match graph {
        POAGraphWithIx::U8(ref g) => print_graph_stats(g),
        POAGraphWithIx::U16(ref g) => print_graph_stats(g),
        POAGraphWithIx::U32(ref g) => print_graph_stats(g),
        POAGraphWithIx::USIZE(ref g) => print_graph_stats(g),
    }

    Ok(())
}

fn print_graph_stats<G: AlignableRefGraph>(graph: &G) {
    eprintln!("node_count: {}", graph.node_count());
    eprintln!(
        "node_count_with_start: {}",
        graph.node_count_with_start_and_end()
    );
    eprintln!("edge_count: {}", graph.edge_count());

    let in_degrees: Vec<_> = graph.all_nodes().map(|n| graph.in_degree(n)).collect();
    let avg_in_degree = in_degrees.iter().sum::<usize>() as f64 / in_degrees.len() as f64;
    let out_degrees: Vec<_> = graph.all_nodes().map(|n| graph.out_degree(n)).collect();
    let avg_out_degree = out_degrees.iter().sum::<usize>() as f64 / out_degrees.len() as f64;

    eprintln!("avg_in_degree: {:.2}", avg_in_degree);
    eprintln!("avg_out_degree: {:.2}", avg_out_degree);
}

fn main() -> Result<()> {
    let args = CliArgs::parse();

    match &args.command {
        Some(PoastaSubcommand::Align(v)) => align_subcommand(v)?,
        Some(PoastaSubcommand::View(v)) => view_subcommand(v)?,
        Some(PoastaSubcommand::Stats(v)) => stats_subcommand(v)?,
        None => return Err(PoastaError::Other).with_context(|| "No subcommand given.".to_string()),
    };

    Ok(())
}
