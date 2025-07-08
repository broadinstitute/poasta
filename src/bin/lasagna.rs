use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use noodles::{fasta, fastq};
use flate2::read::MultiGzDecoder;
use poasta::bubbles::index::BubbleIndex;
use poasta::io::gfa::FieldValue;
use rustc_hash::FxHashMap;

use poasta::aligner::config::{AffineMinGapCost, AlignmentConfig};
use poasta::aligner::scoring::{AlignmentType, GapAffine, Score};
use poasta::aligner::PoastaAligner;
use poasta::errors::PoastaError;
use poasta::graphs::poa::{POAGraph, POANodeIndex};
use poasta::graphs::AlignableRefGraph;
use poasta::io::gfa::Field;
use poasta::io::gaf::{alignment_to_gaf, GAFRecord};
use poasta::io::graph::{load_graph_from_gfa, GraphSegments, POAGraphFromGFA};

/// The various output formats supported by Lasagna
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum OutputType {
    /// Output sequence-to-graph alignments as GAF
    Gaf,
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
struct LasagnaCli {
    /// Set verbosity level. Use multiple times to increase the verbosity level.
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    #[command(subcommand)]
    command: Option<LasagnaSubcommand>,
}

#[derive(Subcommand, Debug)]
enum LasagnaSubcommand {
    /// Perform multiple sequence alignment and create or update POA graphs
    Align(AlignArgs),
}

#[derive(Args, Debug)]
struct AlignArgs {
    /// Input graph to align sequences to. Should be acylic.
    #[clap(help_heading = "Inputs")]
    graph: PathBuf,

    /// Sequences to align in FASTA or FASTQ format.
    #[clap(help_heading = "Inputs")]
    sequences: PathBuf,
    
    #[arg(short = 'j', long, default_value = "1")]
    #[clap(help_heading = "Processing")]
    num_threads: Option<usize>,

    /// Output filename. If not given, defaults to stdout
    #[arg(short, long)]
    #[clap(help_heading = "Outputs")]
    output: Option<PathBuf>,

    /// Output file type.
    #[arg(value_enum, short = 'O', long)]
    #[clap(help_heading = "Outputs")]
    output_type: Option<OutputType>,

    /// Alignment span, either global alignment, semi-global alignment or ends-free alignment.
    #[arg(short = 'm', long, default_value = "global")]
    #[clap(help_heading = "Alignment configuration")]
    alignment_span: AlignmentSpan,

    /// Penalty for mismatching bases
    #[arg(short = 'n', default_value = "4")]
    #[clap(help_heading = "Alignment configuration")]
    cost_mismatch: Option<u8>,

    /// Penalty for opening a new gap
    #[arg(short = 'g', default_value = "6")]
    #[clap(help_heading = "Alignment configuration")]
    cost_gap_open: Option<u8>,

    /// Penalty for extending a gap
    #[arg(short = 'e', default_value = "2")]
    #[clap(help_heading = "Alignment configuration")]
    cost_gap_extend: Option<u8>,
}


struct SequenceRecord(String, Vec<u8>);


fn align_sequence<Ix, C>(
    graph: &POAGraph<Ix>, 
    graph_segments: &GraphSegments<Ix>,
    node_to_segment: &FxHashMap<POANodeIndex<Ix>, (usize, usize)>,
    aligner: &mut PoastaAligner<C>,
    bubble_index: Arc<BubbleIndex<POANodeIndex<Ix>>>,
    seq_name: &str, 
    sequence: &[u8],
) -> Option<GAFRecord>
where 
    Ix: petgraph::graph::IndexType + serde::de::DeserializeOwned,
    C: AlignmentConfig,
{
    let result = aligner.align_with_existing_bubbles::<u32, _>(graph, sequence, bubble_index);
    
    alignment_to_gaf(graph, graph_segments, seq_name, sequence, &result.alignment, node_to_segment)
        .map(|mut r| {
            if let Score::Score(v) = result.score {
                r.additional_fields.push(Field {
                    tag: "AS".to_string(),
                    value: FieldValue::Integer(v.get() as i64),
                });
            }
            
            r
        })
}


fn read_fastq(
    tx: crossbeam_channel::Sender<SequenceRecord>,
    reader_inner: impl BufRead + Send,
) -> Result<()> {
    let mut reader = fastq::io::Reader::new(reader_inner);
    
    for record in reader.records() {
        match record {
            Ok(r) => {
                let seq_name = std::str::from_utf8(r.name()).unwrap();
                tx.send(SequenceRecord(seq_name.to_string(), r.sequence().to_vec()))?;
            },
            Err(err) => {
                eprintln!("Could not parse record: {err}");
            }
        }
    }
    
    Ok(())
}

fn read_fasta(
    tx: crossbeam_channel::Sender<SequenceRecord>,
    reader_inner: impl BufRead + Send,
) -> Result<()> {
    let mut reader = fasta::io::Reader::new(reader_inner);
    
    for record in reader.records() {
        match record {
            Ok(r) => {
                let seq_name = std::str::from_utf8(r.name()).unwrap();
                tx.send(SequenceRecord(seq_name.to_string(), r.sequence().as_ref().to_vec()))?;
            },
            Err(err) => {
                eprintln!("Could not parse record: {err}");
            }
        }
    }
    
    Ok(())
}


fn align_subcommand(args: &AlignArgs) -> Result<()> {
    let POAGraphFromGFA { graph, graph_segments } = load_graph_from_gfa::<u32>(&args.graph)
        .with_context(|| "Could not load graph from GFA.")?;
     
    // Construct bubble index
    let bubble_index = Arc::new(BubbleIndex::new(&graph));
    
    // TODO: this is maybe a bit memory inefficient, storing segment id for every node
    let mut node_to_segment = FxHashMap::default();
    for (segment_ix, _) in graph_segments.names.iter().enumerate() {
        let start_node = graph_segments.start_nodes[segment_ix];
        let end_node = graph_segments.end_nodes[segment_ix];
        
        node_to_segment.insert(start_node, (segment_ix, 0));
        
        if start_node == end_node {
            continue;
        }
        
        let mut curr_node = start_node;
        let mut segment_pos = 1;
        while let Some(succ) = graph.successors(curr_node).next() {
            node_to_segment.insert(succ, (segment_ix, segment_pos));
            curr_node = succ;
            segment_pos += 1;
            
            if curr_node == end_node {
                break;
            }
        }
    }
    
    let scoring = GapAffine::new(
        args.cost_mismatch.unwrap_or(4),
        args.cost_gap_extend.unwrap_or(2),
        args.cost_gap_open.unwrap_or(6),
    );
    
    // Read sequences to align
    let is_gzipped = args.sequences.extension().map_or(false, |ext| ext == "gz");
    let non_gzip_fname = if is_gzipped { args.sequences.with_extension("") } else { args.sequences.clone() };
    let is_fastq = non_gzip_fname.extension().map_or(false, |ext| ext == "fastq" || ext == "fq");
    
    let reader_inner: Box<dyn std::io::BufRead + Send> = if is_gzipped {
        Box::new(
            File::open(&args.sequences)
                .map(MultiGzDecoder::new)
                .map(BufReader::new)?,
        )
    } else {
        Box::new(File::open(&args.sequences).map(BufReader::new)?)
    };
    
    let (tx, rx) = crossbeam_channel::unbounded();
    let (tx_out, rx_out) = crossbeam_channel::unbounded();
    
    thread::scope(|scope| -> anyhow::Result<()> {
        // Spawn thread that reads the file
        if is_fastq {
            scope.spawn(move || {
                read_fastq(tx, reader_inner)
            });
        } else {
            scope.spawn(move || {
                read_fasta(tx, reader_inner)
            });
        }
        
        // Spawn thread that writes the output
        scope.spawn(move || -> anyhow::Result<()> {
            if let Some(output) = &args.output {
                let mut writer = BufWriter::new(File::create(output)?);
                
                while let Ok(record) = rx_out.recv() {
                    writeln!(writer, "{}", record)?;
                }
            } else {
                while let Ok(record) = rx_out.recv() {
                    println!("{}", record);
                }
            }
            
            Ok(())
        });
        
        // Spawn aligner threads
        for _ in 0..args.num_threads.unwrap_or(1) {
            let graph = &graph;
            let graph_segments = &graph_segments;
            let node_to_segment = &node_to_segment;
            let thread_rx = rx.clone();
            let tx_out_thread = tx_out.clone();
            let bubble_index_thread = bubble_index.clone();
            
            scope.spawn(move || -> Result<()> {
                let mut aligner = PoastaAligner::new(AffineMinGapCost(scoring), AlignmentType::Global);
                
                while let Ok(SequenceRecord(seq_name, sequence)) = thread_rx.recv() {
                    let result = align_sequence(graph, graph_segments, node_to_segment, &mut aligner, bubble_index_thread.clone(), &seq_name, &sequence);
                    
                    if let Some(r) = result {
                        tx_out_thread.send(r)?;
                    }
                }
                
                Ok(())
            });
        }
        
        drop(tx_out);
        
        Ok(())
    })?;
    
    Ok(())
}


fn main() -> Result<()> {
    let args = LasagnaCli::parse();

    match &args.command {
        Some(LasagnaSubcommand::Align(v)) => align_subcommand(v)?,
        None => return Err(PoastaError::Other).with_context(|| "No subcommand given.".to_string()),
    };

    Ok(())
}