//! Graph serialization to disk using serde

use std::{char, fmt};
use std::fs::File;
use std::io::{Read, Write, BufReader};
use std::path::Path;
use petgraph::dot::Dot;
use petgraph::graph::IndexType;
use crate::errors::PoastaError;
use crate::graphs::poa::{POAGraph, POAGraphWithIx, POANodeData, POANodeIndex, Sequence};

use noodles::fasta;
use flate2::read::GzDecoder;

pub fn save_graph(graph: &POAGraphWithIx, out: impl Write) -> Result<(), PoastaError> {
    bincode::serialize_into(out, graph)?;

    Ok(())
}

pub fn load_graph(reader: impl Read) -> Result<POAGraphWithIx, PoastaError> {
    let graph: POAGraphWithIx = bincode::deserialize_from(reader)?;

    Ok(graph)
}

pub fn load_graph_from_fasta_msa(path: impl AsRef<Path>) -> Result<POAGraphWithIx, PoastaError> {
    // Let's read the sequences from the given FASTA
    let p = path.as_ref();
    let is_gzipped = p.file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
    let reader_inner: Box<dyn std::io::BufRead> = if is_gzipped {
        Box::new(File::open(p)
            .map(GzDecoder::new)
            .map(BufReader::new)?)
    } else {
        Box::new(File::open(p)
            .map(BufReader::new)?)
    };
    let mut reader = fasta::Reader::new(reader_inner);

    let mut graph = POAGraph::<usize>::new();
    let mut nodes_per_col: Vec<Vec<POANodeIndex<usize>>> = Vec::new();
    for (seq_id, record) in reader.records().enumerate() {
        let seq = record?;

        if seq.sequence().len() > nodes_per_col.len() {
            nodes_per_col.resize(seq.sequence().len(), Vec::default());
        }

        let chars: &[u8] = seq.sequence().as_ref();
        let mut prev_node = None;
        for (col, c) in chars.iter().enumerate() {
            if *c == b'-' {
                continue;
            }

            let node_ix = nodes_per_col[col].iter()
                .find(|v| graph.graph[**v].symbol == *c)
                .copied()
                .unwrap_or_else(|| {
                    // No existing node with current symbol found, create a new one
                    let new_node = graph.graph.add_node(POANodeData::new(*c));

                    for other_node in &nodes_per_col[col] {
                        graph.graph[*other_node].aligned_nodes.push(new_node);
                        graph.graph[new_node].aligned_nodes.push(*other_node);
                    }

                    nodes_per_col[col].push(new_node);

                    new_node
                });

            if let Some(prev) = prev_node {
                graph.add_edge(prev, node_ix, seq_id, 2);
            } else {
                // First node of this sequence
                graph.sequences.push(Sequence(seq.name().to_string(), node_ix));
            }

            prev_node = Some(node_ix);
        }
    }

    graph.post_process()?;

    Ok(POAGraphWithIx::USIZE(graph))
}

pub fn format_as_dot<Ix: IndexType>(writer: &mut impl fmt::Write, graph: &POAGraph<Ix>) -> fmt::Result {
    let transformed = graph.graph.map(
        |ix, data|
            format!("{:?} ({:?})", char::from(data.symbol), ix.index()),
        |_, data|
            format!("{}, {:?}", data.weight, data.sequence_ids)
    );

    let dot = Dot::new(&transformed);

    writeln!(writer, "{}", dot)?;

    Ok(())
}