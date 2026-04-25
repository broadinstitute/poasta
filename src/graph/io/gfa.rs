//! Write POA graphs as GFA 1.0.

use std::io::Write;

use petgraph::Incoming;
use petgraph::graph::IndexType;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};

use crate::errors::PoastaIOError;
use crate::graph::consensus::heaviest_bundle_consensus;
use crate::graph::poa::POAGraph;
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

/// Options controlling GFA output generation.
#[derive(Default, Copy, Clone, Debug)]
pub struct GfaOutputOptions {
    /// When true, append an extra `P` line named "consensus" that traces the graph's
    /// heaviest-bundling consensus path.
    pub include_consensus: bool,
}

const CONSENSUS_NAME: &str = "consensus";

/// Serialise a POA graph to GFA 1.0.
///
/// Emits:
/// - `H\tVN:Z:1.0` header
/// - one `S` segment line per real (non-sentinel) node, with `RC:i:<in-edge weight sum>`
/// - one `L` link line per graph edge that does not touch the start/end sentinels
/// - one `P` path line per input sequence folded into the graph
pub fn poa_graph_to_gfa<Ix, W>(
    graph: &POAGraph<Ix>,
    mut out: W,
    opts: GfaOutputOptions,
) -> Result<(), PoastaIOError>
where
    Ix: IndexType,
    W: Write,
{
    let start = graph.start_node();
    let end = graph.end_node();

    writeln!(out, "H\tVN:Z:1.0").map_err(|source| PoastaIOError::FileWriteError { source })?;

    let mut n_segments = 0usize;
    for node in graph.all_nodes_iter() {
        if node == start || node == end {
            continue;
        }
        let symbol = graph.graph[node].symbol as char;
        let in_weight: usize = graph
            .graph
            .edges_directed(node, Incoming)
            .map(|e| e.weight().weight)
            .sum();

        writeln!(out, "S\tn{}\t{}\tRC:i:{}", node.index(), symbol, in_weight)
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
        n_segments += 1;
    }

    let mut n_links = 0usize;
    for edge in graph.graph.edge_references() {
        let src = edge.source();
        let tgt = edge.target();
        if src == start || src == end || tgt == start || tgt == end {
            continue;
        }
        writeln!(out, "L\tn{}\t+\tn{}\t+\t0M", src.index(), tgt.index())
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
        n_links += 1;
    }

    for (seq_idx, seq) in graph.sequences.iter().enumerate() {
        let mut segments: Vec<String> = Vec::new();
        let mut curr = seq.start_node();
        loop {
            segments.push(format!("n{}+", curr.index()));
            let next = graph
                .graph
                .edges(curr)
                .find(|e| e.weight().sequence_ids.contains(&seq_idx))
                .map(|e| e.target());
            match next {
                Some(n) if n != end => curr = n,
                _ => break,
            }
        }

        writeln!(out, "P\t{}\t{}\t*", seq.name(), segments.join(","))
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
    }

    let mut consensus_emitted = false;
    if opts.include_consensus {
        let consensus = heaviest_bundle_consensus(graph);
        if !consensus.is_empty() {
            let segments: Vec<String> = consensus
                .iter()
                .map(|n| format!("n{}+", n.index()))
                .collect();
            writeln!(out, "P\t{}\t{}\t*", CONSENSUS_NAME, segments.join(","))
                .map_err(|source| PoastaIOError::FileWriteError { source })?;
            consensus_emitted = true;
        }
    }

    tracing::debug!(
        segments = n_segments,
        links = n_links,
        paths = graph.sequences.len(),
        consensus = consensus_emitted,
        "wrote GFA"
    );

    Ok(())
}
