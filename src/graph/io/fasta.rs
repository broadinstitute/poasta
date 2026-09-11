//! Read and write POA graphs as FASTA multiple sequence alignments.

use std::io::{BufRead, Write};

use noodles::fasta::record::{Definition, Sequence as NoodlesSequence};
use noodles::fasta::{self, Record};
use petgraph::graph::IndexType;
use petgraph::visit::EdgeRef;

const NO_COL: u32 = u32::MAX;

use crate::errors::PoastaIOError;
use crate::graph::consensus::heaviest_bundle_consensus;
use crate::graph::io::CONSENSUS_NAME;
use crate::graph::poa::{POAGraph, POANodeIndex, Sequence};
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

use super::super::poa::graph_impl::POANodeData;

/// Options controlling FASTA output generation.
#[derive(Default, Copy, Clone, Debug)]
pub struct FastaOutputOptions {
    /// When true, append an extra FASTA record named "consensus" (with gap padding) as the last
    /// row of the MSA. Ignored if `consensus_only` is also set.
    pub include_consensus: bool,
    /// When true, write only the consensus as one FASTA record without MSA columns or gaps.
    /// Takes precedence over `include_consensus`.
    pub consensus_only: bool,
}

/// Load a POA graph from a FASTA multiple sequence alignment.
///
/// Each record is interpreted as one row of the MSA; rows must be equal length (shorter rows are
/// tolerated by only reading up to their actual length). The character `-` denotes a gap.
///
/// Columns that share the same character across rows collapse to one node with the
/// appropriate `aligned_nodes` cross-links, mirroring the construction used when aligning
/// sequences incrementally.
pub fn load_graph_from_fasta_msa<Ix, R>(reader: R) -> Result<POAGraph<Ix>, PoastaIOError>
where
    Ix: IndexType,
    R: BufRead,
{
    let mut fasta_reader = fasta::io::Reader::new(reader);
    let mut graph = POAGraph::<Ix>::new();
    let mut nodes_per_col: Vec<Vec<POANodeIndex<Ix>>> = Vec::new();

    for (seq_id, record_result) in fasta_reader.records().enumerate() {
        let record = record_result.map_err(|source| PoastaIOError::FileReadError { source })?;
        let seq_name = std::str::from_utf8(record.name())?.to_owned();
        let chars: &[u8] = record.sequence().as_ref();

        if chars.len() > nodes_per_col.len() {
            nodes_per_col.resize(chars.len(), Vec::default());
        }

        let mut prev_node: Option<POANodeIndex<Ix>> = None;
        let mut first_node: Option<POANodeIndex<Ix>> = None;

        for (col, &c) in chars.iter().enumerate() {
            if c == b'-' {
                continue;
            }

            let node_ix = match nodes_per_col[col]
                .iter()
                .find(|v| graph.graph[**v].symbol == c)
                .copied()
            {
                Some(existing) => existing,
                None => {
                    let new_node = graph.graph.add_node(POANodeData::new(c));
                    for other_node in &nodes_per_col[col] {
                        graph.graph[*other_node].aligned_nodes.push(new_node);
                        graph.graph[new_node].aligned_nodes.push(*other_node);
                    }
                    nodes_per_col[col].push(new_node);
                    new_node
                }
            };

            if let Some(p) = prev_node {
                graph.add_edge(p, node_ix, seq_id, 2);
            } else {
                first_node = Some(node_ix);
            }
            prev_node = Some(node_ix);
        }

        if let Some(first) = first_node {
            graph.sequences.push(Sequence(seq_name, first));
        }
    }

    graph
        .post_process()
        .map_err(|_| PoastaIOError::InvalidFormat)?;

    tracing::debug!(
        sequences = graph.sequences.len(),
        columns = nodes_per_col.len(),
        nodes = graph.node_count(),
        "loaded POA graph from FASTA MSA"
    );

    Ok(graph)
}

/// Write a POA graph as a FASTA multiple sequence alignment.
///
/// Each sequence that was folded into the graph is emitted as one record with gaps (`-`) in
/// columns it does not cover. Columns correspond to clusters of mutually-aligned nodes
/// (the `aligned_nodes` cross-links created during alignment).
///
/// When `opts.consensus_only` is set, emits only one FASTA record with the consensus
/// sequence (no gaps). When `opts.include_consensus` is set, appends the consensus as the last
/// record of the MSA with gap padding matching the other rows.
pub fn poa_graph_to_fasta<Ix, W>(
    graph: &POAGraph<Ix>,
    writer: W,
    opts: FastaOutputOptions,
) -> Result<(), PoastaIOError>
where
    Ix: IndexType,
    W: Write,
{
    let mut fasta_writer = fasta::io::Writer::new(writer);

    if opts.consensus_only {
        let consensus = heaviest_bundle_consensus(graph);
        let seq: Vec<u8> = consensus.iter().map(|&n| graph.graph[n].symbol).collect();
        tracing::debug!(length = seq.len(), "writing consensus-only FASTA");
        let record = Record::new(
            Definition::new(CONSENSUS_NAME.to_owned(), None),
            NoodlesSequence::from(seq),
        );
        fasta_writer
            .write_record(&record)
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
        return Ok(());
    }

    let start = graph.start_node();
    let end = graph.end_node();

    // Assign every real node to an MSA column. Nodes in the same `aligned_nodes` cluster share a
    // column. We iterate in topological order so the column counter mirrors the sequence of
    // the underlying DAG.
    //
    // Indexed by `node.index()` rather than a hash map: petgraph node indices are dense and
    // bounded, so a Vec lookup is one load instead of a hash + probe.
    let topo = graph.topological_order();
    let max_node_idx = graph
        .all_nodes_iter()
        .map(|n| n.index())
        .max()
        .unwrap_or(0);
    let mut col_of: Vec<u32> = vec![NO_COL; max_node_idx + 1];
    let mut n_cols = 0u32;

    for node in topo {
        if *node == start || *node == end {
            continue;
        }
        if col_of[node.index()] != NO_COL {
            continue;
        }

        let col = n_cols;
        n_cols += 1;
        col_of[node.index()] = col;
        for other in &graph.graph[*node].aligned_nodes {
            col_of[other.index()] = col;
        }
    }
    let n_cols = n_cols as usize;

    tracing::debug!(
        columns = n_cols,
        sequences = graph.sequences.len(),
        "writing FASTA MSA"
    );

    // Walk each sequence one at a time, picking the next edge by binary-searching
    // `seq_idx` in the edge's `sequence_ids`. `POAGraph::add_edge` appends sequence ids in
    // ingest order, so `sequence_ids` is sorted ascending, so binary search is sound.
    // Peak extra memory is one row buffer instead of the full `O(total_path_length)`
    // path prebuild we used to do.
    for (seq_idx, seq) in graph.sequences.iter().enumerate() {
        let mut row = vec![b'-'; n_cols];
        let mut cur = seq.start_node();
        while cur != end {
            let col = col_of[cur.index()];
            if col == NO_COL {
                return Err(PoastaIOError::InvalidFormat);
            }
            row[col as usize] = graph.graph[cur].symbol;

            let next = graph
                .graph
                .edges(cur)
                .find(|e| e.weight().sequence_ids.binary_search(&seq_idx).is_ok())
                .map(|e| e.target());
            match next {
                Some(v) => cur = v,
                None => break,
            }
        }

        let definition = Definition::new(seq.name().clone(), None);
        let sequence = NoodlesSequence::from(row);
        let record = Record::new(definition, sequence);
        fasta_writer
            .write_record(&record)
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
    }

    if opts.include_consensus {
        let consensus = heaviest_bundle_consensus(graph);
        let mut row = vec![b'-'; n_cols];
        for node in &consensus {
            let col = col_of[node.index()];
            if col == NO_COL {
                return Err(PoastaIOError::InvalidFormat);
            }
            row[col as usize] = graph.graph[*node].symbol;
        }
        let record = Record::new(
            Definition::new(CONSENSUS_NAME.to_owned(), None),
            NoodlesSequence::from(row),
        );
        fasta_writer
            .write_record(&record)
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
    }

    Ok(())
}
