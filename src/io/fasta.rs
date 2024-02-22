use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::iter;
use petgraph::graph::IndexType;
use petgraph::visit::EdgeRef;
use serde::de::DeserializeOwned;
use noodles::fasta::{self as fasta, record::{Definition, Sequence}, Record};
use rustc_hash::FxHashSet;

use crate::errors::PoastaError;
use crate::graphs::AlignableRefGraph;
use crate::graphs::poa::{POAGraph, POANodeIndex};


fn fasta_aln_for_seq<Ix>(
    graph: &POAGraph<Ix>,
    node_to_column: &HashMap<POANodeIndex<Ix>, usize>,
    seq_id: usize,
    start_node: POANodeIndex<Ix>
) -> Vec<u8>
    where Ix: IndexType + DeserializeOwned
{
    let mut seq = Vec::new();
    let mut curr = Some(start_node);
    let mut last_col = 0;
    while let Some(n) = curr {
        let node_col = node_to_column[&n];

        let gap_length = node_col.saturating_sub(1) - last_col;
        seq.extend(iter::repeat(b'-').take(gap_length));
        seq.push(graph.get_symbol(n));

        curr = None;
        for out_edge in graph.graph.edges(n) {
            if out_edge.weight().sequence_ids.binary_search(&seq_id).is_ok() {
                curr = Some(out_edge.target())
            }
        }

        last_col = node_col;
    }

    if let Some(max_col) = node_to_column.values().max() {
        let gap_length = *max_col - last_col;
        seq.extend(iter::repeat(b'-').take(gap_length));
    }

    seq
}


pub fn poa_graph_to_fasta<Ix, W>(graph: &POAGraph<Ix>, output: W) -> Result<(), PoastaError>
    where Ix: IndexType + DeserializeOwned,
          W: Write
{
    // First assign a column in the final MSA output for every node in the graph,
    // by visiting nodes in topological order, taking into account 'aligned_nodes'
    let mut node_to_column = HashMap::new();
    let mut stack = vec![
        (graph.start_node(), RefCell::new(Vec::from_iter(graph.successors(graph.start_node()))))
    ];

    let next_valid_child = |
        succ_iter: &RefCell<Vec<POANodeIndex<Ix>>>,
        visited: &FxHashSet<POANodeIndex<Ix>>
    | -> Option<POANodeIndex<Ix>> {
        while let Some(succ) = succ_iter.borrow_mut().pop() {
            if !visited.contains(&succ) {
                return Some(succ)
            }
        }

        None
    };

    let mut visited = FxHashSet::default();
    let mut rev_postorder = vec![];
    while !stack.is_empty() {
        let (_, succ_iter) = stack.last().unwrap();

        if let Some(child) = next_valid_child(succ_iter, &visited) {
            visited.insert(child);
            let mut successors = Vec::from_iter(graph.successors(child));
            for aln_node in graph.get_aligned_nodes(child) {
                if !visited.contains(aln_node) {
                    visited.insert(*aln_node);
                    successors.extend(graph.successors(*aln_node));
                }
            }
            stack.push((child, RefCell::new(successors)));
        } else {
            let (last, _) = stack.pop().unwrap();
            rev_postorder.push(last);
        }

    }

    rev_postorder.reverse();

    let mut curr_col = 0;
    for n in &rev_postorder {
        if *n == graph.start_node() || *n == graph.end_node() {
            continue;
        }

        if !node_to_column.contains_key(n) {
            node_to_column.insert(*n, curr_col);

            for aligned_node in graph.get_aligned_nodes(*n) {
                node_to_column.insert(*aligned_node, curr_col);
            }

            curr_col += 1;
        }
    }

    let mut writer = fasta::Writer::new(output);

    for (seq_id, seq) in graph.sequences.iter().enumerate() {
        let header = Definition::new(seq.name(), None);

        let seq = Sequence::from_iter(fasta_aln_for_seq(graph, &node_to_column, seq_id, seq.start_node()));
        let record = Record::new(header, seq);

        writer.write_record(&record)?;
    }


    Ok(())
}