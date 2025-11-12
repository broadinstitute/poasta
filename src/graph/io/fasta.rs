//! Importing and exporting POA graphs from and to FASTA files.

use std::collections::HashMap;
use std::io::Write;
use std::ops::Range;

use petgraph::graph::IndexType;
use rustc_hash::FxHashMap;

use crate::align::traits::{AlignableGraph, AlignableGraphNodePos};
use crate::errors::PoastaError;
use crate::graph::alignment::POANodePos;
use crate::graph::poa::{AlignedInterval, POANodeIndex, POASeqGraph};
use crate::graph::traits::{GraphBase, GraphWithNodeLengths};
use crate::graph::utils::rev_postorder_nodes;

/// During import of a POA graph from a FASTA file, this struct keeps track
/// which nodes cover which columns of the MSA.
#[derive(Debug, Default)]
pub struct MSANodeCover<Ix>
where
    Ix: IndexType,
{
    col_to_nodes: Vec<Vec<POANodePos<Ix>>>,
    node_to_cols: HashMap<POANodeIndex<Ix>, Range<usize>>,
}

impl<Ix> MSANodeCover<Ix>
where
    Ix: IndexType,
{
    pub fn new() -> Self {
        MSANodeCover {
            col_to_nodes: Vec::default(),
            node_to_cols: HashMap::default(),
        }
    }

    pub fn msa_length(&self) -> usize {
        self.col_to_nodes.len()
    }

    pub fn resize(&mut self, msa_length: usize) {
        self.col_to_nodes.resize_with(msa_length, Vec::default);
    }

    pub fn add_node(&mut self, node: POANodeIndex<Ix>, cols: Range<usize>) {
        self.node_to_cols.insert(node, cols.clone());
        for (i, col) in cols.enumerate() {
            self.col_to_nodes[col].push(POANodePos::new(node, i));
        }
    }

    pub fn split_node(
        &mut self,
        node: POANodeIndex<Ix>,
        split_pos: usize,
        left_node: POANodeIndex<Ix>,
        right_node: POANodeIndex<Ix>,
    ) {
        let orig_range = self.node_to_cols.remove(&node).unwrap();

        self.node_to_cols
            .insert(left_node, orig_range.start..(orig_range.start + split_pos));
        self.node_to_cols
            .insert(right_node, (orig_range.start + split_pos)..orig_range.end);

        for col in orig_range.clone() {
            self.col_to_nodes[col]
                .iter_mut()
                .filter(|node_pos| node_pos.node() == node)
                .for_each(|node_pos| {
                    if node_pos.pos() < split_pos {
                        *node_pos = POANodePos::new(left_node, node_pos.pos());
                    } else {
                        *node_pos = POANodePos::new(right_node, node_pos.pos() - split_pos);
                    }
                });
        }
    }

    pub fn has_match<G>(&self, graph: &G, col: usize, symbol: u8) -> Option<POANodePos<Ix>>
    where
        G: AlignableGraph<NodePosType = POANodePos<Ix>>,
    {
        self.col_to_nodes[col]
            .iter()
            .find(|node_pos| graph.get_node_symbol(**node_pos) == symbol)
            .copied()
    }

    pub fn has_nodes_covering_col(&self, col: usize) -> bool {
        !self.col_to_nodes[col].is_empty()
    }

    pub fn get_nodes_for_col(&self, col: usize) -> &[POANodePos<Ix>] {
        &self.col_to_nodes[col]
    }
}


/// Divide each node into one or more intervals
///
/// Intervals are defined by a node's associated 'aligned intervals', i.e., ranges on of sequences aligned
/// to another sequence range on another node.
///
/// The 'aligned intervals', however, don't necessarily cover an entire node, and might include redundant
/// intervals, e.g., when a range is aligned to multiple different nodes.
///
/// This struct processes the list of aligned intervals of each node and 1) removes redundant intervals, i.e.,
/// those within another larger interval, and 2) adds additional intervals covering node ranges not covered
/// by the aligned intervals.
struct IntervalsOnNodes<Ix>
where
    Ix: IndexType,
{
    node_intervals: FxHashMap<POANodeIndex<Ix>, Vec<(usize, usize)>>,
}

impl<Ix: IndexType> IntervalsOnNodes<Ix> {
    fn new(graph: &POASeqGraph<Ix>) -> Self {
        let mut node_intervals = FxHashMap::default();

        for n in graph.all_nodes_iter() {
            let mut aln_ivals = graph.get_node_aligned_intervals(n).to_owned();
            // Sort by node start, longer intervals first.
            aln_ivals.sort_by_key(|e| (e.node_start(), -(e.length() as isize)));

            // Ignore intervals that completely fall within another
            // Easy to check because of sorting above
            let mut node_ivals = aln_ivals.into_iter()
                .fold(Vec::default(), |mut ivals, elem| {
                    if let Some((_, ie)) = ivals.last().copied() {
                        // Check for gap
                        if elem.node_start() > ie {
                            let gap = elem.node_start() - ie;
                            ivals.push((ie, ie + gap))
                        }

                        // Check if ival starts after the previous interval (i.e., not sub-interval of another)
                        if elem.node_start() >= ie {
                            ivals.push((elem.node_start(), elem.node_end()));
                        }
                    } else {
                        if elem.node_start() > 0 {
                            ivals.push((0, elem.node_start()))
                        }

                        ivals.push((elem.node_start(), elem.node_end()));
                    };

                    ivals
                });

            // Ensure interval 
            if let Some((_, ie)) = node_ivals.last().copied() {
                if ie != graph.node_length(n) {
                    node_ivals.push((ie, graph.node_length(n)));
                }
            }

            node_intervals.entry(n)
                .or_insert_with(|| node_ivals);
        }

        Self { node_intervals }
    }
}


pub fn poa_graph_to_fasta<Ix, W>(graph: &POASeqGraph<Ix>, output: W) -> Result<(), PoastaError>
where
    Ix: IndexType,
    W: Write,
{
    // For each aligned interval on a node, assign a range i..j corresponding to the columns in the FASTA MSA.
    // To determine the range, we visit nodes in topological order, and inspect each aligned interval.

    let mut col_assignments = FxHashMap::default();

    let mut curr_col = 0;
    for n in rev_postorder_nodes(graph) {
        let mut aln_ivals = graph.get_node_aligned_intervals(n).to_owned();
        // Sort by node start, longer intervals first.
        aln_ivals.sort_by_key(|e| (e.node_start(), -(e.length() as isize)));

        // First iteration we assign column i..j to each interval on this node
        let mut last_ival_end = 0;
        for aln_ival in aln_ivals {
            if aln_ival.node_end() <= last_ival_end {
                // We sorted 
                continue;
            }
            // col_assignments.entry((n, aln_ival))
            //     .or_insert_with(|| {
            //         let assigned_range = curr_col..curr_col+aln_ival.length();

            //         if 
            //     });

            // Assign same col range for intervals on aligned nodes
            let ival_on_other = aln_ival.flip_nodes(n);
            col_assignments.entry((aln_ival.other(), ival_on_other))
                .or_insert_with(|| curr_col..curr_col+aln_ival.length());

            curr_col += aln_ival.length();
        }

    }

    Ok(())
}
