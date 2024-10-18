//! Importing and exporting POA graphs from and to FASTA files.

use std::ops::Range;

use foldhash::HashMap;
use petgraph::graph::IndexType;

use crate::aligner::astar::{AlignableGraph, AlignableGraphNodePos};
use crate::graph::alignment::POANodePos;
use crate::graph::poa::POANodeIndex;

/// During import of a POA graph from a FASTA file, this struct keeps track
/// which nodes cover which columns of the MSA.
#[derive(Debug, Default)]
pub struct MSANodeCover<Ix>
where
    Ix: IndexType
{
    col_to_nodes: Vec<Vec<POANodePos<Ix>>>,
    node_to_cols: HashMap<POANodeIndex<Ix>, Range<usize>>
}


impl<Ix> MSANodeCover<Ix>
where
    Ix: IndexType
{
    pub fn new() -> Self {
        MSANodeCover {
            col_to_nodes: Vec::default(),
            node_to_cols: HashMap::default()
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
    
    pub fn split_node(&mut self, node: POANodeIndex<Ix>, split_pos: usize, left_node: POANodeIndex<Ix>, right_node: POANodeIndex<Ix>) {
        let orig_range = self.node_to_cols.remove(&node).unwrap();
        
        self.node_to_cols.insert(left_node, orig_range.start..(orig_range.start + split_pos));
        self.node_to_cols.insert(right_node, (orig_range.start + split_pos)..orig_range.end);
        
        for col in orig_range.clone() {
            self.col_to_nodes[col].iter_mut()
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
        G: AlignableGraph<NodePosType = POANodePos<Ix>>
    {
        self.col_to_nodes[col].iter()
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
