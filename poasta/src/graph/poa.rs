use std::fs::File;
use std::io::BufRead;
use std::ops::{Range, RangeBounds, Bound};

use petgraph::algo::toposort;
pub use petgraph::graph::IndexType;
use petgraph::stable_graph::Neighbors;
use petgraph::visit::EdgeRef;
use petgraph::{Incoming, Outgoing, visit::NodeIndexable};

use tracing::{debug_span, debug};
use foldhash::HashMap;
use noodles::fasta;

use crate::aligner::astar::{AlignableGraphNodeId, AlignableGraph, AlignableGraphNodePos};
use crate::aligner::utils::AlignedPair;
use crate::errors::{GraphError, PoastaIOError};
use crate::graph::alignment::AlignmentBlockType;
use crate::graph::io::dot::graph_to_dot;

use super::alignment::{AlignmentBlocks, AlignmentClassification, POANodePos};
use super::io::fasta::MSANodeCover;

pub(crate) mod graph_impl {
    use petgraph::graph::IndexType;
    use petgraph::{stable_graph::StableDiGraph, visit::GraphBase};
    
    use tracing::debug;
    
    pub type POASeqGraph<Ix> = StableDiGraph<POANodeData<Ix>, POAEdgeData, Ix>;
    pub type POANodeIndex<Ix> = <POASeqGraph<Ix> as GraphBase>::NodeId;
    
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct AlignedInterval<Ix>
    where
        Ix: IndexType,
    {
        node_start: usize,
        length: usize,
        other: POANodeIndex<Ix>,
        other_start: usize,
    }
    
    impl<Ix> AlignedInterval<Ix>
    where 
        Ix: IndexType,
    {
        pub fn new(node_start: usize, length: usize, other: POANodeIndex<Ix>, other_start: usize) -> Self {
            AlignedInterval {
                node_start,
                length,
                other,
                other_start,
            }
        }
        
        #[inline(always)]
        pub fn length(&self) -> usize {
            self.length
        }
        
        #[inline]
        pub fn increase_length(&mut self, v: usize) {
            self.length += v;
        }
        
        #[inline(always)]
        pub fn other(&self) -> POANodeIndex<Ix> {
            self.other
        }
        
        #[inline(always)]
        pub fn node_start(&self) -> usize {
            self.node_start
        }
        
        #[inline(always)]
        pub fn node_end(&self) -> usize {
            self.node_start + self.length
        }
        
        #[inline(always)]
        pub fn other_start(&self) -> usize {
            self.other_start
        }
        
        #[inline(always)]
        pub fn other_end(&self) -> usize {
            self.other_start + self.length
        }
        
        #[inline]
        pub fn with_new_other(&self, other: POANodeIndex<Ix>) -> Self {
            AlignedInterval {
                node_start: self.node_start,
                length: self.length,
                other,
                other_start: self.other_start,
            }
        }
        
        #[inline]
        pub fn pos_to_other(&self, pos: usize) -> Option<usize> {
            if pos < self.node_start || pos >= self.node_end() {
                return None;
            }
            
            Some(self.other_start + (pos - self.node_start))
        }
        
        #[inline]
        pub fn move_node_start(&self, diff: usize) -> Self {
            AlignedInterval {
                node_start: self.node_start - diff,
                length: self.length,
                other: self.other,
                other_start: self.other_start,
            }
        }
        
        #[inline]
        pub fn move_other_start(&self, other: POANodeIndex<Ix>, diff: usize) -> Self {
            AlignedInterval {
                node_start: self.node_start,
                length: self.length,
                other,
                other_start: self.other_start - diff,
            }
        }
        
        #[inline]
        pub fn split_at_left(&self, split_pos: usize) -> Self {
            assert!(split_pos > self.node_start);
            let length_diff = self.node_end().saturating_sub(split_pos);
            
            AlignedInterval {
                node_start: self.node_start,
                length: self.length() - length_diff,
                other: self.other,
                other_start: self.other_start,
            }
        }
        
        
        #[inline]
        pub fn split_at_right(&self, split_pos: usize) -> Self {
            assert!(split_pos < self.node_end());
            let length_diff = split_pos.saturating_sub(self.node_start);
            
            AlignedInterval {
                node_start: self.node_start.saturating_sub(split_pos),
                length: self.length() - length_diff,
                other: self.other,
                other_start: self.other_start + length_diff,
            }
        }
        
        #[inline]
        pub fn other_split_at_left(&self, new_node: POANodeIndex<Ix>, split_pos: usize) -> Self {
            debug!("{:?}, splitpos: {:?}", self, split_pos);
            assert!(split_pos > self.other_start);
            
            if split_pos >= self.other_end() {
                return *self;
            }
            
            AlignedInterval {
                node_start: self.node_start,
                length: split_pos - self.other_start,
                other: new_node,
                other_start: self.other_start,
            }
        }
        
        #[inline]
        pub fn other_split_at_right(&self, new_node: POANodeIndex<Ix>, split_pos: usize) -> Self {
            assert!(split_pos <= self.other_end());
            
            if split_pos <= self.other_start {
                let diff = self.other_start() - split_pos;
                AlignedInterval {
                    node_start: self.node_start(),
                    length: self.length(),
                    other: new_node,
                    other_start: diff,
                }
            } else {
                let length_diff = split_pos - self.other_start();
                AlignedInterval {
                    node_start: self.node_start + length_diff,
                    length: self.other_end() - split_pos,
                    other: new_node,
                    other_start: 0,
                }
            }
        }
        
        #[inline]
        pub fn intersect(&self, aln_ival: &Self) -> Option<Self> {
            let start = std::cmp::max(self.node_start, aln_ival.node_start);
            let end = std::cmp::min(self.node_end(), aln_ival.node_end());
            
            if start >= end {
                return None;
            }
            
            let other_start = self.pos_to_other(start).unwrap();
            
            Some(AlignedInterval {
                node_start: start,
                length: end - start,
                other: self.other,
                other_start,
            })
        }
        
        #[inline]
        pub fn transitive_intersect(&self, aln_ival: &Self) -> Option<(Self, Self)> {
            let intersection = self.intersect(aln_ival)?;
            
            let ival1 = AlignedInterval {
                node_start: intersection.other_start(),
                length: intersection.length(),
                other: aln_ival.other(),
                other_start: aln_ival.pos_to_other(intersection.node_start()).unwrap(),
            };
            
            let intersection = aln_ival.intersect(self)?;
            let ival2 = AlignedInterval {
                node_start: intersection.other_start(),
                length: intersection.length(),
                other: self.other(),
                other_start: self.pos_to_other(intersection.node_start()).unwrap(),
            };
            
            Some((ival1, ival2))
        }
        
        #[inline]
        pub fn flip_nodes(&self, node: POANodeIndex<Ix>) -> Self {
            AlignedInterval {
                node_start: self.other_start,
                length: self.length,
                other: node,
                other_start: self.node_start,
            }
        }
        
    }

    #[derive(Debug)]
    pub struct POANodeData<Ix>
    where 
        Ix: IndexType,
    {
        pub sequence: Vec<u8>,
        pub weights: Vec<usize>,
        pub aligned_intervals: Vec<AlignedInterval<Ix>>,
        pub rank: Ix,
    }

    impl<Ix> POANodeData<Ix>
    where
        Ix: IndexType,
    {
        pub fn new(sequence: Vec<u8>) -> Self {
            let seq_len = sequence.len();
            POANodeData {
                sequence,
                weights: vec![1; seq_len],
                aligned_intervals: Vec::default(),
                rank: Ix::new(0)
            }
        }

        pub fn new_with_weights(sequence: Vec<u8>, weights: Vec<usize>) -> Self {
            POANodeData {
                sequence,
                weights,
                aligned_intervals: Vec::default(),
                rank: Ix::new(0),
            }
        }
    }

    #[derive(Debug, Default)]
    pub struct POAEdgeData {
        pub sequence_ids: Vec<usize>,
    }

    impl POAEdgeData {
        pub fn new_with_seq_id(seq_id: usize) -> Self {
            POAEdgeData {
                sequence_ids: vec![seq_id],
            }
        }
        
        pub fn new_with_seq_ids(seq_ids: Vec<usize>) -> Self {
            POAEdgeData {
                sequence_ids: seq_ids,
            }
        }
    }
    
    #[cfg(test)]
    mod tests {
        use super::*;
        
        #[test]
        fn test_aligned_interval_basic() {
            let aln_ival = AlignedInterval::<usize>::new(
                10,
                5,
                POANodeIndex::<usize>::new(0usize),
                20,
            );
            
            assert_eq!(aln_ival.node_start(), 10);
            assert_eq!(aln_ival.node_end(), 15);
            assert_eq!(aln_ival.other_start(), 20);
            assert_eq!(aln_ival.other_end(), 25);
            
            assert_eq!(aln_ival.pos_to_other(10), Some(20));
            assert_eq!(aln_ival.pos_to_other(14), Some(24));
            assert_eq!(aln_ival.pos_to_other(15), None);
            assert_eq!(aln_ival.pos_to_other(9), None);
            
            let new_aln_ival = aln_ival.split_at_left(13);
            assert_eq!(new_aln_ival.node_start(), 10);
            assert_eq!(new_aln_ival.node_end(), 13);
            assert_eq!(new_aln_ival.other_start(), 20);
            assert_eq!(new_aln_ival.other_end(), 23);
            
            let new_aln_ival = aln_ival.split_at_right(12);
            assert_eq!(new_aln_ival.node_start(), 0);
            assert_eq!(new_aln_ival.node_end(), 3);
            assert_eq!(new_aln_ival.other_start(), 22);
            assert_eq!(new_aln_ival.other_end(), 25);
        }
        
        #[test]
        fn test_aligned_interval_other_split() {
            let aln_ival = AlignedInterval::<usize>::new(
                10,
                5,
                POANodeIndex::<usize>::new(1), 
                20
            );
            
            let new_aln_ival = aln_ival.other_split_at_left(POANodeIndex::<usize>::new(1), 22);
            assert_eq!(new_aln_ival.node_start(), 10);
            assert_eq!(new_aln_ival.node_end(), 12);
            assert_eq!(new_aln_ival.other_start(), 20);
            assert_eq!(new_aln_ival.other_end(), 22);
            
            let new_aln_ival = aln_ival.other_split_at_left(POANodeIndex::<usize>::new(1), 26);
            assert_eq!(new_aln_ival.node_start(), 10);
            assert_eq!(new_aln_ival.node_end(), 15);
            assert_eq!(new_aln_ival.other_start(), 20);
            assert_eq!(new_aln_ival.other_end(), 25);
            
            let new_aln_ival = aln_ival.other_split_at_right(POANodeIndex::<usize>::new(1), 24);
            assert_eq!(new_aln_ival.node_start(), 14);
            assert_eq!(new_aln_ival.node_end(), 15);
            assert_eq!(new_aln_ival.other_start(), 0);
            assert_eq!(new_aln_ival.other_end(), 1);
            
            let new_aln_ival = aln_ival.other_split_at_right(POANodeIndex::<usize>::new(1), 18);
            assert_eq!(new_aln_ival.node_start(), 10);
            assert_eq!(new_aln_ival.node_end(), 15);
            assert_eq!(new_aln_ival.other_start(), 2);
            assert_eq!(new_aln_ival.other_end(), 7);
            
        }
        
        #[test]
        fn test_intersect() {
            let aln_ival = AlignedInterval::<usize>::new(
                10,
                5,
                POANodeIndex::<usize>::new(1), 
                20
            );
            
            let other_aln_ival = AlignedInterval::<usize>::new(
                12,
                6,
                POANodeIndex::<usize>::new(1), 
                22
            );
            
            let Some(intersect) = aln_ival.intersect(&other_aln_ival) else {
                panic!("No intersection found")
            };
            
            assert_eq!(intersect.node_start(), 12);
            assert_eq!(intersect.node_end(), 15);
            assert_eq!(intersect.other_start(), 22);
            assert_eq!(intersect.other_end(), 25);
        }
    }
}

use graph_impl::{AlignedInterval, POANodeData, POAEdgeData};
pub use graph_impl::POANodeIndex;


/// The partial order alignment graph with sequence-labeled nodes
///
/// In contrast to most other POA tools, we use a graph with sequence-labeled nodes
/// as opposed to character-labeled nodes. This is reduces memory usage
/// and increased efficiency of most graph algorithms since there are fewer nodes.
///
/// The downsize is increased complexity to track which nucleotides are aligned to each other.
pub struct POASeqGraph<Ix>
where
    Ix: IndexType,
{
    /// Petgraph graph object holding nodes and edges
    pub(crate) graph: graph_impl::POASeqGraph<Ix>,
    
    /// List of sequences added to the graph
    sequences: Vec<Sequence<Ix>>,
    
    /// POA graph special start node
    pub(crate) start_node: POANodeIndex<Ix>,
    
    /// POA graph special end node
    pub(crate) end_node: POANodeIndex<Ix>,
    
    /// Toplogical ordering of nodes
    toposorted: Vec<POANodeIndex<Ix>>,
}

impl<Ix> POASeqGraph<Ix>
where
    Ix: IndexType,
{
    pub fn new() -> Self {
        let mut graph = graph_impl::POASeqGraph::<Ix>::default();
        let start_node = graph.add_node(graph_impl::POANodeData::new(b"#".to_vec()));
        let end_node = graph.add_node(graph_impl::POANodeData::new(b"$".to_vec()));
        graph.add_edge(start_node, end_node, POAEdgeData::default());

        POASeqGraph {
            graph,
            sequences: Vec::new(),
            start_node,
            end_node,
            toposorted: vec![start_node, end_node],
        }
    }
    
    pub fn try_from_fasta_msa(ifile: &mut impl BufRead) -> Result<Self, PoastaIOError> {
        let mut graph = Self::new();
        
        let mut msa_nodes_cover = MSANodeCover::new();
        let mut reader = fasta::Reader::new(ifile);
        
        let mut last_len = None;
        for record in reader.records() {
            let record = record?;
            let seq = record.sequence().as_ref();
            
            if let Some(len) = last_len {
                if len != seq.len() {
                    return Err(PoastaIOError::InvalidFormat);
                }
            }            
            
            if seq.len() > msa_nodes_cover.msa_length() {
                msa_nodes_cover.resize(seq.len());
            }
            
            let seq_name = std::str::from_utf8(record.name())?;
            graph.add_msa_sequence(seq_name, seq, &mut msa_nodes_cover)?;
            
            last_len = Some(seq.len());
        }
        
        Ok(graph)
    }

    pub fn is_empty(&self) -> bool {
        // An empty graph still contains the special start and end node
        self.graph.node_count() == 2
    }
    
    pub fn num_nodes(&self) -> usize {
        self.graph.node_count() - 2
    }
    
    pub fn num_nodes_with_start_end(&self) -> usize {
        self.graph.node_count()
    }
    
    pub fn node_data(&self, node: POANodeIndex<Ix>) -> &POANodeData<Ix> {
        self.graph.node_weight(node).unwrap()
    }
    
    pub fn node_data_mut(&mut self, node: POANodeIndex<Ix>) -> &mut POANodeData<Ix> {
        self.graph.node_weight_mut(node).unwrap()
    }
    
    pub fn get_nodes(&self) -> impl Iterator<Item = POANodeIndex<Ix>> + '_ {
        self.graph.node_indices()
    }
    
    pub fn get_sequences(&self) -> &[Sequence<Ix>] {
        &self.sequences
    }
    
    /// Get the path through the graph for a particular sequence
    pub fn path_for_sequence(&self, seq_id: usize) -> Vec<POANodeIndex<Ix>> {
        // TODO: make this into an iterator
        let start_node = self.sequences[seq_id].start_node();
        let mut path = vec![start_node];
        
        loop {
            let last_node = path.last().unwrap();
            let next_node = self.graph.neighbors(*last_node)
                .find(|v| {
                    let edge = self.graph.find_edge(*last_node, *v).unwrap();
                    
                    self.graph.edge_weight(edge).unwrap()
                        .sequence_ids
                        .binary_search(&seq_id)
                        .is_ok()
                });
            
            let Some(next_node) = next_node else {
                break;
            };
            
            if next_node == self.end_node {
                break;
            }
            
            path.push(next_node);
        }
        
        path
    }
    
    /// Concatenate the sequence of a path through the graph
    pub fn path_sequence(&self, path: impl IntoIterator<Item=POANodeIndex<Ix>>) -> Result<Vec<u8>, GraphError<Ix>> {
        let mut seq = Vec::new();
        
        let mut last = None;
        for node in path.into_iter() {
            if let Some(p) = last {
                if self.graph.find_edge(p, node).is_none() {
                    return Err(GraphError::InvalidEdge(p, node));
                }
            }
            
            let node_data = self.node_data(node);
            seq.extend(node_data.sequence.iter());
            
            last = Some(node);
        }
        
        Ok(seq)
    }
    
    pub fn add_node(&mut self, sequence: Vec<u8>) -> graph_impl::POANodeIndex<Ix> {
        let node_data = graph_impl::POANodeData::new(sequence);
        self.add_node_with_data(node_data)
    }

    pub fn add_node_with_weights(
        &mut self,
        sequence: Vec<u8>,
        weights: Vec<usize>,
    ) -> POANodeIndex<Ix> {
        let node_data = POANodeData::new_with_weights(sequence, weights);
        self.add_node_with_data(node_data)
    }
    
    fn add_node_with_data(
        &mut self,
        node_data: POANodeData<Ix>,
    ) -> POANodeIndex<Ix> {
        let new_ix = self.graph.add_node(node_data);
        
        if self.is_empty() {
            self.graph.remove_edge(self.graph.find_edge(self.start_node, self.end_node).unwrap());
            self.graph.add_edge(self.start_node, new_ix, POAEdgeData::new_with_seq_id(self.sequences.len()));
            self.graph.add_edge(new_ix, self.end_node, POAEdgeData::new_with_seq_id(self.sequences.len()));
        }
        
        new_ix
    }
    
    pub fn add_aligned_sequence(
        &mut self,
        sequence_name: &str,
        sequence: impl AsRef<[u8]>,
        weights: impl AsRef<[usize]>,
        alignment: Option<&[AlignedPair<POANodePos<Ix>>]>,
    ) -> Result<(), GraphError<Ix>> {
        let span = debug_span!("add_aligned_sequence");
        let _enter = span.enter();
        
        let seq = sequence.as_ref();
        let w = weights.as_ref();
        let mut nodes_split = SplitTracker::default();

        assert_eq!(seq.len(), w.len());

        if self.is_empty() {
            // Make sure mutable borrow of the graph ends at the block, before calling post process.
            {
                let new_node = self.add_node_with_weights(sequence.as_ref().to_vec(), weights.as_ref().to_vec());
                self.sequences.push(Sequence(sequence_name.to_string(), new_node));
            }
            
            self.post_process(&nodes_split)?;

            return Ok(());
        }

        if !self.is_empty() && alignment.is_none() {
            return Err(GraphError::EmptyAlignment);
        }

        let aln = alignment.unwrap();

        // Process the alignment and break it up in multiple blocks. A block is either a
        // sequence of matches, mismatches, insertions, or deletions.
        let mut blocks = AlignmentBlocks::new(aln, self.start_node);
        let mut pred = None;
        let mut start_node = None;
        while let Some(block) = blocks.next_block(self, seq, &nodes_split) {
            match block.block_type {
                AlignmentBlockType::Match { qry_range, node_ivals } => {
                    if start_node.is_none() {
                        start_node = Some(node_ivals.first().unwrap().0.node());
                    }
                    
                    pred = Some(self.process_matches(w, qry_range, &node_ivals, pred, &mut nodes_split));
                },
                AlignmentBlockType::MismatchNew { qry_range, node_ivals } => {
                    pred = Some(self.process_mismatches_new(seq, w, qry_range, &node_ivals, pred));
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                AlignmentBlockType::MismatchPresent { qry_range, other_node } => {
                    pred = Some(self.process_mismatches_present(w, qry_range, other_node, pred, &mut nodes_split));
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                AlignmentBlockType::Insertion { qry_range } => {
                    pred = Some(self.process_insertion(seq, w, qry_range, pred));
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                AlignmentBlockType::Deletion { node_ivals } => {
                    pred = self.process_deletion(&node_ivals, pred, &mut nodes_split);
                }
            }
        }
        
        self.sequences.push(Sequence(sequence_name.to_string(), start_node.unwrap()));
        self.post_process(&nodes_split)?;

        Ok(())
    }
    
    fn process_matches(
        &mut self, 
        weights: &[usize], 
        qry_range: Range<usize>,
        node_ivals: &[(POANodePos<Ix>, usize)], 
        pred: Option<POANodeIndex<Ix>>,
        nodes_split: &mut SplitTracker<Ix>,
    ) -> POANodeIndex<Ix> {
        let new_seq_id = self.sequences.len();
        
        let Bound::Included(qry_start) = qry_range.start_bound() else {
            panic!("Invalid query start position");
        };
        
        // Update node weights with new alignment
        let mut curr_qry_pos = *qry_start;
        for (node_pos, length) in node_ivals {
            debug!("Adding weights for {:?}, {}..{} (length: {})", node_pos.node(), node_pos.pos(), node_pos.pos() + *length, length);
            let node_weights = &mut self.graph.node_weight_mut(node_pos.node()).unwrap().weights;
            
            let range = node_pos.pos()..node_pos.pos()+*length;
            node_weights[range].iter_mut()
                .zip(weights[curr_qry_pos..curr_qry_pos+*length].iter())
                .for_each(|(nw, w)| *nw += *w);
            
            curr_qry_pos += *length;
        }
        
        // Update edge sequence IDs if the aligned block spans multiple nodes in the graph
        for window in node_ivals.windows(2) {
            let (n1, _) = window[0];
            let (n2, _) = window[1];
            
            let edge = self.graph.find_edge(n1.node(), n2.node()).unwrap();
            self.graph.edge_weight_mut(edge).unwrap().sequence_ids.push(new_seq_id);
        }
        
        let (first_node_pos, _) = node_ivals.first().unwrap();
        if first_node_pos.pos() > 0 {
            // Block of matches starts in the middle of a node, split the node
            let (left_node, right_node) = self.split_node_at(first_node_pos.node(), first_node_pos.pos());
            nodes_split.add_split(first_node_pos.node(), first_node_pos.pos(), left_node, right_node);
        } 
        
        // Add edge from predecessor if given
        if let Some(p) = pred {
            let new = nodes_split.to_new_node_pos(*first_node_pos);
            self.graph.add_edge(p, new.node(), POAEdgeData::new_with_seq_id(new_seq_id));
        }
        
        let (last_node_pos, ival_length) = node_ivals.last().unwrap();
        let last_poa_node_pos = nodes_split.to_new_node_pos(POANodePos(last_node_pos.node(), last_node_pos.pos() + *ival_length));
        let last_node_seq_len = self.graph.node_weight(last_poa_node_pos.node()).unwrap().sequence.len();
        
        if last_poa_node_pos.pos() < last_node_seq_len {
            // Block of matches ends in the middle of a node, split the node
            let (left_node, right_node) = self.split_node_at(last_poa_node_pos.node(), last_poa_node_pos.pos());
            nodes_split.add_split(last_poa_node_pos.node(), last_poa_node_pos.pos(), left_node, right_node);
            
            left_node
        } else {
            // Block of matches ends exactly at the end of a node, no need to split
            last_poa_node_pos.node()
        }
    }
    
    fn process_mismatches_new(
        &mut self, 
        seq: &[u8], 
        weights: &[usize], 
        qry_range: Range<usize>,
        node_ivals: &[(POANodePos<Ix>, usize)], 
        pred: Option<POANodeIndex<Ix>>,
    ) -> POANodeIndex<Ix> {
        debug!("Process mismatches, node_ivals: {:?}", node_ivals);
        
        let new_seq_id = self.sequences.len();
        let mut new_qry_node_data = POANodeData::new_with_weights(
            seq[qry_range.clone()].to_vec(),   // Freakin' Range not being Copy...
            weights[qry_range.clone()].to_vec()
        );
        
        // Store aligned intervals to other nodes
        let mut curr_node_pos = 0;
        for (node_pos, length) in node_ivals {
            let aligned_interval = AlignedInterval::new(curr_node_pos, *length, node_pos.node(), node_pos.pos());
            debug!("Add aligned interval new node to -> {:?}, {}", node_pos.node(), node_pos.pos());
            debug!("{:?}", aligned_interval);
            new_qry_node_data.aligned_intervals.push(aligned_interval);
            
            curr_node_pos += *length;
        }
        
        let new_qry_node_ix = self.add_node_with_data(new_qry_node_data);
        debug!("Added node with sequence length {} (new node ix: {:?})", qry_range.len(), new_qry_node_ix);
        
        // Make sure to add aligned intervals pointing to the new node in other nodes
        curr_node_pos = 0;
        let mut additional_aln_ivals_for_qry = Vec::new();
        for (node_pos, length) in node_ivals {
            let node_data = self.graph.node_weight(node_pos.node()).unwrap();
            let aln_ival_to_qry = AlignedInterval::new(node_pos.pos(), *length, new_qry_node_ix, curr_node_pos);
            
            debug!("Add aligned interval aligned node {:?}, {} to -> query {:?}, {}", node_pos.node(), node_pos.pos(), new_qry_node_ix, curr_node_pos);
            debug!("{:?}", aln_ival_to_qry);
            
            // Our query is aligned to `node`, but `node` is potentially already aligned to other nodes.
            // Make sure we update their `aligned_intervals` too.
            let mut other_aln_ivals_to_new = Vec::new();
            for aln_ival_to_other in &node_data.aligned_intervals {
                debug!("Checking for intersection of: {:?}", aln_ival_to_other);
                
                if let Some((ival_other_to_query, ival_query_to_other)) = aln_ival_to_other.transitive_intersect(&aln_ival_to_qry) {
                    // Node was already aligned to 'other', and the aligned interval intersects with the 
                    // aligned interval to the new query.
                    debug!("Transitive alignment other {:?} -> qry {:?}", aln_ival_to_other.other(), new_qry_node_ix);
                    debug!("{:?}", ival_other_to_query);
                    other_aln_ivals_to_new.push((aln_ival_to_other.other(), ival_other_to_query));
                    
                    debug!("Transitive alignment qry {:?} -> other {:?}", new_qry_node_ix, aln_ival_to_other.other());
                    debug!("{:?}", ival_query_to_other);
                    additional_aln_ivals_for_qry.push(ival_query_to_other);
                    
                }
            }
            
            for (other_node, new_aln_ival_for_other) in other_aln_ivals_to_new {
                let other_node_data = self.graph.node_weight_mut(other_node).unwrap();
                other_node_data.aligned_intervals.push(new_aln_ival_for_other);
            }
            
            let node_data = self.graph.node_weight_mut(node_pos.node()).unwrap();
            node_data.aligned_intervals.push(aln_ival_to_qry);
            curr_node_pos += *length;
        }
        
        let qry_node_data = self.graph.node_weight_mut(new_qry_node_ix).unwrap();
        for aln_ival in additional_aln_ivals_for_qry {
            debug!("Add additional aligned interval query {:?}, {} to -> other {:?}, {}", new_qry_node_ix, aln_ival.node_start(), aln_ival.other(), aln_ival.other_start());
            qry_node_data.aligned_intervals.push(aln_ival);
        }
        
        // Create edge to new node from predecessor.
        if let Some(from_node) = pred {
            self.graph.add_edge(from_node, new_qry_node_ix, POAEdgeData::new_with_seq_id(new_seq_id));
        }
        
        new_qry_node_ix
    }
    
    fn process_mismatches_present(
        &mut self,
        weights: &[usize],
        qry_range: Range<usize>,
        other_node: POANodePos<Ix>,
        pred: Option<POANodeIndex<Ix>>,
        nodes_split: &mut SplitTracker<Ix>,
    ) -> POANodeIndex<Ix> {
        let new_seq_id = self.sequences.len();
        
        // Update node weights
        let query_len = qry_range.len();
        let node_range = other_node.pos()..other_node.pos() + query_len;
        
        let node_data = self.graph.node_weight_mut(other_node.node()).unwrap();
        node_data.weights[node_range].iter_mut()
            .zip(weights[qry_range.clone()].iter())
            .for_each(|(v, w)| *v += *w);
        
        if let Some(p) = pred {
            if other_node.pos() == 0 {
                // Matching sequencing at the beginning of the node, no need to split
                if let Some(e) = self.graph.find_edge(p, other_node.node()) {
                    let edge_data = self.graph.edge_weight_mut(e).unwrap();
                    edge_data.sequence_ids.push(new_seq_id);
                } else {
                    let edge_data = POAEdgeData::new_with_seq_id(new_seq_id);
                    self.graph.add_edge(p, other_node.node(), edge_data);
                }
            } else {
                // Split the node at the mismatch start position
                let (left_node, right_node) = self.split_node_at(other_node.node(), other_node.pos());
                nodes_split.add_split(other_node.node(), other_node.pos(), left_node, right_node);
                
                // Add edge from predecessor to the new right node
                let edge_data = POAEdgeData::new_with_seq_id(new_seq_id);
                self.graph.add_edge(p, right_node, edge_data);
            }
        }
        
        let end_node = nodes_split.to_new_node_pos(other_node);
        let end_node_seq_len = self.graph.node_weight_mut(end_node.node()).unwrap().sequence.len();
        let mut to_return = end_node.node();
        if end_node.pos() < end_node_seq_len - 1 {
            // End of the matching sequence is not at the end of the node, split it
            let (left_node, right_node) = self.split_node_at(end_node.node(), end_node.pos()+1);
            nodes_split.add_split(end_node.node(), end_node.pos()+1, left_node, right_node);
            
            to_return = left_node;
        }
        
        to_return
    }
    
    fn process_insertion(
        &mut self,
        seq: &[u8],
        weights: &[usize],
        qry_range: Range<usize>,
        pred: Option<POANodeIndex<Ix>>,
    ) -> POANodeIndex<Ix> {
        let new_seq_id = self.sequences.len();
        let new_node_data = POANodeData::new_with_weights(
            seq[qry_range.clone()].to_vec(),
            weights[qry_range.clone()].to_vec(),
        );
        let new_node = self.add_node_with_data(new_node_data);
        
        if let Some(p) = pred {
            let edge_data = POAEdgeData::new_with_seq_id(new_seq_id);
            self.graph.add_edge(p, new_node, edge_data);
        }
        
        new_node
    }
    
    fn process_deletion(
        &mut self,
        node_ivals: &[(POANodePos<Ix>, usize)],
        pred: Option<POANodeIndex<Ix>>,
        nodes_split: &mut SplitTracker<Ix>,
    ) -> Option<POANodeIndex<Ix>> {
        let (del_end_pos, length) = node_ivals.last().copied().unwrap();
        let del_end_pos = POANodePos::<Ix>(del_end_pos.node(), del_end_pos.pos() + length);
        
        let end_node_data = self.graph.node_weight(del_end_pos.node()).unwrap();
        let end_node_seq_len = end_node_data.sequence.len();
        
        if del_end_pos.pos() < end_node_seq_len {
            // End of the deletion is not at the end of the node, split it
            let (left_node, right_node) = self.split_node_at(del_end_pos.node(), del_end_pos.pos());
            nodes_split.add_split(del_end_pos.node(), del_end_pos.pos(), left_node, right_node);
        }
        
        pred
    }
    
    /// Add a new sequence to graph originating from a FASTA MSA
    fn add_msa_sequence(
        &mut self, 
        seq_name: &str, 
        seq: &[u8], 
        msa_nodes_cover: &mut MSANodeCover<Ix>,
    ) -> Result<(), GraphError<Ix>> {
        // First transform aligned sequence to blocks
        let aln: Vec<_> = seq.iter().enumerate()
            .map(|(col, &nuc)| {
                if nuc == b'-' {
                    (col, None, AlignmentClassification::Deletion)
                } else if !msa_nodes_cover.has_nodes_covering_col(col) {
                    (col, None, AlignmentClassification::Insertion)
                } else {
                    msa_nodes_cover.has_match(self, col, nuc)
                        .map(|node_pos| (col, Some(node_pos), AlignmentClassification::Match))
                        .unwrap_or((col, None, AlignmentClassification::Mismatch(None)))
                }
            })
            .collect();
        
        let mut curr_block_start = 0;
        let mut prev_state = AlignmentClassification::Start;
        let mut pred = None;
        let mut prev_match_node: Option<POANodeIndex<Ix>> = None;
        let mut start_node = None;
        let mut nodes_added = Vec::new();
        let mut nodes_split = SplitTracker::default();
        for (qpos, npos, curr_state) in &aln {
            match (prev_state, *curr_state) {
                (AlignmentClassification::Start, _) => {
                    // Start of the sequence, do nothing
                },
                
                (AlignmentClassification::Match, AlignmentClassification::Match) => {
                    // Check if the match is on the same node
                    if prev_match_node.is_some() && prev_match_node != npos
                        .map(|p| nodes_split.to_new_node_pos(p).node())
                    {
                        // Matches on a different node, start a new block
                        pred = Some(self.process_msa_matches(
                            &aln[curr_block_start..*qpos],
                            pred, 
                            &mut nodes_split,
                            msa_nodes_cover,
                        ));
                        
                        curr_block_start = *qpos;
                        
                        if start_node.is_none() {
                            start_node = pred;
                        }
                    }
                },
                
                (AlignmentClassification::Mismatch(_), AlignmentClassification::Mismatch(_)) 
                | (AlignmentClassification::Deletion, AlignmentClassification::Deletion)
                | (AlignmentClassification::Insertion, AlignmentClassification::Insertion)
                => {
                    // Simply continue
                }
                
                // Processed a block of matches
                (AlignmentClassification::Match, AlignmentClassification::Mismatch(_))
                | (AlignmentClassification::Match, AlignmentClassification::Deletion)
                | (AlignmentClassification::Match, AlignmentClassification::Insertion)
                => {
                    pred = Some(self.process_msa_matches(&aln[curr_block_start..*qpos], pred, &mut nodes_split, msa_nodes_cover));
                    curr_block_start = *qpos;
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                
                // Processed a block of mismatches
                (AlignmentClassification::Mismatch(_), AlignmentClassification::Match)
                | (AlignmentClassification::Mismatch(_), AlignmentClassification::Deletion)
                | (AlignmentClassification::Mismatch(_), AlignmentClassification::Insertion)
                => {
                    pred = Some(self.process_msa_mismatches(
                        &seq[curr_block_start..*qpos], 
                        &aln[curr_block_start..*qpos],
                        pred, 
                        &mut nodes_added,
                        msa_nodes_cover
                    ));
                    curr_block_start = *qpos;
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                
                (AlignmentClassification::Insertion, AlignmentClassification::Match)
                | (AlignmentClassification::Insertion, AlignmentClassification::Mismatch(_))
                | (AlignmentClassification::Insertion, AlignmentClassification::Deletion)
                => {
                    pred = Some(self.process_msa_insertion(
                        &seq[curr_block_start..*qpos], 
                        &aln[curr_block_start..*qpos], 
                        pred,
                        &mut nodes_added
                    ));
                    curr_block_start = *qpos;
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                
                (AlignmentClassification::Deletion, AlignmentClassification::Match)
                | (AlignmentClassification::Deletion, AlignmentClassification::Mismatch(_))
                | (AlignmentClassification::Deletion, AlignmentClassification::Insertion)
                => {
                    // Ignore deletions
                    curr_block_start = *qpos;
                },
                    
                (_, AlignmentClassification::Start) => 
                    panic!("Invalid state transition to AlignmentClassification::Start.")
            }
            
            if *curr_state == AlignmentClassification::Match {
                prev_match_node = Some(nodes_split.to_new_node_pos(npos.unwrap()).node());
            } else {
                prev_match_node = None;
            }
            
            prev_state = *curr_state;
        }
        
        // Process the last block
        if curr_block_start < aln.len() {
            match aln[curr_block_start].2 {
                AlignmentClassification::Match => {
                    pred = Some(self.process_msa_matches(
                        &aln[curr_block_start..], 
                        pred, 
                        &mut nodes_split,
                        msa_nodes_cover
                    ));
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                
                AlignmentClassification::Mismatch(_) => {
                    pred = Some(self.process_msa_mismatches(
                        &seq[curr_block_start..], 
                        &aln[curr_block_start..], 
                        pred, 
                        &mut nodes_added,
                        msa_nodes_cover
                    ));
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                
                AlignmentClassification::Insertion => {
                    pred = Some(self.process_msa_insertion(
                        &seq[curr_block_start..], 
                        &aln[curr_block_start..], 
                        pred, 
                        &mut nodes_added
                    ));
                    
                    if start_node.is_none() {
                        start_node = pred;
                    }
                },
                
                AlignmentClassification::Deletion => {
                    // Ignore deletions
                },
                
                AlignmentClassification::Start => 
                    panic!("Invalid state at the end of the alignment.")
            }
        }
        
        self.sequences.push(Sequence(seq_name.to_string(), start_node.unwrap()));
        
        // Update the nodes per MSA col with newly added nodes
        for (astart, aend, node) in nodes_added {
            msa_nodes_cover.add_node(node, astart..aend);
        }
        
        self.post_process(&nodes_split)?;
        
        Ok(())
    }
    
    fn process_msa_matches(
        &mut self, 
        aln: &[(usize, Option<POANodePos<Ix>>, AlignmentClassification<Ix>)], 
        pred: Option<POANodeIndex<Ix>>,
        nodes_split: &mut SplitTracker<Ix>,
        msa_nodes_cover: &mut MSANodeCover<Ix>
    ) -> POANodeIndex<Ix> {
        let new_seq_id = self.sequences.len();
        
        let node_pos_start = nodes_split.to_new_node_pos(aln[0].1.unwrap());
        if node_pos_start.pos() > 0 {
            let (left_node, right_node) = self.split_node_at(node_pos_start.node(), node_pos_start.pos());
            nodes_split.add_split(node_pos_start.node(), node_pos_start.pos(), left_node, right_node);
            msa_nodes_cover.split_node(node_pos_start.node(), node_pos_start.pos(), left_node, right_node);
        }
        
        if let Some(p) = pred {
            let new = nodes_split.to_new_node_pos(node_pos_start);
            if let Some(e) = self.graph.find_edge(p, new.node()) {
                self.graph.edge_weight_mut(e).unwrap().sequence_ids.push(new_seq_id);
            } else {
                self.graph.add_edge(p, new.node(), POAEdgeData::new_with_seq_id(new_seq_id));
            }
        }
        
        let last_node_pos = nodes_split.to_new_node_pos(aln.last().unwrap().1.unwrap());
        let last_node_seq_len = self.graph.node_weight(last_node_pos.node()).unwrap().sequence.len();
        if last_node_pos.pos() < last_node_seq_len - 1 {
            let (left_node, right_node) = self.split_node_at(last_node_pos.node(), last_node_pos.pos()+1);
            nodes_split.add_split(last_node_pos.node(), last_node_pos.pos() + 1, left_node, right_node);
            msa_nodes_cover.split_node(last_node_pos.node(), last_node_pos.pos() + 1, left_node, right_node);
            
            left_node
        } else {
            last_node_pos.node()
        }
    }
    
    fn process_msa_mismatches(
        &mut self, 
        seq: &[u8],
        aln: &[(usize, Option<POANodePos<Ix>>, AlignmentClassification<Ix>)], 
        pred: Option<POANodeIndex<Ix>>,
        nodes_added: &mut Vec<(usize, usize, POANodeIndex<Ix>)>,
        msa_nodes_cover: &mut MSANodeCover<Ix>,
    ) -> POANodeIndex<Ix> {
        let new_seq_id = self.sequences.len();
        let mut new_node_data = POANodeData::new(seq.to_vec());
        
        // Construct aligned intervals pointing to other nodes
        let mut aln_ivals: HashMap<POANodeIndex<Ix>, AlignedInterval<Ix>> = HashMap::default();
        for (i, (qpos, _, _)) in aln.iter().enumerate() {
            for node_pos in msa_nodes_cover.get_nodes_for_col(*qpos) {
                aln_ivals.entry(node_pos.node())
                    .and_modify(|v| v.increase_length(1))
                    .or_insert_with(|| {
                        AlignedInterval::new(i, 1, node_pos.node(), node_pos.pos())
                    });
            }
        }
        
        for ival in aln_ivals.values() {
            new_node_data.aligned_intervals.push(*ival);
        }
        
        let new_node_ix = self.add_node_with_data(new_node_data);
        
        // Add aligned intervals in other nodes pointing to the new node
        for (other, ival) in aln_ivals.iter() {
            let other_node_data = self.graph.node_weight_mut(*other).unwrap();
            
            let flipped = ival.flip_nodes(new_node_ix);
            other_node_data.aligned_intervals.push(flipped);
        }

        // Update which nodes cover which portions of the MSA
        nodes_added.push((aln[0].0, aln.last().unwrap().0 + 1, new_node_ix));
        
        if let Some(p) = pred {
            self.graph.add_edge(p, new_node_ix, POAEdgeData::new_with_seq_id(new_seq_id));
        }
        
        new_node_ix
    }
    
    fn process_msa_insertion(
        &mut self,
        seq: &[u8],
        aln: &[(usize, Option<POANodePos<Ix>>, AlignmentClassification<Ix>)],
        pred: Option<POANodeIndex<Ix>>,
        nodes_added: &mut Vec<(usize, usize, POANodeIndex<Ix>)>
    ) -> POANodeIndex<Ix> {
        let new_seq_id = self.sequences.len();
        let new_node_data = POANodeData::new(seq.to_vec());
        let new_node = self.add_node_with_data(new_node_data);
        
        if let Some(p) = pred {
            self.graph.add_edge(p, new_node, POAEdgeData::new_with_seq_id(new_seq_id));
        }
        
        nodes_added.push((aln[0].0, aln.last().unwrap().0 + 1, new_node));
        
        new_node
    }
    
    fn post_process(&mut self, nodes_split: &SplitTracker<Ix>) -> Result<(), GraphError<Ix>> {
        let mut to_remove = Vec::new();
        for e in self.graph.edges(self.start_node) {
            let in_degree_target = self.graph.edges_directed(e.target(), Incoming).count();
            
            if in_degree_target > 1 {
                // Remove edges with multiple incoming edges
                to_remove.push(e.id());
            }
        }
        
        for e in self.graph.edges_directed(self.end_node, Incoming) {
            let out_degree_src = self.graph.edges_directed(e.source(), Outgoing).count();
            
            if out_degree_src > 1 {
                to_remove.push(e.id());
            }
        }
        
        for e in to_remove {
            self.graph.remove_edge(e);
        }
        
        let mut to_add = Vec::new();
        for n in self.graph.node_indices() {
            if n == self.start_node || n == self.end_node {
                continue;
            }
            
            // Add new edges from start node to nodes with no incoming edges
            if self.graph.edges_directed(n, Incoming).count() == 0 {
                let mut all_outgoing_ids: Vec<_> = self.graph.edges_directed(n, Outgoing)
                    .flat_map(|e| self.graph[e.id()].sequence_ids.iter().copied())
                    .collect();
                
                all_outgoing_ids.sort();
                
                to_add.push((self.start_node, n, all_outgoing_ids));
            }
            
            // Add new edges from nodes with no outgoing edges to the end node
            if self.graph.edges_directed(n, Outgoing).count() == 0 {
                let mut all_incoming_ids: Vec<_> = self.graph.edges_directed(n, Incoming)
                    .flat_map(|e| self.graph[e.id()].sequence_ids.iter().copied())
                    .collect();
                
                all_incoming_ids.sort();
                to_add.push((n, self.end_node, all_incoming_ids));
            }
        }
        
        for e in to_add {
            debug!("Adding edge {:?}", e);
            self.graph.add_edge(e.0, e.1, POAEdgeData::new_with_seq_ids(e.2));
        }
        
        self.toposorted.clear();
        self.toposorted = toposort(&self.graph, None)?;
        
        for (rank, n) in self.toposorted.iter().enumerate() {
            self.graph.node_weight_mut(*n).unwrap().rank = Ix::new(rank);
        }
        
        // Fix stored start node IDs to account for possible node splits
        for seq in self.sequences.iter_mut() {
            seq.1 = nodes_split.to_new_node_pos(POANodePos(seq.1, 0)).node();
        }
        
        Ok(())
    }
        
    fn split_node_at(
        &mut self,
        node: POANodeIndex<Ix>,
        pos: usize,
    ) -> (POANodeIndex<Ix>, POANodeIndex<Ix>) {
        let node_data = self.graph.node_weight(node).unwrap();
        let (left_seq, right_seq) = node_data.sequence.split_at(pos);
        let (left_weights, right_weights) = node_data.weights.split_at(pos);
        
        let (left_len, right_len) = (left_seq.len(), right_seq.len());
        
        let (left_seq, right_seq) = (left_seq.to_vec(), right_seq.to_vec());
        let (left_weights, right_weights) = (left_weights.to_vec(), right_weights.to_vec());
        
        let left_node = self.add_node_with_weights(left_seq.to_vec(), left_weights.to_vec());
        let right_node = self.add_node_with_weights(right_seq.to_vec(), right_weights.to_vec());
        
        debug!("Split node {:?} at position {}. New nodes: {:?} (len: {}), {:?} (len: {})", 
            node, pos, left_node, left_len, right_node, right_len);
        
        // Re-borrow node data such that the borrow checker is happy. Otherwise, the lifetime of the previous `node_data` borrow
        // would also span the above `add_node_with_weights` calls, which require mutable borrows to the graph.
        let node_data = self.graph.node_weight(node).unwrap();
        
        // Find intervals that are completely left of the split position
        let mut intervals_left: Vec<_> = node_data.aligned_intervals
            .iter()
            .filter(|ival| ival.node_end() <= pos)
            .copied()
            .collect();
        
        // Extend with overlapping intervals, modify interval end position
        intervals_left.extend(
            node_data.aligned_intervals
                .iter()
                .filter(|ival| ival.node_start() < pos && ival.node_end() > pos)
                .map(|ival| ival.split_at_left(pos))
        );
        
        // Find overlapping intervals, modify interval start position
        let mut intervals_right: Vec<_> = node_data.aligned_intervals
            .iter()
            .filter(|ival| ival.node_start() < pos && ival.node_end() > pos)
            .map(|ival| ival.split_at_right(pos))
            .collect();
        
        // Extend with intervals completely to the right of the split position
        intervals_right.extend(
            node_data.aligned_intervals
                .iter()
                .filter(|ival| ival.node_start() >= pos)
                .map(|ival| ival.move_node_start(pos))
        );
        
        // Fix intervals of other nodes, pointing to the split node
        let mut all_other_nodes_updates = Vec::default();
        for ival in &node_data.aligned_intervals {
            let other_node_data = self.graph.node_weight(ival.other()).unwrap();
            
            let mut new_other_ivals = Vec::default();
            other_node_data.aligned_intervals
                .iter()
                .for_each(
                    |other_ival| {
                        if other_ival.other() == node {
                            if other_ival.other_end() <= pos {
                                // Interval in the other node pointing to the old split node, completely left of the split position
                                new_other_ivals.push(other_ival.with_new_other(left_node))
                            } else if other_ival.other_start() >= pos {
                                // Interval in the other node pointing to the old split node, completely right of the split position
                                new_other_ivals.push(other_ival.move_other_start(right_node, pos))
                            } else {
                                // Interval in the other node pointing to the old split node, overlapping the split position
                                new_other_ivals.push(other_ival.other_split_at_left(left_node, pos));
                                new_other_ivals.push(other_ival.other_split_at_right(right_node, pos));
                            }
                        } else {
                            new_other_ivals.push(*other_ival);
                        }
                    }
                );
            
            all_other_nodes_updates.push((ival.other(), new_other_ivals));
        }
        
        self.graph.node_weight_mut(left_node).unwrap().aligned_intervals = intervals_left;
        self.graph.node_weight_mut(right_node).unwrap().aligned_intervals = intervals_right;
        
        for (other_node, new_aln_ivals) in all_other_nodes_updates {
            let other_node_data = self.graph.node_weight_mut(other_node).unwrap();
            other_node_data.aligned_intervals = new_aln_ivals;
        }
        
        // Add internal edge from left node -> right node retaining all incoming sequence ids except the newly added sequence
        // Use itertools::kmerge to ensure edge sequence IDs remain sorted
        let seq_ids: Vec<_> = itertools::kmerge(
            self.graph.edges_directed(node, Incoming)
                .map(|e| {
                    e.weight().sequence_ids.iter()
                        .copied()
                        .filter(|v| *v != self.sequences.len())
                })
        ).collect();
        debug!("SPLIT: seq ids in {:?}", seq_ids);
        
        self.graph
            .add_edge(left_node, right_node, POAEdgeData::new_with_seq_ids(seq_ids));
        
        // Add original incoming and outgoing edges to the left and right node, respectively
        let mut new_edges = Vec::default();
        for in_edge in self.graph.edges_directed(node, Incoming) {
            let edge_data = in_edge.weight();
            let new_edge_data = POAEdgeData::new_with_seq_ids(edge_data.sequence_ids.clone());
            debug!("OLD: {:?} -> {:?} NEW {:?} -> {:?} {:?}", in_edge.source(), node, in_edge.source(), left_node, new_edge_data);
            
            new_edges.push((in_edge.source(), left_node, new_edge_data));
        }
        
        for out_edge in self.graph.edges_directed(node, Outgoing) {
            let edge_data = out_edge.weight();
            let new_edge_data = POAEdgeData::new_with_seq_ids(edge_data.sequence_ids.clone());
            debug!("OLD: {:?} -> {:?} NEW {:?} -> {:?} {:?}", node, out_edge.target(), right_node, out_edge.target(), new_edge_data);
            
            new_edges.push((right_node, out_edge.target(), new_edge_data));
        }
        
        for (source, target, edge_data) in new_edges {
            self.graph.add_edge(source, target, edge_data);
        }
        
        // Delete the original node
        self.graph.remove_node(node);

        (left_node, right_node)
    }
}

impl<Ix> Default for POASeqGraph<Ix>
where
    Ix: petgraph::graph::IndexType,
{
    fn default() -> Self {
        Self::new()
    }
}

/// A sequence aligned to the POA graph.
///
/// Stores the sequence name and the start node in the graph.
#[derive(Debug, Clone)]
pub struct Sequence<Ix>(pub(crate) String, pub(crate) POANodeIndex<Ix>)
    where Ix: IndexType;

impl<Ix> Sequence<Ix>
    where Ix: IndexType
{
    pub fn name(&self) -> &str {
        &self.0
    }

    pub fn start_node(&self) -> POANodeIndex<Ix> {
        self.1
    }
}


impl<Ix> AlignableGraph for POASeqGraph<Ix> 
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    type NodePosType = POANodePos<Ix>;
    type Successors<'a> = Neighbors<'a, POAEdgeData, Ix>;
    type Predecessors<'a> = Neighbors<'a, POAEdgeData, Ix>;
    
    fn start_node(&self) -> Self::NodeType {
        self.start_node
    }
    
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }
    
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }
    
    fn node_bound(&self) -> usize {
        self.graph.node_bound()
    }
    
    fn node_seq(&self, node: Self::NodeType) -> &[u8] {
        &self.graph.node_weight(node).unwrap().sequence
    }
    
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.graph.neighbors_directed(node, Outgoing)
    }
    
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }
}

impl<'a, Ix> AlignableGraph for &'a POASeqGraph<Ix> 
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    type NodePosType = POANodePos<Ix>;
    type Successors<'b> = Neighbors<'b, POAEdgeData, Ix>
        where 'a: 'b;
    type Predecessors<'b> = Neighbors<'b, POAEdgeData, Ix>
        where 'a: 'b;
    
    fn start_node(&self) -> Self::NodeType {
        self.start_node
    }
    
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }
    
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }
    
    fn node_bound(&self) -> usize {
        self.graph.node_bound()
    }
    
    fn node_seq(&self, node: Self::NodeType) -> &[u8] {
        &self.graph.node_weight(node).unwrap().sequence
    }
    
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.graph.neighbors_directed(node, Outgoing)
    }
    
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }
}

impl<'a, Ix> AlignableGraph for &'a mut POASeqGraph<Ix> 
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    type NodePosType = POANodePos<Ix>;
    type Successors<'b> = Neighbors<'b, POAEdgeData, Ix>
        where 'a: 'b;
    type Predecessors<'b> = Neighbors<'b, POAEdgeData, Ix>
        where 'a: 'b;
    
    fn start_node(&self) -> Self::NodeType {
        self.start_node
    }
    
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }
    
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }
    
    fn node_bound(&self) -> usize {
        self.graph.node_bound()
    }
    
    fn node_seq(&self, node: Self::NodeType) -> &[u8] {
        &self.graph.node_weight(node).unwrap().sequence
    }
    
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.graph.neighbors_directed(node, Outgoing)
    }
    
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }
}


impl<T> AlignableGraphNodeId for T
where
    T: IndexType,
{
    #[inline(always)]
    fn index(&self) -> usize {
        self.index()
    }
}

/// Keep track of split nodes while processing a new alignment.
///
/// This enables to convert node positions in the original graph to the new graph after a split.
#[derive(Default)]
pub(crate) struct SplitTracker<Ix>
where 
    Ix: IndexType
{
    split_nodes: HashMap<POANodeIndex<Ix>, (usize, POANodeIndex<Ix>, POANodeIndex<Ix>)>
}

impl<Ix> SplitTracker<Ix> 
where 
    Ix: IndexType,
{
    pub(crate) fn add_split(&mut self, node: POANodeIndex<Ix>, split_pos: usize, left: POANodeIndex<Ix>, right: POANodeIndex<Ix>) {
        self.split_nodes.insert(node, (split_pos, left, right));
    }

    pub(crate) fn to_new_node_pos(&self, mut node_pos: POANodePos<Ix>) -> POANodePos<Ix> {
        while let Some((split_pos, left, right)) = self.split_nodes.get(&node_pos.node()) {
            if node_pos.pos() < *split_pos {
                node_pos = POANodePos(*left, node_pos.pos());
            } else {
                node_pos = POANodePos(*right, node_pos.pos() - *split_pos);
            }
        }
        
        node_pos
    }
}


#[cfg(test)]
mod tests {
    use std::{fs::File, io::{self, BufReader}};

    use tracing::{span, Level};
    use noodles::fasta;

    use crate::{aligner::utils::AlignedPair, graph::{alignment::POANodePos, poa::{graph_impl::AlignedInterval, SplitTracker}}};

    use super::POASeqGraph;
    
    #[test]
    fn test_split_node() {
        let mut graph = POASeqGraph::<usize>::new();
        let orig_node_seq = b"GTCTGCTATACTGCGTACGTCGT";
        let orig_node = graph.add_node(orig_node_seq.to_vec());
        let orig_node2 = graph.add_node(orig_node_seq.to_vec());
        
        // Add some test aligned intervals
        let aln_ivals = [
            AlignedInterval::new(0, 5, orig_node2, 2),
            AlignedInterval::new(8, 6, orig_node2, 10),
            AlignedInterval::new(19, 3, orig_node2, 20),
        ];
        
        let aln_ivals2 = aln_ivals.iter()
            .map(|v| v.flip_nodes(orig_node))
            .collect::<Vec<_>>();
        
        {
            let node_data = graph.node_data_mut(orig_node);
            node_data.aligned_intervals.extend(&aln_ivals);
        }
        
        {
            let node_data2 = graph.node_data_mut(orig_node2);
            node_data2.aligned_intervals.extend(aln_ivals2);
        }
        
        // Do the split
        let (lnode, rnode) = {
            graph.split_node_at(orig_node, 10)
        };
        
        assert_eq!(graph.num_nodes(), 3);
        
        // Test whether the sequences are split correctly
        let lnode_data = graph.node_data(lnode);
        let rnode_data = graph.node_data(rnode);
        
        assert_eq!(&lnode_data.sequence, &orig_node_seq[..10]);
        assert_eq!(&rnode_data.sequence, &orig_node_seq[10..]);
        
        // Test whether the aligned intervals are split correctly
        assert_eq!(&lnode_data.aligned_intervals, &[
            AlignedInterval::new(0, 5, orig_node2, 2),
            AlignedInterval::new(8, 2, orig_node2, 10),
        ]);
        
        assert_eq!(&rnode_data.aligned_intervals, &[
            AlignedInterval::new(0, 4, orig_node2, 12),
            AlignedInterval::new(9, 3, orig_node2, 20),
        ]);
        
        // Test updated aligned intervals in the other node, pointing to the new split nodes
        let orig_node2_data = graph.node_data(orig_node2);
        
        assert_eq!(&orig_node2_data.aligned_intervals, &[
            AlignedInterval::new(2, 5, lnode, 0),
            AlignedInterval::new(10, 2, lnode, 8),
            AlignedInterval::new(12, 4, rnode, 0),
            AlignedInterval::new(20, 3, rnode, 9),
        ]);
    }
    
    #[test]
    fn test_add_alignment_mismatches() {
        let mut graph = POASeqGraph::<usize>::new();
        let orig_node_seq = b"GTCTGCTATACTGCGTACGTCGT";
        let orig_node = graph.add_node(orig_node_seq.to_vec());
        let nodes_split = SplitTracker::default();
        graph.post_process(&nodes_split).unwrap();
        
        let seq2 = b"GTCTGCTATGGGGCGTACGTCGT";
        let weights2 = vec![1; seq2.len()];
        let alignment = vec![
            AlignedPair::new(Some(POANodePos(orig_node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(orig_node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(orig_node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(orig_node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(orig_node, 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(orig_node, 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(orig_node, 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(orig_node, 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(orig_node, 8)), Some(8)),
            AlignedPair::new(Some(POANodePos(orig_node, 9)), Some(9)),
            AlignedPair::new(Some(POANodePos(orig_node, 10)), Some(10)),
            AlignedPair::new(Some(POANodePos(orig_node, 11)), Some(11)),
            AlignedPair::new(Some(POANodePos(orig_node, 12)), Some(12)),
            AlignedPair::new(Some(POANodePos(orig_node, 13)), Some(13)),
            AlignedPair::new(Some(POANodePos(orig_node, 14)), Some(14)),
            AlignedPair::new(Some(POANodePos(orig_node, 15)), Some(15)),
            AlignedPair::new(Some(POANodePos(orig_node, 16)), Some(16)),
            AlignedPair::new(Some(POANodePos(orig_node, 17)), Some(17)),
            AlignedPair::new(Some(POANodePos(orig_node, 18)), Some(18)),
            AlignedPair::new(Some(POANodePos(orig_node, 19)), Some(19)),
            AlignedPair::new(Some(POANodePos(orig_node, 20)), Some(20)),
            AlignedPair::new(Some(POANodePos(orig_node, 21)), Some(21)),
            AlignedPair::new(Some(POANodePos(orig_node, 22)), Some(22)),
        ];
        
        graph.add_aligned_sequence("seq2", seq2, &weights2, Some(&alignment)).unwrap();
        
        let seq_truth: Vec<&[u8]> = vec![
            b"#",
            b"GTCTGCTAT",
            b"ACT",
            b"GGG",
            b"GCGTACGTCGT",
            b"$"
        ];
        
        let w_truth: Vec<&[usize]> = vec![
            &[1],
            &[2, 2, 2, 2, 2, 2, 2, 2, 2],
            &[1, 1, 1],
            &[1, 1, 1],
            &[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            &[1],
        ];
        
        for (rank, n) in graph.toposorted.iter().enumerate() {
            let node_data = graph.node_data(*n);
            assert_eq!(&node_data.sequence, seq_truth[rank]);
            assert_eq!(&node_data.weights, w_truth[rank]);
        }
        
        let seq3 = b"GTCTGCTATGCGGCGTACGTCGT";
        let weights3 = vec![1; seq2.len()];
        let alignment = vec![
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 8)), Some(8)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[3], 0)), Some(9)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[3], 1)), Some(10)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[3], 2)), Some(11)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 0)), Some(12)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 1)), Some(13)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 2)), Some(14)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 3)), Some(15)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 4)), Some(16)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 5)), Some(17)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 6)), Some(18)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 7)), Some(19)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 8)), Some(20)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 9)), Some(21)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 10)), Some(22)),
        ];
        
        graph.add_aligned_sequence("seq3", seq3, &weights3, Some(&alignment)).unwrap();
        
        let seq_truth: Vec<&[u8]> = vec![
            b"#",
            b"GTCTGCTAT",
            b"A",
            b"G",
            b"G",
            b"C",
            b"G",
            b"T",
            b"GCGTACGTCGT",
            b"$",
        ];
        
        let w_truth: Vec<&[usize]> = vec![
            &[1],
            &[3, 3, 3, 3, 3, 3, 3, 3, 3],
            &[1],
            &[2],
            &[1],
            &[2],
            &[2],
            &[1],
            &[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            &[1],
        ];
        
        let aln_ival_truth: Vec<Vec<AlignedInterval<usize>>> = vec![
            vec![],
            vec![],
            vec![AlignedInterval::new(0, 1, graph.toposorted[3], 0)],
            vec![AlignedInterval::new(0, 1, graph.toposorted[2], 0)],
            vec![AlignedInterval::new(0, 1, graph.toposorted[5], 0)],
            vec![AlignedInterval::new(0, 1, graph.toposorted[4], 0)],
            vec![AlignedInterval::new(0, 1, graph.toposorted[7], 0)],
            vec![AlignedInterval::new(0, 1, graph.toposorted[6], 0)],
            vec![],
            vec![]
        ];
        
        for (rank, n) in graph.toposorted.iter().enumerate() {
            let node_data = graph.node_data(*n);
            assert_eq!(&node_data.sequence, seq_truth[rank]);
            assert_eq!(&node_data.weights, w_truth[rank]);
            assert_eq!(&node_data.aligned_intervals, &aln_ival_truth[rank]);
        }
    }
    
    #[test]
    fn test_add_alignment_mismatches_multinode() {
        let mut graph = POASeqGraph::<usize>::new();
        let seq1 = b"GTCTGCTATACTGCGTACGTCGT";
        let seq2 = b"GTCTGCTATGGGGCGTACGTCGT";
        let seq3 = b"GTCTGCTATGGAAAATACGTCGT";
        
        let orig_node = graph.add_node(seq1.to_vec());
        let nodes_split = SplitTracker::default();
        graph.post_process(&nodes_split).unwrap();
        
        let weights2 = vec![1; seq2.len()];
        let alignment = vec![
            AlignedPair::new(Some(POANodePos(orig_node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(orig_node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(orig_node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(orig_node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(orig_node, 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(orig_node, 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(orig_node, 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(orig_node, 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(orig_node, 8)), Some(8)),
            AlignedPair::new(Some(POANodePos(orig_node, 9)), Some(9)),
            AlignedPair::new(Some(POANodePos(orig_node, 10)), Some(10)),
            AlignedPair::new(Some(POANodePos(orig_node, 11)), Some(11)),
            AlignedPair::new(Some(POANodePos(orig_node, 12)), Some(12)),
            AlignedPair::new(Some(POANodePos(orig_node, 13)), Some(13)),
            AlignedPair::new(Some(POANodePos(orig_node, 14)), Some(14)),
            AlignedPair::new(Some(POANodePos(orig_node, 15)), Some(15)),
            AlignedPair::new(Some(POANodePos(orig_node, 16)), Some(16)),
            AlignedPair::new(Some(POANodePos(orig_node, 17)), Some(17)),
            AlignedPair::new(Some(POANodePos(orig_node, 18)), Some(18)),
            AlignedPair::new(Some(POANodePos(orig_node, 19)), Some(19)),
            AlignedPair::new(Some(POANodePos(orig_node, 20)), Some(20)),
            AlignedPair::new(Some(POANodePos(orig_node, 21)), Some(21)),
            AlignedPair::new(Some(POANodePos(orig_node, 22)), Some(22)),
        ];
        
        graph.add_aligned_sequence("seq2", seq2, &weights2, Some(&alignment)).unwrap();
        
        let weights3 = vec![1; seq2.len()];
        let alignment = vec![
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[1], 8)), Some(8)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[3], 0)), Some(9)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[3], 1)), Some(10)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[3], 2)), Some(11)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 0)), Some(12)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 1)), Some(13)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 2)), Some(14)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 3)), Some(15)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 4)), Some(16)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 5)), Some(17)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 6)), Some(18)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 7)), Some(19)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 8)), Some(20)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 9)), Some(21)),
            AlignedPair::new(Some(POANodePos(graph.toposorted[4], 10)), Some(22)),
        ];
        
        graph.add_aligned_sequence("seq3", seq3, &weights3, Some(&alignment)).unwrap();
        
        let seq_truth: Vec<&[u8]> = vec![
            b"#",
            b"GTCTGCTAT",
            b"GG",
            b"AAAA",
            b"G",
            b"ACT",
            b"GCG",
            b"TACGTCGT",
            b"$",
        ];
        
        let w_truth: Vec<&[usize]> = vec![
            &[1],
            &[3, 3, 3, 3, 3, 3, 3, 3, 3],
            &[2, 2],
            &[1, 1, 1, 1],
            &[1],
            &[1, 1, 1],
            &[2, 2, 2],
            &[3, 3, 3, 3, 3, 3, 3, 3],
            &[1],
        ];
        
        let aln_ival_truth: Vec<Vec<AlignedInterval<usize>>> = vec![
            vec![],
            vec![],
            vec![AlignedInterval::new(0, 2, graph.toposorted[5], 0)],
            vec![
                AlignedInterval::new(0, 1, graph.toposorted[4], 0),
                AlignedInterval::new(1, 3, graph.toposorted[6], 0),
                AlignedInterval::new(0, 1, graph.toposorted[5], 2),
            ],
            vec![
                AlignedInterval::new(0, 1, graph.toposorted[5], 2),
                AlignedInterval::new(0, 1, graph.toposorted[3], 0),
            ],
            vec![
                AlignedInterval::new(0, 2, graph.toposorted[2], 0),
                AlignedInterval::new(2, 1, graph.toposorted[4], 0),
                AlignedInterval::new(2, 1, graph.toposorted[3], 0),
            ],
            vec![AlignedInterval::new(0, 3, graph.toposorted[3], 1)],
            vec![],
            vec![],
        ];
        
        for (rank, n) in graph.toposorted.iter().enumerate() {
            let node_data = graph.node_data(*n);
            
            assert_eq!(&node_data.sequence, seq_truth[rank]);
            assert_eq!(&node_data.weights, w_truth[rank]);
            assert_eq!(&node_data.aligned_intervals, &aln_ival_truth[rank]);
            
        }
    }
    
    #[test]
    fn test_insertion() {
        let mut graph = POASeqGraph::<usize>::new();
        let seq1 = b"GTCTGCTATACTGCGTACGTCGT";
        let seq2 = b"GTCTGCTATTTAATACTGCGTACGTCGT";
        let orig_node = graph.add_node(seq1.to_vec());
        let nodes_split = SplitTracker::default();
        graph.post_process(&nodes_split).unwrap();
        
        let weights2 = vec![1; seq2.len()];
        let alignment = vec![
            AlignedPair::new(Some(POANodePos(orig_node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(orig_node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(orig_node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(orig_node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(orig_node, 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(orig_node, 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(orig_node, 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(orig_node, 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(orig_node, 8)), Some(8)),
            AlignedPair::new(None, Some(9)),
            AlignedPair::new(None, Some(10)),
            AlignedPair::new(None, Some(11)),
            AlignedPair::new(None, Some(12)),
            AlignedPair::new(None, Some(13)),
            AlignedPair::new(Some(POANodePos(orig_node, 9)), Some(14)),
            AlignedPair::new(Some(POANodePos(orig_node, 10)), Some(15)),
            AlignedPair::new(Some(POANodePos(orig_node, 11)), Some(16)),
            AlignedPair::new(Some(POANodePos(orig_node, 12)), Some(17)),
            AlignedPair::new(Some(POANodePos(orig_node, 13)), Some(18)),
            AlignedPair::new(Some(POANodePos(orig_node, 14)), Some(19)),
            AlignedPair::new(Some(POANodePos(orig_node, 15)), Some(20)),
            AlignedPair::new(Some(POANodePos(orig_node, 16)), Some(21)),
            AlignedPair::new(Some(POANodePos(orig_node, 17)), Some(22)),
            AlignedPair::new(Some(POANodePos(orig_node, 18)), Some(23)),
            AlignedPair::new(Some(POANodePos(orig_node, 19)), Some(24)),
            AlignedPair::new(Some(POANodePos(orig_node, 20)), Some(25)),
            AlignedPair::new(Some(POANodePos(orig_node, 21)), Some(26)),
            AlignedPair::new(Some(POANodePos(orig_node, 22)), Some(27)),
        ];
        
        graph.add_aligned_sequence("seq2", seq2, &weights2, Some(&alignment)).unwrap();
        
        let seq_truth: Vec<&[u8]> = vec![
            b"#",
            b"GTCTGCTAT",
            b"TTAAT",
            b"ACTGCGTACGTCGT",
            b"$",
        ];
        
        let w_truth: Vec<&[usize]> = vec![
            &[1],
            &[2, 2, 2, 2, 2, 2, 2, 2, 2],
            &[1, 1, 1, 1, 1],
            &[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            &[1],
        ];
        
        for (rank, n) in graph.toposorted.iter().enumerate() {
            let node_data = graph.node_data(*n);
            // println!("{} {:?} seq: {:?}", rank, n, std::str::from_utf8(node_data.sequence()).unwrap());
            // println!("{} {:?}   w: {:?}", rank, n, node_data.weights());
            // println!("{} {:?} aln: {:?}", rank, n, node_data.aligned_intervals());
            
            assert_eq!(&node_data.sequence, seq_truth[rank]);
            assert_eq!(&node_data.weights, w_truth[rank]);
        }
    }
    
    #[test]
    fn test_deletion() {
        let mut graph = POASeqGraph::<usize>::new();
        let seq1 = b"GTCTGCTATACTGCGTACGTCGT";
        let seq2 = b"GTCTGCTATGCGTACGTCGT"; 
        let orig_node = graph.add_node(seq1.to_vec());
        let nodes_split = SplitTracker::default();
        graph.post_process(&nodes_split).unwrap();
        
        let weights2 = vec![1; seq2.len()];
        let alignment = vec![
            AlignedPair::new(Some(POANodePos(orig_node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(orig_node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(orig_node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(orig_node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(orig_node, 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(orig_node, 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(orig_node, 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(orig_node, 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(orig_node, 8)), Some(8)),
            AlignedPair::new(Some(POANodePos(orig_node, 9)), None),
            AlignedPair::new(Some(POANodePos(orig_node, 10)), None),
            AlignedPair::new(Some(POANodePos(orig_node, 11)), None),
            AlignedPair::new(Some(POANodePos(orig_node, 12)), Some(9)),
            AlignedPair::new(Some(POANodePos(orig_node, 13)), Some(10)),
            AlignedPair::new(Some(POANodePos(orig_node, 14)), Some(11)),
            AlignedPair::new(Some(POANodePos(orig_node, 15)), Some(12)),
            AlignedPair::new(Some(POANodePos(orig_node, 16)), Some(13)),
            AlignedPair::new(Some(POANodePos(orig_node, 17)), Some(14)),
            AlignedPair::new(Some(POANodePos(orig_node, 18)), Some(15)),
            AlignedPair::new(Some(POANodePos(orig_node, 19)), Some(16)),
            AlignedPair::new(Some(POANodePos(orig_node, 20)), Some(17)),
            AlignedPair::new(Some(POANodePos(orig_node, 21)), Some(18)),
            AlignedPair::new(Some(POANodePos(orig_node, 22)), Some(19)),
        ];
        
        graph.add_aligned_sequence("seq2", seq2, &weights2, Some(&alignment)).unwrap();
        
        let seq_truth: Vec<&[u8]> = vec![
            b"#",
            b"GTCTGCTAT",
            b"ACT",
            b"GCGTACGTCGT",
            b"$",
        ];
        
        let w_truth: Vec<&[usize]> = vec![
            &[1],
            &[2, 2, 2, 2, 2, 2, 2, 2, 2],
            &[1, 1, 1],
            &[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            &[1],
        ];
        
        for (rank, n) in graph.toposorted.iter().enumerate() {
            let node_data = graph.node_data(*n);
            println!("{} {:?} seq: {:?}", rank, n, std::str::from_utf8(&node_data.sequence).unwrap());
            println!("{} {:?}   w: {:?}", rank, n, node_data.weights);
            println!("{} {:?} aln: {:?}", rank, n, node_data.aligned_intervals);
            
            assert_eq!(&node_data.sequence, seq_truth[rank]);
            assert_eq!(&node_data.weights, w_truth[rank]);
        }
    }
    
    #[test]
    fn test_from_msa() {
        tracing_subscriber::FmtSubscriber::builder()
            .with_writer(io::stderr)
            .with_file(false)
            .init();
        
        span!(Level::DEBUG, "test_from_msa");
        
        let mut file = File::open("../tests/test2_from_abpoa.truth.fa")
            .map(BufReader::new).unwrap();
        
        let graph = POASeqGraph::<u32>::try_from_fasta_msa(&mut file).unwrap();
        drop(file);
        
        let mut reader = File::open("../tests/test2_from_abpoa.fa")
            .map(BufReader::new)
            .map(fasta::Reader::new)
            .unwrap();
        
        for (seq_id, record) in reader.records().enumerate() {
            let r = record.unwrap();
            
            let seq_path = graph.path_for_sequence(seq_id);
            let seq = graph.path_sequence(seq_path).unwrap();
            
            assert_eq!(&seq, r.sequence().as_ref());
        }
        
        let path_seq1 = graph.path_for_sequence(0);
        let node_data = graph.node_data(path_seq1[4]);
        assert_eq!(node_data.aligned_intervals.len(), 1);
        assert_eq!(node_data.aligned_intervals[0].other_start(), 0);
        assert_eq!(graph.node_data(node_data.aligned_intervals[0].other()).sequence, b"C");
    }
}

