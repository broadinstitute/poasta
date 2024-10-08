//! Contains types for representing sequence-to-graph alignments.
use std::ops::Range;

use petgraph::graph::IndexType;

use tracing::debug;

use crate::aligner::astar::AlignableGraphNodePos;
use crate::aligner::utils::AlignedPair;

use super::poa::{POANodeIndex, POASeqGraph, SplitTracker};

/// Refers to a specific position with a POA graph node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct POANodePos<Ix>(pub POANodeIndex<Ix>, pub usize)
where
    Ix: IndexType;

impl<Ix> POANodePos<Ix>
where 
    Ix: IndexType,
{
    #[inline]
    pub fn node_equal(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<Ix> AlignableGraphNodePos for POANodePos<Ix>
where 
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    
    #[inline]
    fn new(node: Self::NodeType, pos: usize) -> Self {
        POANodePos(node, pos)
    }
    
    #[inline]
    fn node(&self) -> Self::NodeType {
        self.0
    }
    
    #[inline]
    fn pos(&self) -> usize {
        self.1
    }
}


/// Represents the aligned state of a pair of positions in a sequence-to-graph alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentClassification<Ix>
where 
    Ix: IndexType
{
    Start,
    Match,
    Mismatch(Option<POANodePos<Ix>>),
    Insertion,
    Deletion,
}


/// Divide a sequence-to-graph alignment into multiple "blocks".
/// 
/// Blocks are defined as a sequence of aligned pairs that are classified as either:
/// - All matches
/// - All mismatches, where the mismatching sequence is not present in other aligned nodes
/// - All mismatches, where the mismatching sequence is present in another aligned node, and this node is the same for a block
/// - All deletions
/// - All insertions
pub(crate) struct AlignmentBlocks<'a, Ix>
where
    Ix: IndexType
{
    /// The alignment to process
    alignment: &'a [AlignedPair<POANodePos<Ix>>],
    
    /// Iterator over the alignment
    aln_iter: std::iter::Enumerate<std::slice::Iter<'a, AlignedPair<POANodePos<Ix>>>>,
    
    // Previous and current alignment block type
    prev_state: AlignmentClassification<Ix>,
    curr_state: AlignmentClassification<Ix>,
    
    /// The intervals in this vector represent indices in the alignment marking ranges spanning 
    /// a single node in the POA graph.
    ///
    /// The ranges are denoted as half-open [start, end) intervals, where `start` is inclusive and `end` is exclusive.
    ///
    /// A visual example, where the block spans a range in the alignment:
    ///
    /// ```text
    /// |-----------------------------|   <- alignment block (e.g., all matches)
    /// |----------|------------------|   <- example node intervals.
    ///     N1              N2            <- POA graph nodes
    /// ```
    aln_node_ivals: Vec<(usize, usize)>,
    
    /// The start index of the current node interval in the alignment block
    node_ival_start: usize,
    
    /// The previous alignment position was aligned to this POA graph node and position
    prev_node: POANodePos<Ix>,
}


impl<'a, Ix> AlignmentBlocks<'a, Ix>
where 
    Ix: IndexType
{
    pub fn new(alignment: &'a [AlignedPair<POANodePos<Ix>>], start_node: POANodeIndex<Ix>) -> Self {
        AlignmentBlocks {
            alignment,
            aln_iter: alignment.iter().enumerate(),
            prev_state: AlignmentClassification::Start,
            curr_state: AlignmentClassification::Start,
            aln_node_ivals: Vec::new(),
            node_ival_start: 0,
            prev_node: POANodePos(start_node, 0),
        }
    }
    
    fn start_new_block(&mut self, curr_pos: usize) {
        self.aln_node_ivals.clear();
        self.node_ival_start = curr_pos;
    }
    
    pub fn next_block(&mut self, graph: &POASeqGraph<Ix>, seq: &[u8], node_splits: &SplitTracker<Ix>) -> Option<AlignmentBlock<Ix>> {
        let mut to_return = None;
        
        while let Some((aln_pos, pair)) = self.aln_iter.next() {
            // Update any node references in the alignment that refer to old, split, nodes
            let pair = AlignedPair::new(
                pair.node_pos().map(|pos| node_splits.to_new_node_pos(pos)),
                pair.query_pos(),
            );
            self.prev_node = node_splits.to_new_node_pos(self.prev_node);
            
            self.prev_state = match self.curr_state {
                AlignmentClassification::Mismatch(Some(pos)) => {
                    AlignmentClassification::Mismatch(Some(node_splits.to_new_node_pos(pos)))
                },
                _ => self.curr_state
            };
            
            self.curr_state = self.classify_aln_pair(&pair, graph, seq);
            
            debug!("State: {:?} -> {:?}", self.prev_state, self.curr_state);
            
            match (self.prev_state, self.curr_state) {
                (AlignmentClassification::Start, _) => {
                    self.start_new_block(aln_pos);
                },
                
                // Both matches and deletions can traverse multiple nodes in the graph.
                // If traversing to a new node, update node_ivals accordingly.
                (AlignmentClassification::Match, AlignmentClassification::Match) 
                | (AlignmentClassification::Deletion, AlignmentClassification::Deletion) => {
                    if let Some(n) = pair.node() {
                        if n != self.prev_node.node() {
                            self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                            self.node_ival_start = aln_pos;
                        }
                    }
                },
                
                // Mismatches can also traverse multiple nodes in the graph, but we first check
                // if the mismatching residue is present in another aligned node.
                (AlignmentClassification::Mismatch(other_node_aln1), AlignmentClassification::Mismatch(other_node_aln2)) => {
                    if other_node_aln1.map(|v| v.node()) != other_node_aln2.map(|v| v.node()) {
                        // Mismatch, but the mismatching residue is present in another aligned node that is 
                        // different from the previous. End the current block and start a new one.
                        self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                        self.node_ival_start = aln_pos;
                        
                        if let Some(other_node) = other_node_aln1 {
                            to_return = Some(self.make_block_mismatches_present(other_node));
                        } else {
                            to_return = Some(self.make_block_mismatches_new(node_splits));
                        }
                        
                        self.start_new_block(aln_pos);
                    } else if let Some(node) = pair.node() {
                        // Else, if the mismatch remains in the same node (or is not present in the graph), 
                        // just check if we need to update the node intervals.
                        if node != self.prev_node.node() {
                            self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                            self.node_ival_start = aln_pos;
                        }
                    }
                },
                
                // For insertions we don't have to do anything (alignment doesn't traverse additional graph nodes)
                (AlignmentClassification::Insertion, AlignmentClassification::Insertion) => (),
                
                // A block of matches ends
                (AlignmentClassification::Match, AlignmentClassification::Mismatch(_))
                | (AlignmentClassification::Match, AlignmentClassification::Deletion)
                | (AlignmentClassification::Match, AlignmentClassification::Insertion)
                => {
                    self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                    to_return = Some(self.make_block_matches(node_splits));
                    
                    self.start_new_block(aln_pos);
                },
                
                // A block of mismatches ends
                (AlignmentClassification::Mismatch(aligned_to), AlignmentClassification::Match)
                | (AlignmentClassification::Mismatch(aligned_to), AlignmentClassification::Deletion)
                | (AlignmentClassification::Mismatch(aligned_to), AlignmentClassification::Insertion)
                => {
                    self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                    
                    if let Some(other_node) = aligned_to {
                        to_return = Some(self.make_block_mismatches_present(other_node));
                    } else {
                        to_return = Some(self.make_block_mismatches_new(node_splits));
                    }
                    
                    self.start_new_block(aln_pos);
                },
                
                // A block of insertions ends
                (AlignmentClassification::Insertion, AlignmentClassification::Match)
                | (AlignmentClassification::Insertion, AlignmentClassification::Mismatch(_))
                => {
                    self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                    to_return = Some(self.make_block_insertion());
                    
                    self.start_new_block(aln_pos);
                },
                
                (AlignmentClassification::Insertion, AlignmentClassification::Deletion) => 
                    panic!("Invalid alignment, deletion immediately after insertion."),
                    
                // A block of deletions ends
                (AlignmentClassification::Deletion, AlignmentClassification::Match)
                | (AlignmentClassification::Deletion, AlignmentClassification::Mismatch(_))
                => {
                    self.aln_node_ivals.push((self.node_ival_start, aln_pos));
                    to_return = Some(self.make_block_deletion(node_splits));
                    
                    self.start_new_block(aln_pos);
                },
                
                (AlignmentClassification::Deletion, AlignmentClassification::Insertion) => 
                    panic!("Invalid alignment, insertion immediately after deletion."),
                
                (_, AlignmentClassification::Start) => 
                    panic!("Invalid state transition to AlignmentClassification::Start.")
                
            }
            
            if let Some(n) = pair.node_pos() {
                self.prev_node = n;
            }
            
            if to_return.is_some() {
                break;
            }
        }
        
        if to_return.is_none() && self.node_ival_start < self.alignment.len() {
            self.aln_node_ivals.push((self.node_ival_start, self.alignment.len()));
            
            to_return = match self.curr_state {
                AlignmentClassification::Match => Some(self.make_block_matches(node_splits)),
                AlignmentClassification::Mismatch(aligned_to) => {
                    if let Some(other_node) = aligned_to {
                        Some(self.make_block_mismatches_present(other_node))
                    } else {
                        Some(self.make_block_mismatches_new(node_splits))
                    }
                },
                AlignmentClassification::Deletion => Some(self.make_block_deletion(node_splits)),
                AlignmentClassification::Insertion => Some(self.make_block_insertion()),
                _ => None,
            };
            
            self.node_ival_start = self.alignment.len();
        }
        
        to_return
    }
    
    /// Classifies a specific pair of aligned residues in the sequence-to-graph alignment.
    /// 
    /// The aligned pair can either be a match, mismatch, insertion or deletion. Mismatches are further
    /// classified as whether the mismatching residue is present in another aligned node.
    fn classify_aln_pair(&self, pair: &AlignedPair<POANodePos<Ix>>, graph: &POASeqGraph<Ix>, seq: &[u8]) -> AlignmentClassification<Ix> {
        let node = pair.node();
        let query_pos = pair.query_pos();
        
        match (node, query_pos) {
            (Some(n), Some(q)) => {
                let node_pos = pair.within_node_pos().unwrap();
                let node_symbol = graph.get_node_symbol(n, node_pos);
                let query_symbol = seq[q];
                
                if node_symbol == query_symbol {
                    AlignmentClassification::Match
                } else {
                    // Check if mismatching symbol is present in another aligned node
                    let mut match_in_other_node = None;
                    for aligned_ival in &graph.node_data(n).aligned_intervals {
                        if let Some(other_pos) = aligned_ival.pos_to_other(node_pos) {
                            let other_node = aligned_ival.other();
                            let other_node_data = graph.node_data(other_node);
                            
                            if other_node_data.sequence[other_pos] == seq[q] {
                                match_in_other_node = Some(POANodePos(other_node, other_pos));
                                break;
                            }
                        }
                    }
                    
                    AlignmentClassification::Mismatch(match_in_other_node)
                }
            },
            (Some(_), None) => AlignmentClassification::Deletion,
            (None, Some(_)) => AlignmentClassification::Insertion,
            (None, None) => panic!("Invalid alignment pair: {:?}", pair),
        }
    }
    
    fn get_query_range(&self, aln_start: usize, aln_end: usize) -> Result<Range<usize>, ()> {
        let Some(qstart) = self.alignment[aln_start].query_pos() else {
            return Err(());
        };
        
        if aln_end > self.alignment.len() {
            return Err(());
        }
        
        let ix = aln_end - 1;
        
        let qry_end = if let Some(qend) = self.alignment[ix].query_pos() {
            qend + 1
        } else {
            return Err(());
        };
        
        Ok(qstart..qry_end)
    }
    
    fn aln_ivals_to_node_ivals(&self, node_splits: &SplitTracker<Ix>) -> Vec<(POANodePos<Ix>, usize)> {
        self.aln_node_ivals.iter()
            .map(|(start, end)| {
                let ival_length = *end - *start;
                let node = self.alignment[*start].node().unwrap();
                let pos = self.alignment[*start].within_node_pos().unwrap();
                let node_pos = node_splits.to_new_node_pos(POANodePos(node, pos));
                
                (node_pos, ival_length)
            })
            .collect()
    }
    
    fn make_block_matches(&self, node_splits: &SplitTracker<Ix>) -> AlignmentBlock<Ix> {
        let first_ival = self.aln_node_ivals.first().unwrap();
        let last_ival = self.aln_node_ivals.last().unwrap();
        
        let aln_start = first_ival.0;
        let aln_end = last_ival.1;
        
        let qry_range = self.get_query_range(aln_start, aln_end).unwrap();
        
        AlignmentBlock {
            aln_range: aln_start..aln_end,
            block_type: AlignmentBlockType::Match {
                qry_range,
                node_ivals: self.aln_ivals_to_node_ivals(node_splits),
            },
        }
    }
    
    fn make_block_mismatches_new(&self, node_splits: &SplitTracker<Ix>) -> AlignmentBlock<Ix> {
        let first_ival = self.aln_node_ivals.first().unwrap();
        let last_ival = self.aln_node_ivals.last().unwrap();
        
        let aln_start = first_ival.0;
        let aln_end = last_ival.1;
        
        let qry_range = self.get_query_range(aln_start, aln_end).unwrap();
        
        AlignmentBlock {
            aln_range: aln_start..aln_end,
            block_type: AlignmentBlockType::MismatchNew {
                qry_range,
                node_ivals: self.aln_ivals_to_node_ivals(node_splits),
            },
        }
    }
    
    fn make_block_mismatches_present(&self, other_node: POANodePos<Ix>) -> AlignmentBlock<Ix> {
        let first_ival = self.aln_node_ivals.first().unwrap();
        let last_ival = self.aln_node_ivals.last().unwrap();
        
        let aln_start = first_ival.0;
        let aln_end = last_ival.1;
        
        let qry_range = self.get_query_range(aln_start, aln_end).unwrap();
        
        AlignmentBlock {
            aln_range: aln_start..aln_end,
            block_type: AlignmentBlockType::MismatchPresent {
                qry_range,
                other_node,
            },
        }
    }
    
    fn make_block_insertion(&self) -> AlignmentBlock<Ix> {
        let first_ival = self.aln_node_ivals.first().unwrap();
        let last_ival = self.aln_node_ivals.last().unwrap();
        
        let aln_start = first_ival.0;
        let aln_end = last_ival.1;
        
        let qry_range = self.get_query_range(aln_start, aln_end).unwrap();
        
        AlignmentBlock {
            aln_range: aln_start..aln_end,
            block_type: AlignmentBlockType::Insertion {
                qry_range
            },
        }
    }
    
    fn make_block_deletion(&self, node_splits: &SplitTracker<Ix>) -> AlignmentBlock<Ix> {
        let first_ival = self.aln_node_ivals.first().unwrap();
        let last_ival = self.aln_node_ivals.last().unwrap();
        
        let aln_start = first_ival.0;
        let aln_end = last_ival.1;
        
        AlignmentBlock {
            aln_range: aln_start..aln_end,
            block_type: AlignmentBlockType::Deletion {
                node_ivals: self.aln_ivals_to_node_ivals(node_splits),
            },
        }
    }
}


#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlignmentBlockType<Ix>
where 
    Ix: IndexType
{
    Match { 
        qry_range: Range<usize>,
        
        /// The alignment intervals of each traversed node in the graph. Encoded as `(POANodePos, length)`, 
        /// representing the start position and the length of an interval.
        node_ivals: Vec<(POANodePos<Ix>, usize)>
    },
    
    /// Mismatching sequence not present in other aligned nodes, so requires a new graph node
    MismatchNew { 
        qry_range: Range<usize>,
        
        /// The alignment intervals of each traversed node in the graph. Encoded as `(POANodePos, length)`, 
        /// representing the start position and the length of an interval.
        node_ivals: Vec<(POANodePos<Ix>, usize)>
    },
    
    /// Mismatching sequence present in another aligned node
    MismatchPresent {
        qry_range: Range<usize>,
        
        /// The node and position of the other aligned node that contains the mismatching sequence
        other_node: POANodePos<Ix>,
    },
    
    /// Block of newly inserted sequence
    Insertion {
        qry_range: Range<usize>,
    },
    
    /// Block of deleted sequence
    Deletion {
        /// The alignment intervals of each traversed node in the graph. Encoded as `(POANodePos, length)`, 
        /// representing the start position and the length of an interval.
        node_ivals: Vec<(POANodePos<Ix>, usize)>
    }
    
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentBlock<Ix>
where 
    Ix: IndexType,
{
    /// The range spanning this block in the alignment string
    pub aln_range: Range<usize>,
    
    /// Type of alignment block, with additional data
    pub block_type: AlignmentBlockType<Ix>,
}


#[cfg(test)]
mod tests {
    use crate::graph::poa::POASeqGraph;
    use super::{AlignedPair, AlignmentBlocks, AlignmentBlock, AlignmentBlockType, POANodePos, SplitTracker};

    #[test]
    fn test_alignment_block_mismatches() {
        let mut test_graph = POASeqGraph::<usize>::new();
        let node = test_graph.add_node(b"CCGCTTTTCGCG".to_vec());
        
        let qry_seq = b"CCGCAAAACGCG".to_vec();
        let aln: Vec<AlignedPair<POANodePos<usize>>> = vec![
            AlignedPair::new(Some(POANodePos(node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(node, 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(node, 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(node, 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(node, 7)), Some(7)),
            AlignedPair::new(Some(POANodePos(node, 8)), Some(8)),
            AlignedPair::new(Some(POANodePos(node, 9)), Some(9)),
            AlignedPair::new(Some(POANodePos(node, 10)), Some(10)),
            AlignedPair::new(Some(POANodePos(node, 11)), Some(11)),
        ];
        
        let mut blocks = AlignmentBlocks::new(&aln, test_graph.start_node);
        
        let mut all_blocks = Vec::default();
        let nodes_split = SplitTracker::default();
        while let Some(block) = blocks.next_block(&test_graph, &qry_seq, &nodes_split) {
            all_blocks.push(block);
        }
        
        assert_eq!(all_blocks, vec![
            AlignmentBlock {
                aln_range: 0..4,
                block_type: AlignmentBlockType::Match {
                    qry_range: 0..4,
                    node_ivals: vec![(POANodePos(node, 0), 4)]
                }
            },
            AlignmentBlock {
                aln_range: 4..8,
                block_type: AlignmentBlockType::MismatchNew {
                    qry_range: 4..8,
                    node_ivals: vec![(POANodePos(node, 4), 4)]
                }
            },
            AlignmentBlock {
                aln_range: 8..12,
                block_type: AlignmentBlockType::Match {
                    qry_range: 8..12,
                    node_ivals: vec![(POANodePos(node, 8), 4)]
                }
            },
        ]);
    }
    
    #[test]
    fn test_alignment_block_insertion() {
        let mut test_graph = POASeqGraph::<usize>::new();
        let node = test_graph.add_node(b"CCGCTTTTCGCG".to_vec());
        
        let qry_seq = b"CCGCAAAATTTCGCG".to_vec();
        let aln: Vec<AlignedPair<POANodePos<usize>>> = vec![
            AlignedPair::new(Some(POANodePos(node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(node, 4)), Some(4)),
            AlignedPair::new(Some(POANodePos(node, 5)), Some(5)),
            AlignedPair::new(Some(POANodePos(node, 6)), Some(6)),
            AlignedPair::new(Some(POANodePos(node, 7)), Some(7)),
            AlignedPair::new(None, Some(8)),
            AlignedPair::new(None, Some(9)),
            AlignedPair::new(None, Some(10)),
            AlignedPair::new(Some(POANodePos(node, 8)), Some(11)),
            AlignedPair::new(Some(POANodePos(node, 9)), Some(12)),
            AlignedPair::new(Some(POANodePos(node, 10)), Some(13)),
            AlignedPair::new(Some(POANodePos(node, 11)), Some(14)),
        ];
        
        let mut blocks = AlignmentBlocks::new(&aln, test_graph.start_node);
        
        let mut all_blocks = Vec::default();
        let nodes_split = SplitTracker::default();
        while let Some(block) = blocks.next_block(&test_graph, &qry_seq, &nodes_split) {
            all_blocks.push(block);
        }
        
        assert_eq!(all_blocks, vec![
            AlignmentBlock {
                aln_range: 0..4,
                block_type: AlignmentBlockType::Match {
                    qry_range: 0..4,
                    node_ivals: vec![(POANodePos(node, 0), 4)]
                }
            },
            AlignmentBlock {
                aln_range: 4..8,
                block_type: AlignmentBlockType::MismatchNew {
                    qry_range: 4..8,
                    node_ivals: vec![(POANodePos(node, 4), 4)]
                }
            },
            AlignmentBlock {
                aln_range: 8..11,
                block_type: AlignmentBlockType::Insertion {
                    qry_range: 8..11,
                }
            },
            AlignmentBlock {
                aln_range: 11..15,
                block_type: AlignmentBlockType::Match {
                    qry_range: 11..15,
                    node_ivals: vec![(POANodePos(node, 8), 4)]
                }
            },
        ]);
    }
    
    #[test]
    fn test_alignment_block_deletion() {
        let mut test_graph = POASeqGraph::<usize>::new();
        let node = test_graph.add_node(b"CCGCTTTTCGCG".to_vec());
        
        let qry_seq = b"CCGCCGCG".to_vec();
        let aln: Vec<AlignedPair<POANodePos<usize>>> = vec![
            AlignedPair::new(Some(POANodePos(node, 0)), Some(0)),
            AlignedPair::new(Some(POANodePos(node, 1)), Some(1)),
            AlignedPair::new(Some(POANodePos(node, 2)), Some(2)),
            AlignedPair::new(Some(POANodePos(node, 3)), Some(3)),
            AlignedPair::new(Some(POANodePos(node, 4)), None),
            AlignedPair::new(Some(POANodePos(node, 5)), None),
            AlignedPair::new(Some(POANodePos(node, 6)), None),
            AlignedPair::new(Some(POANodePos(node, 7)), None),
            AlignedPair::new(Some(POANodePos(node, 8)), Some(4)),
            AlignedPair::new(Some(POANodePos(node, 9)), Some(5)),
            AlignedPair::new(Some(POANodePos(node, 10)), Some(6)),
            AlignedPair::new(Some(POANodePos(node, 11)), Some(7)),
        ];
        
        let mut blocks = AlignmentBlocks::new(&aln, test_graph.start_node);
        
        let mut all_blocks = Vec::default();
        let nodes_split = SplitTracker::default();
        while let Some(block) = blocks.next_block(&test_graph, &qry_seq, &nodes_split) {
            all_blocks.push(block);
        }
        
        assert_eq!(all_blocks, vec![
            AlignmentBlock {
                aln_range: 0..4,
                block_type: AlignmentBlockType::Match {
                    qry_range: 0..4,
                    node_ivals: vec![(POANodePos(node, 0), 4)]
                }
            },
            AlignmentBlock {
                aln_range: 4..8,
                block_type: AlignmentBlockType::Deletion {
                    node_ivals: vec![(POANodePos(node, 4), 4)]
                }
            },
            AlignmentBlock {
                aln_range: 8..12,
                block_type: AlignmentBlockType::Match {
                    qry_range: 4..8,
                    node_ivals: vec![(POANodePos(node, 8), 4)]
                }
            },
        ]);
    }
}

