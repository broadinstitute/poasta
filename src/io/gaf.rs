use std::fmt::{self, Display, Write};

use rustc_hash::FxHashMap;
use itertools::Itertools;

use crate::aligner::Alignment;
use crate::graphs::poa::{POAGraph, POANodeIndex};
use crate::graphs::AlignableRefGraph;

use super::gfa::{Field, FieldValue};
use super::graph::GraphSegments;


pub struct GAFRecord {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub graph_path: String,
    pub path_length: usize,
    pub path_aln_start: usize,
    pub path_aln_end: usize,
    pub num_matches: usize,
    pub aln_block_len: usize,
    pub mapping_quality: usize,
    pub additional_fields: Vec<Field>,
}

impl Display for GAFRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let fields_str = self.additional_fields.iter()
            .fold(String::new(), |mut output, f| {
                let _ = write!(output, "\t{}", f);
                
                output
            });
        
        write!(
            f, 
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            self.graph_path,
            self.path_length,
            self.path_aln_start,
            self.path_aln_end,
            self.num_matches,
            self.aln_block_len,
            self.mapping_quality,
            fields_str.trim()
        )?;
        
        Ok(())
    }
}
            


pub fn alignment_to_gaf<Ix>(
    graph: &POAGraph<Ix>,
    graph_segments: &GraphSegments<Ix>,
    seq_name: &str,
    sequence: &[u8],
    alignment: &Alignment<POANodeIndex<Ix>>, 
    node_to_segment: &FxHashMap<POANodeIndex<Ix>, (usize, usize)>
) -> Option<GAFRecord>
where 
    Ix: petgraph::graph::IndexType + serde::de::DeserializeOwned,
{
    if alignment.is_empty() {
        return None;
    }
    
    let mut query_start = 0;
    let mut path_aln_start = 0;
    let mut path_segments = Vec::new();
    let mut cigar_ops = Vec::new();
    
    let mut at_aln_start = true;
    
    let mut last_match_segment_ix = 0;
    let mut last_match_segment_pos = 0;
    let mut num_matches = 0;
    for aln_pair in alignment.iter() {
        if at_aln_start {
            if aln_pair.is_insertion() {
                query_start += 1;
            } else if aln_pair.is_aligned() {
                let (segment_ix, segment_pos) = node_to_segment.get(&aln_pair.rpos.unwrap()).unwrap();
                path_aln_start = *segment_pos;
                path_segments.push(*segment_ix);
                cigar_ops.push(if graph.is_symbol_equal(aln_pair.rpos.unwrap(), sequence[aln_pair.qpos.unwrap()]) {
                    num_matches += 1;
                    '='
                } else {
                    'X'
                });
                    
                at_aln_start = false;
                last_match_segment_ix = path_segments.len() - 1;
                last_match_segment_pos = *segment_pos;
            }
        } else {
            match (aln_pair.rpos, aln_pair.qpos) {
                (Some(node), Some(qpos)) => {
                    let (segment_ix, segment_pos) = node_to_segment.get(&node).unwrap();
                    
                    if path_segments.last().copied() != Some(*segment_ix) {
                        path_segments.push(*segment_ix);
                    }
                    
                    cigar_ops.push(if graph.is_symbol_equal(node, sequence[qpos]) {
                        num_matches += 1;
                        '='
                    } else {
                        'X'
                    });
                    
                    last_match_segment_ix = path_segments.len() - 1;
                    last_match_segment_pos = *segment_pos;
                },
                
                (Some(node), None) => {
                    let (segment_ix, _) = node_to_segment.get(&node).unwrap();
                    
                    if path_segments.last().copied() != Some(*segment_ix) {
                        path_segments.push(*segment_ix);
                    }
                    
                    cigar_ops.push('D');
                },
                
                (None, Some(_)) => {
                    cigar_ops.push('I');
                },
                _ => unreachable!(),
            }
        }
    }
    
    let graph_path = path_segments[..=last_match_segment_ix].iter()
        .fold(String::new(), |mut output, segment_ix| {
            let _ = write!(output, ">{}", graph_segments.names[*segment_ix]);
            output
        });
    
    eprintln!("Path segments: {:?}", path_segments);
    eprintln!("Path segments (cut): {:?}", &path_segments[..=last_match_segment_ix]);
    let path_length = path_segments[..=last_match_segment_ix].iter()
        .map(|segment_ix| graph_segments.segment_lengths[*segment_ix])
        .sum();
    
    let path_aln_end = path_length - graph_segments.segment_lengths[last_match_segment_ix] + last_match_segment_pos;
    
    let query_end = alignment.iter().rev()
        .find(|aln_pair| aln_pair.is_aligned())
        .unwrap()
        .qpos.unwrap();
    
    let mut cigar_rle = cigar_ops.iter()
        .chunk_by(|op| **op)
        .into_iter()
        .map(|(op, group)| (op, group.count()))
        .collect::<Vec<_>>();
    
    match cigar_rle.last() {
        Some(('I', _)) => {
            let removed = cigar_rle.pop().unwrap();
            eprintln!("Removed insertion from end {}", removed.1);
        },
        Some(('D', _)) => {
            let removed = cigar_rle.pop().unwrap();
            eprintln!("Removed insertion from end {}", removed.1);
        }
        _ => (),
    }
    
    let aln_block_len = cigar_rle.iter()
        .map(|(_, count)| *count)
        .sum();
    
    let cigar_string = cigar_rle.iter()
        .fold(String::new(), |mut output, (op, count)| {
            let _ = write!(output, "{}{}", *count, *op);
            output
        });
    
    Some(GAFRecord {
        query_name: seq_name.to_string(),
        query_length: sequence.len(),
        query_start,
        query_end,
        strand: '+',
        graph_path,
        path_length,
        path_aln_start,
        path_aln_end,
        num_matches,
        aln_block_len,
        mapping_quality: 60,
        additional_fields: vec![
            Field { tag: "cg".to_string(), value: FieldValue::String(cigar_string) },
        ],
    })
}