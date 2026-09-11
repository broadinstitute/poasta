//! Write sequence-to-graph alignments as GAF (Graph Alignment Format).
//!
//! The input is a `Vec<AlignedPair<POANodeIndex<Ix>>>` produced by an
//! `AlignmentEngine` and a `GraphSegments<Ix>` produced by
//! `load_graph_from_gfa`. Leading and trailing indels are trimmed so the
//! emitted path covers only matched/mismatched positions.

use std::fmt::{self, Display, Write};

use itertools::Itertools;
use petgraph::graph::IndexType;
use rustc_hash::FxHashMap;

use crate::align::traits::AlignedPair;
use crate::graph::io::gfa::{Field, FieldValue, GraphSegments};
use crate::graph::poa::{POAGraph, POANodeIndex};

#[derive(Debug)]
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
        let mut fields_str = String::new();
        for field in &self.additional_fields {
            let _ = write!(fields_str, "\t{}", field);
        }

        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
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
            fields_str,
        )
    }
}

/// Convert an alignment produced by an `AlignmentEngine` into a GAF record.
/// Returns `None` if the alignment is empty or contains no aligned pairs
/// (i.e. no segment of the query actually matched the graph).
///
/// Behaviour vs. the stable lasagna port:
/// - Leading insertions *and* leading deletions are skipped during the
///   `at_aln_start` window (stable only skipped leading insertions, which
///   inflated `path_aln_start` to 0 when the optimum starts with a deletion).
/// - The trailing-indel trim no longer prints a stray "Removed insertion from
///   end" log; the deletion arm in stable also mislabeled itself "insertion".
pub fn alignment_to_gaf<Ix>(
    graph: &POAGraph<Ix>,
    graph_segments: &GraphSegments<Ix>,
    seq_name: &str,
    sequence: &[u8],
    alignment: &[AlignedPair<POANodeIndex<Ix>>],
    node_to_segment: &FxHashMap<POANodeIndex<Ix>, (usize, usize)>,
    score: u32,
) -> Option<GAFRecord>
where
    Ix: IndexType,
{
    if alignment.is_empty() {
        return None;
    }

    let mut query_start: usize = 0;
    let mut path_aln_start: usize = 0;
    let mut path_segments: Vec<usize> = Vec::new();
    let mut cigar_ops: Vec<char> = Vec::new();

    let mut at_aln_start = true;

    let mut last_match_segment_ix: usize = 0;
    let mut last_match_segment_pos: usize = 0;
    let mut num_matches: usize = 0;

    for aln_pair in alignment.iter() {
        debug_assert!(
            aln_pair.is_aligned() || aln_pair.is_insertion() || aln_pair.is_deletion(),
            "AlignedPair must be aligned, insertion, or deletion"
        );

        if at_aln_start {
            // Skip both leading insertions and leading deletions. Leading
            // insertions advance query_start; leading deletions are dropped
            // silently — they would otherwise inflate the path with segments
            // the alignment doesn't actually cover.
            if aln_pair.is_insertion() {
                query_start += 1;
                continue;
            }
            if aln_pair.is_deletion() {
                continue;
            }

            // First aligned pair anchors the path.
            let node = aln_pair.node().unwrap();
            let qpos = aln_pair.query_pos().unwrap();
            let (segment_ix, segment_pos) = *node_to_segment
                .get(&node)
                .expect("aligned graph node missing from node_to_segment");
            path_aln_start = segment_pos;
            path_segments.push(segment_ix);
            cigar_ops.push(if graph.node_symbol(node) == sequence[qpos] {
                num_matches += 1;
                '='
            } else {
                'X'
            });

            at_aln_start = false;
            last_match_segment_ix = path_segments.len() - 1;
            last_match_segment_pos = segment_pos;
            continue;
        }

        match (aln_pair.node(), aln_pair.query_pos()) {
            (Some(node), Some(qpos)) => {
                let (segment_ix, segment_pos) = *node_to_segment
                    .get(&node)
                    .expect("aligned graph node missing from node_to_segment");

                if path_segments.last().copied() != Some(segment_ix) {
                    path_segments.push(segment_ix);
                }

                cigar_ops.push(if graph.node_symbol(node) == sequence[qpos] {
                    num_matches += 1;
                    '='
                } else {
                    'X'
                });

                last_match_segment_ix = path_segments.len() - 1;
                last_match_segment_pos = segment_pos;
            }
            (Some(node), None) => {
                let (segment_ix, _) = *node_to_segment
                    .get(&node)
                    .expect("aligned graph node missing from node_to_segment");

                if path_segments.last().copied() != Some(segment_ix) {
                    path_segments.push(segment_ix);
                }

                cigar_ops.push('D');
            }
            (None, Some(_)) => {
                cigar_ops.push('I');
            }
            (None, None) => unreachable!("debug_assert guards this branch"),
        }
    }

    if path_segments.is_empty() {
        // No aligned pair was ever seen — every operation was an indel. Don't
        // emit a record; the caller can decide whether to log/skip.
        return None;
    }

    let mut graph_path = String::new();
    for &segment_ix in &path_segments[..=last_match_segment_ix] {
        let _ = write!(graph_path, ">{}", graph_segments.names[segment_ix]);
    }

    let path_length: usize = path_segments[..=last_match_segment_ix]
        .iter()
        .map(|&segment_ix| graph_segments.segment_lengths[segment_ix])
        .sum();

    let path_aln_end =
        path_length - graph_segments.segment_lengths[path_segments[last_match_segment_ix]]
            + last_match_segment_pos;

    let query_end = alignment
        .iter()
        .rev()
        .find(|p| p.is_aligned())
        .and_then(|p| p.query_pos())
        .expect("at least one aligned pair (path_segments non-empty)");

    // Run-length encode CIGAR ops, then trim a trailing I or D.
    let mut cigar_rle: Vec<(char, usize)> = cigar_ops
        .iter()
        .chunk_by(|op| **op)
        .into_iter()
        .map(|(op, group)| (op, group.count()))
        .collect();

    match cigar_rle.last() {
        Some(('I', n)) => {
            tracing::trace!(query = %seq_name, n = n, "trimmed trailing insertion");
            cigar_rle.pop();
        }
        Some(('D', n)) => {
            tracing::trace!(query = %seq_name, n = n, "trimmed trailing deletion");
            cigar_rle.pop();
        }
        _ => (),
    }

    let aln_block_len: usize = cigar_rle.iter().map(|(_, count)| *count).sum();

    let mut cigar_string = String::new();
    for (op, count) in &cigar_rle {
        let _ = write!(cigar_string, "{}{}", count, op);
    }

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
            Field {
                tag: "AS".to_string(),
                value: FieldValue::Integer(score as i64),
            },
            Field {
                tag: "cg".to_string(),
                value: FieldValue::String(cigar_string),
            },
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::io::gfa::load_graph_from_gfa;
    use crate::graph::traits::GraphBase;

    /// Build `node_to_segment` from `GraphSegments` by walking each segment
    /// using its recorded length (NOT `successors().next()` until end_node,
    /// which can run off the end of a tip-attached segment).
    fn build_node_to_segment<Ix: IndexType>(
        graph: &POAGraph<Ix>,
        segments: &GraphSegments<Ix>,
    ) -> FxHashMap<POANodeIndex<Ix>, (usize, usize)> {
        let mut map = FxHashMap::default();
        for (segment_ix, _) in segments.names.iter().enumerate() {
            let start_node = segments.start_nodes[segment_ix];
            let end_node = segments.end_nodes[segment_ix];
            let len = segments.segment_lengths[segment_ix];

            map.insert(start_node, (segment_ix, 0));
            if len == 1 {
                continue;
            }

            let mut curr = start_node;
            for pos in 1..len {
                // Within a GFA segment, internal nodes have exactly one
                // successor (the next character in the segment).
                let succ = graph
                    .successors(curr)
                    .find(|&s| {
                        // Filter out edges to tip-attached end-sentinel if any.
                        // For internal positions, the segment continues; at
                        // pos == len-1 we want end_node.
                        if pos < len - 1 {
                            s != end_node || s == start_node
                        } else {
                            s == end_node
                        }
                    })
                    .expect("segment walk found no internal successor");
                map.insert(succ, (segment_ix, pos));
                curr = succ;
            }
        }
        map
    }

    #[test]
    fn test_alignment_to_gaf_perfect_match() {
        let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test.gfa");
        let loaded = load_graph_from_gfa::<u32>(&p).expect("load test.gfa");
        let node_to_segment = build_node_to_segment(&loaded.graph, &loaded.graph_segments);

        // Build a synthetic alignment that perfectly matches segment "1"
        // (AAACTGC, len 7) — exactly all-`=`, no indels.
        let mut alignment: Vec<AlignedPair<POANodeIndex<u32>>> = Vec::new();
        let segment_1_start = loaded.graph_segments.start_nodes[0];
        let mut curr = segment_1_start;
        let segment_1_end = loaded.graph_segments.end_nodes[0];
        for q in 0..7 {
            alignment.push(AlignedPair::new(Some(curr), Some(q)));
            if curr == segment_1_end {
                break;
            }
            curr = loaded
                .graph
                .successors(curr)
                .next()
                .expect("walk segment 1");
        }

        let seq = b"AAACTGC";
        let rec = alignment_to_gaf(
            &loaded.graph,
            &loaded.graph_segments,
            "q1",
            seq,
            &alignment,
            &node_to_segment,
            0,
        )
        .expect("non-empty record");

        assert_eq!(rec.query_name, "q1");
        assert_eq!(rec.query_length, 7);
        assert_eq!(rec.query_start, 0);
        assert_eq!(rec.query_end, 6);
        assert_eq!(rec.num_matches, 7);
        assert_eq!(rec.aln_block_len, 7);
        assert_eq!(rec.graph_path, ">1");
        // CIGAR field should be `7=`.
        let cg = rec
            .additional_fields
            .iter()
            .find(|f| f.tag == "cg")
            .expect("cg field");
        assert!(matches!(&cg.value, FieldValue::String(s) if s == "7="));
        // AS:i:0 for a perfect match.
        let as_tag = rec
            .additional_fields
            .iter()
            .find(|f| f.tag == "AS")
            .expect("AS field");
        assert!(matches!(&as_tag.value, FieldValue::Integer(0)));
    }

    #[test]
    fn test_alignment_to_gaf_skips_leading_and_trailing_indels() {
        let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test.gfa");
        let loaded = load_graph_from_gfa::<u32>(&p).expect("load test.gfa");
        let node_to_segment = build_node_to_segment(&loaded.graph, &loaded.graph_segments);

        let segment_1_start = loaded.graph_segments.start_nodes[0];
        let mut alignment: Vec<AlignedPair<POANodeIndex<u32>>> = Vec::new();
        // Leading insertion (query base before any graph match)
        alignment.push(AlignedPair::new(None, Some(0)));
        // A single aligned pair at segment 1's first node, query position 1
        alignment.push(AlignedPair::new(Some(segment_1_start), Some(1)));
        // Trailing deletion: graph advances without consuming query
        let segment_1_next = loaded
            .graph
            .successors(segment_1_start)
            .next()
            .expect("walk segment 1");
        alignment.push(AlignedPair::new(Some(segment_1_next), None));

        let seq = b"XA";
        let rec = alignment_to_gaf(
            &loaded.graph,
            &loaded.graph_segments,
            "q",
            seq,
            &alignment,
            &node_to_segment,
            42,
        )
        .expect("non-empty record");

        assert_eq!(rec.query_start, 1, "leading insertion should advance query_start");
        assert_eq!(
            rec.aln_block_len, 1,
            "trailing deletion should be trimmed from aln_block_len"
        );
        // Score round-trips into AS tag.
        let as_tag = rec
            .additional_fields
            .iter()
            .find(|f| f.tag == "AS")
            .expect("AS field");
        assert!(matches!(&as_tag.value, FieldValue::Integer(42)));
    }
}
