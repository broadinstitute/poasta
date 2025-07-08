//! Graph serialization to disk using serde

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Read, Write, BufRead};
use std::path::Path;
use std::{char, fmt};

use flate2::read::MultiGzDecoder;
use noodles::fasta;
use petgraph::dot::Dot;
use petgraph::graph::IndexType;
use petgraph::visit::EdgeRef;
use petgraph::visit::IntoEdgeReferences;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::de::DeserializeOwned;

use crate::errors::PoastaError;
use crate::graphs::poa::{POAGraph, POAGraphWithIx, POANodeData, POANodeIndex, Sequence};
use crate::graphs::AlignableRefGraph;

use super::gfa::{GfaLine, Strand};

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
    let is_gzipped = p
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
    let reader_inner: Box<dyn std::io::BufRead> = if is_gzipped {
        Box::new(File::open(p).map(MultiGzDecoder::new).map(BufReader::new)?)
    } else {
        Box::new(File::open(p).map(BufReader::new)?)
    };
    let mut reader = fasta::io::Reader::new(reader_inner);

    let mut graph = POAGraph::<u32>::new();
    let mut nodes_per_col: Vec<Vec<POANodeIndex<u32>>> = Vec::new();
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

            let node_ix = nodes_per_col[col]
                .iter()
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
                graph
                    .sequences
                    .push(Sequence(std::str::from_utf8(seq.name()).unwrap().to_string(), node_ix));
            }

            prev_node = Some(node_ix);
        }
    }

    graph.post_process()?;

    Ok(POAGraphWithIx::U32(graph))
}

pub struct POAGraphFromGFA<Ix>
where 
    Ix: petgraph::graph::IndexType,
{
    pub graph: POAGraph<Ix>,
    pub graph_segments: GraphSegments<Ix>,
}

#[derive(Debug, Default)]
pub struct GraphSegments<Ix> 
where 
    Ix: petgraph::graph::IndexType,
{
    pub names: Vec<String>,
    pub start_nodes: Vec<POANodeIndex<Ix>>,
    pub end_nodes: Vec<POANodeIndex<Ix>>,
    pub segment_lengths: Vec<usize>,
}


pub fn load_graph_from_gfa<Ix>(path: impl AsRef<Path>) -> Result<POAGraphFromGFA<Ix>, PoastaError> 
where
    Ix: petgraph::graph::IndexType + serde::de::DeserializeOwned,
{
    let p = path.as_ref();
    
    let is_gzipped = p
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
    let reader_inner: BufReader<Box<dyn std::io::Read>> = BufReader::new(if is_gzipped {
        Box::new(File::open(p).map(MultiGzDecoder::new)?) as Box<dyn std::io::Read>
    } else {
        Box::new(File::open(p)?) as Box<dyn std::io::Read>
    });
    
    let mut graph = POAGraph::new();
    let mut graph_segments = GraphSegments::default();
    let mut name_to_ix = FxHashMap::default();
    let mut links_to_add = Vec::new();
    
    for line in reader_inner.lines() {
        let line = line?;
        let trimmed = line.trim();
        
        if trimmed.is_empty() {
            continue;
        }
        
        match GfaLine::try_from(trimmed) {
            Ok(gfa_line) => {
                match gfa_line {
                    GfaLine::Segment(segment) => {
                        if let Some(seq) = segment.sequence {
                            let weights = vec![1; seq.len()];
                            let (start, end) = graph.add_nodes_for_sequence(seq.as_bytes(), &weights, 0, seq.len()).unwrap();                        
                            
                            name_to_ix.insert(segment.sid.clone(), graph_segments.names.len());
                            graph_segments.names.push(segment.sid.clone());
                            graph_segments.start_nodes.push(start);
                            graph_segments.end_nodes.push(end);
                            graph_segments.segment_lengths.push(seq.len());
                        } else {
                            eprintln!("Omitting segment {:?} because it has no sequence.", segment.sid);
                        }
                    },
                    GfaLine::Link(link) => {
                        if link.strand1 == Strand::Reverse || link.strand2 == Strand::Reverse {
                            eprintln!("Links using the reverse strand of a segment are not supported!");
                            eprintln!("Link: {link:?}");
                            return Err(PoastaError::GraphError);
                        }
                        
                        let from_ix = name_to_ix.get(&link.sid1);
                        let to_ix = name_to_ix.get(&link.sid2);
                        
                        // No corresponding segment yet, maybe later after parsing the entire file?
                        if from_ix.is_none() || to_ix.is_none() {
                            links_to_add.push(link);
                            continue;
                        }
                        
                        let from = graph_segments.end_nodes[*from_ix.unwrap()];
                        let to = graph_segments.start_nodes[*to_ix.unwrap()];
                        
                        graph.add_edge(from, to, 0, 1);
                    },
                    _ => (),
                }
            },
            Err(e) => {
                eprintln!("Failed to parse line: {trimmed}");
                eprintln!("Error: {e}")
            },
        }
    }
    
    for link in links_to_add {
        let from_ix = name_to_ix.get(&link.sid1);
        let to_ix = name_to_ix.get(&link.sid2);
        
        if from_ix.is_none() || to_ix.is_none() {
            eprintln!("Omitting link {} -> {} since at least one segment ID does not exists.", link.sid1, link.sid2);
            continue;
        }
        
        let from = graph_segments.end_nodes[*from_ix.unwrap()];
        let to = graph_segments.start_nodes[*to_ix.unwrap()];
        
        graph.add_edge(from, to, 0, 1);
    }
    
    graph.post_process()?;
    
    Ok(POAGraphFromGFA { graph, graph_segments })
}

pub fn format_as_dot<Ix: IndexType>(
    writer: &mut impl fmt::Write,
    graph: &POAGraph<Ix>,
) -> fmt::Result {
    let transformed = graph.graph.map(
        |ix, data| format!("{:?} ({:?})", char::from(data.symbol), ix.index()),
        |_, data| format!("{}, {:?}", data.weight, data.sequence_ids),
    );

    let dot = Dot::new(&transformed);

    writeln!(writer, "{}", dot)?;

    Ok(())
}

pub fn graph_to_gfa<Ix>(writer: &mut impl Write, graph: &POAGraph<Ix>) -> Result<(), PoastaError>
where
    Ix: IndexType + DeserializeOwned,
{
    let mut visited = FxHashSet::default();
    let mut queue = VecDeque::default();
    queue.push_back(graph.start_node());
    visited.insert(graph.start_node());

    writeln!(writer, "H\tVN:Z:1.1")?;

    // Compress non-branching paths into GFA segments
    let mut node_to_segment = FxHashMap::default();
    let mut segment_starts = FxHashMap::default();
    let mut segment_ends = FxHashMap::default();
    let mut segment_lengths = FxHashMap::default();
    let mut curr_segment_id = 0;
    while let Some(front) = queue.pop_front() {
        if front == graph.start_node() {
            for succ in graph.successors(front) {
                if !visited.contains(&succ) {
                    queue.push_back(succ);
                    visited.insert(succ);
                }
            }
        } else {
            let mut segment = vec![graph.get_symbol(front)];
            let mut curr_node = front;
            let mut curr_out_degree = graph.out_degree(front);

            let mut seg_pos = 0usize;
            node_to_segment.insert(front, (curr_segment_id, seg_pos));
            segment_starts.insert(front, curr_segment_id);
            while curr_out_degree == 1 {
                let next_node = graph.successors(curr_node).next().unwrap();
                let in_degree_next = graph.in_degree(next_node);

                if in_degree_next == 1 && next_node != graph.end_node() {
                    segment.push(graph.get_symbol(next_node));
                    node_to_segment.insert(next_node, (curr_segment_id, seg_pos));
                } else {
                    break;
                }

                curr_node = next_node;
                curr_out_degree = graph.out_degree(curr_node);
                seg_pos += 1;
            }

            writeln!(
                writer,
                "S\ts{curr_segment_id}\t{}",
                String::from_utf8_lossy(&segment)
            )?;
            segment_ends.insert(curr_node, curr_segment_id);
            segment_lengths.insert(curr_segment_id, segment.len());
            visited.insert(curr_node);

            for succ in graph.successors(curr_node) {
                if !visited.contains(&succ) && succ != graph.end_node() {
                    visited.insert(succ);
                    queue.push_back(succ);
                }
            }

            curr_segment_id += 1;
        }
    }

    // Add links between segments
    for edge in graph.graph.edge_references() {
        if segment_ends.contains_key(&edge.source()) && segment_starts.contains_key(&edge.target()) {
            let src = segment_ends[&edge.source()];
            let target = segment_starts[&edge.target()];
            writeln!(writer, "L\ts{src}\t+\ts{target}\t+\t0M")?;
        }
    }

    // Add walks indicating each individual aligned sequence
    for (seq_id, seq) in graph.sequences.iter().enumerate() {
        let mut curr = Some(seq.start_node());
        let (mut prev_segment, start_pos) = node_to_segment[&seq.start_node()];
        let mut walk_segments = vec![prev_segment];
        let mut last_pos = 0;
        let mut total_segments_length = segment_lengths[&prev_segment];
        while let Some(n) = curr {
            let node_segment;
            (node_segment, last_pos) = node_to_segment[&n];

            if node_segment != prev_segment {
                walk_segments.push(node_segment);
                total_segments_length += segment_lengths[&node_segment];
            }

            curr = None;
            for out_edge in graph.graph.edges(n) {
                if out_edge
                    .weight()
                    .sequence_ids
                    .binary_search(&seq_id)
                    .is_ok()
                {
                    curr = Some(out_edge.target())
                }
            }

            prev_segment = node_segment;
        }

        // End position with respect to total path length of all segments concatenated
        let end_pos = total_segments_length - segment_lengths[&prev_segment] + last_pos;
        writeln!(
            writer,
            "W\t*\t0\t{}\t{start_pos}\t{end_pos}\t{}",
            seq.name(),
            walk_segments
                .into_iter()
                .map(|v| format!(">s{v}"))
                .collect::<Vec<String>>()
                .join("")
        )?
    }

    Ok(())
}

pub fn graph_to_dot<Ix>(writer: &mut impl Write, graph: &POAGraph<Ix>) -> Result<(), PoastaError>
where
    Ix: IndexType + DeserializeOwned,
{
    let seq_names_str = graph
        .sequences
        .iter()
        .map(|v| format!("{}:{}", v.name(), v.start_node().index()))
        .collect::<Vec<String>>()
        .join("\t");

    writeln!(writer, "# seq:\t{seq_names_str}")?;

    writeln!(writer, "digraph {{")?;
    writeln!(writer, "rankdir=\"LR\"")?;
    writeln!(
        writer,
        "node [shape=square, style=filled, fillcolor=\"#e3e3e3\", penwidth=0]"
    )?;
    writeln!(writer)?;

    for n in graph.all_nodes() {
        writeln!(
            writer,
            "{} [label=\"{}\"; fontcolor=\"{}\"]",
            n.index(),
            graph.get_symbol_char(n),
            graphviz_node_color(graph.get_symbol(n))
        )?;
    }

    let mut aligned_nodes_processed = FxHashSet::default();
    for n in graph.all_nodes() {
        if aligned_nodes_processed.contains(&n) {
            continue;
        }

        let mut node_list = vec![n];
        node_list.extend(graph.graph[n].aligned_nodes.iter().copied());

        if node_list.len() > 1 {
            let node_list_str = node_list
                .iter()
                .map(|v| format!("{}", v.index()))
                .collect::<Vec<String>>()
                .join("; ");

            writeln!(writer, "{{rank=same; {node_list_str}}}")?;
        }

        aligned_nodes_processed.extend(node_list);
    }

    let max_num_seq = graph
        .graph
        .edge_references()
        .map(|e| e.weight().sequence_ids.len())
        .max()
        .unwrap_or(1);
    let min_weight = 1.0;
    let max_weight = 40.0;
    let min_penwidth = 0.5;
    let max_penwidth = 3.5;

    for e in graph.graph.edge_references() {
        let seq_list_str = e
            .weight()
            .sequence_ids
            .iter()
            .map(|v| format!("s{v}"))
            .collect::<Vec<String>>()
            .join(" ");

        let num_seq = e.weight().sequence_ids.len();
        let scaled_weight = (min_weight
            + ((num_seq as f64 / max_num_seq as f64) * (max_weight - min_weight)))
            .round() as i64;
        let scaled_penwidth =
            min_penwidth + ((num_seq as f64 / max_num_seq as f64) * (max_penwidth - min_penwidth));

        writeln!(
            writer,
            "{} -> {} [weight={}; penwidth={}; label={}; class=\"{}\"]",
            e.source().index(),
            e.target().index(),
            scaled_weight,
            scaled_penwidth,
            num_seq,
            seq_list_str
        )?;
    }

    writeln!(writer, "}}")?;
    Ok(())
}

fn graphviz_node_color(label: u8) -> &'static str {
    match label {
        b'A' => "#80BC42",
        b'C' => "#006DB6",
        b'G' => "#F36C3E",
        b'T' => "#B12028",
        _ => "#939393",
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_load_graph_from_gfa() {
        let gfa_path = Path::new("tests/test.gfa");
        
        let POAGraphFromGFA { graph, graph_segments: _ } = load_graph_from_gfa::<u32>(gfa_path).unwrap();
        
        // The test GFA has 4 segments with a total of 35 characters
        // s1→s2→s3→s4 and s2→s4
        // Only s4 has no outgoing edges, so only it connects to end node
        assert_eq!(graph.successors(graph.start_node()).count(), 1);
        assert_eq!(graph.predecessors(graph.end_node()).count(), 1);
        assert_eq!(graph.node_count(), 35);
    }
}