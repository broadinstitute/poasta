//! Read and write POA graphs in GFA 1.0.

use std::fmt::{self, Display, Write as FmtWrite};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use petgraph::Incoming;
use petgraph::graph::IndexType;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use rustc_hash::FxHashMap;

use crate::errors::PoastaIOError;
use crate::graph::consensus::heaviest_bundle_consensus;
use crate::graph::io::CONSENSUS_NAME;
use crate::graph::poa::{POAGraph, POANodeIndex};
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

/// Options controlling GFA output generation.
#[derive(Default, Copy, Clone, Debug)]
pub struct GfaOutputOptions {
    /// When true, append an extra `P` line named "consensus" that traces the graph's
    /// heaviest-bundling consensus path.
    pub include_consensus: bool,
}

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

    // Walk each sequence one node at a time, picking the next edge by binary-searching
    // `seq_idx` in the edge's `sequence_ids` (sorted ascending by `add_edge`).
    // Segments are streamed straight to the writer so we never hold a Vec<String> of
    // the full path in memory.
    for (seq_idx, seq) in graph.sequences.iter().enumerate() {
        write!(out, "P\t{}\t", seq.name())
            .map_err(|source| PoastaIOError::FileWriteError { source })?;

        let mut curr = seq.start_node();
        let mut first = true;
        loop {
            if first {
                write!(out, "n{}+", curr.index())
            } else {
                write!(out, ",n{}+", curr.index())
            }
            .map_err(|source| PoastaIOError::FileWriteError { source })?;
            first = false;

            let next = graph
                .graph
                .edges(curr)
                .find(|e| e.weight().sequence_ids.binary_search(&seq_idx).is_ok())
                .map(|e| e.target());
            match next {
                Some(n) if n != end => curr = n,
                _ => break,
            }
        }

        writeln!(out, "\t*").map_err(|source| PoastaIOError::FileWriteError { source })?;
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

// ─── GFA 1.0 parsing types ────────────────────────────────────────────────────

/// DNA strand orientation used by `L` (link) records.
#[derive(Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// A GFA tag of the form `TAG:TYPE:VALUE`.
#[derive(Debug)]
pub struct Field {
    pub tag: String,
    pub value: FieldValue,
}

impl TryFrom<&str> for Field {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(3, ':');

        let tag = parts.next().ok_or("No tag")?;
        let value_type = parts.next().ok_or("No type")?;
        let value = parts.next().ok_or("No value")?;

        let value = match value_type {
            "A" | "Z" => FieldValue::String(value.to_string()),
            "i" => FieldValue::Integer(value.parse().map_err(|_| "Could not parse integer")?),
            "f" => FieldValue::Float(value.parse().map_err(|_| "Could not parse float")?),
            "J" => FieldValue::Json(value.to_string()),
            "H" => {
                let hex_data: Result<Vec<u8>, _> = value
                    .as_bytes()
                    .chunks(2)
                    .map(|v| {
                        u8::from_str_radix(
                            std::str::from_utf8(v)
                                .map_err(|_| "Could not convert bytes to UTF-8")?,
                            16,
                        )
                        .map_err(|_| "Could not parse hex")
                    })
                    .collect();

                FieldValue::ByteArray(hex_data?)
            }
            _ => return Err("Unsupported tag type"),
        };

        Ok(Field {
            tag: tag.to_string(),
            value,
        })
    }
}

impl Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value_type_tag = match self.value {
            FieldValue::String(_) => "Z",
            FieldValue::Integer(_) => "i",
            FieldValue::Float(_) => "f",
            FieldValue::Json(_) => "J",
            FieldValue::ByteArray(_) => "H",
        };
        write!(f, "{}:{}:{}", self.tag, value_type_tag, self.value)
    }
}

/// A subset of the GFA value types we actually round-trip. The full GFA spec
/// also allows `B` (numeric arrays), which are not needed here.
#[derive(Debug, PartialEq)]
pub enum FieldValue {
    String(String),
    Integer(i64),
    Float(f32),
    Json(String),
    ByteArray(Vec<u8>),
}

impl Display for FieldValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FieldValue::String(s) => write!(f, "{}", s),
            FieldValue::Integer(v) => write!(f, "{}", v),
            FieldValue::Float(v) => write!(f, "{}", v),
            FieldValue::Json(v) => write!(f, "{}", v),
            FieldValue::ByteArray(v) => {
                let mut byte_str = String::with_capacity(v.len() * 2);
                for b in v {
                    let _ = write!(byte_str, "{:02X}", b);
                }
                write!(f, "{}", byte_str)
            }
        }
    }
}

#[derive(Debug)]
pub struct Header {
    pub version: Option<String>,
    pub fields: Vec<Field>,
}

impl TryFrom<&str> for Header {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(3, '\t');

        let start = parts.next().ok_or("Empty line?")?;
        if start != "H" {
            return Err("Not a header line");
        }

        let version = parts
            .next()
            .and_then(|v| Field::try_from(v).ok())
            .and_then(|v| {
                if v.tag == "VN" {
                    match v.value {
                        FieldValue::String(s) => Some(s),
                        _ => None,
                    }
                } else {
                    None
                }
            });

        let fields = if let Some(extra_fields) = parts.next() {
            extra_fields
                .split('\t')
                .filter_map(|v| Field::try_from(v).ok())
                .collect::<Vec<_>>()
        } else {
            Vec::default()
        };

        Ok(Header { version, fields })
    }
}

/// A GFA `S` (segment) record.
pub struct Segment {
    pub sid: String,
    pub sequence: Option<String>,
    pub fields: Vec<Field>,
}

impl TryFrom<&str> for Segment {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(4, '\t');

        let start = parts.next().ok_or("Empty line?")?;
        if start != "S" {
            return Err("Not a segment line");
        }

        let sid = parts.next().ok_or("No sid")?;

        let sequence = parts.next().and_then(|v| {
            if v != "*" {
                Some(v.to_ascii_uppercase())
            } else {
                None
            }
        });

        let fields = if let Some(extra_fields) = parts.next() {
            extra_fields
                .split('\t')
                .filter_map(|v| Field::try_from(v).ok())
                .collect::<Vec<_>>()
        } else {
            Vec::default()
        };

        Ok(Segment {
            sid: sid.to_string(),
            sequence,
            fields,
        })
    }
}

/// A GFA `L` (link) record.
#[derive(Debug)]
pub struct Link {
    pub sid1: String,
    pub sid2: String,
    pub strand1: Strand,
    pub strand2: Strand,
    pub overlap: Option<String>,
    pub fields: Vec<Field>,
}

impl TryFrom<&str> for Link {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(6, '\t');

        let start = parts.next().ok_or("Empty line?")?;
        if start != "L" {
            return Err("Not a link line");
        }

        let sid1 = parts.next().ok_or("Missing first segment ID")?;
        let strand1_str = parts.next().ok_or("Missing first strand")?;
        let sid2 = parts.next().ok_or("Missing second segment ID")?;
        let strand2_str = parts.next().ok_or("Missing second strand")?;
        let overlap_str = parts.next().ok_or("Missing overlap")?;

        let strand1 = match strand1_str {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => return Err("Invalid first strand orientation"),
        };

        let strand2 = match strand2_str {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => return Err("Invalid second strand orientation"),
        };

        let overlap = if overlap_str == "*" {
            None
        } else {
            Some(overlap_str.to_string())
        };

        let fields = if let Some(extra_fields) = parts.next() {
            extra_fields
                .split('\t')
                .filter_map(|v| Field::try_from(v).ok())
                .collect::<Vec<_>>()
        } else {
            Vec::default()
        };

        Ok(Link {
            sid1: sid1.to_string(),
            sid2: sid2.to_string(),
            strand1,
            strand2,
            overlap,
            fields,
        })
    }
}

pub enum GfaLine {
    Header(Header),
    Segment(Segment),
    Link(Link),
    Other(String),
}

impl TryFrom<&str> for GfaLine {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let first_char = value.chars().next().ok_or("Empty line")?;

        match first_char {
            'H' => Ok(GfaLine::Header(Header::try_from(value)?)),
            'S' => Ok(GfaLine::Segment(Segment::try_from(value)?)),
            'L' => Ok(GfaLine::Link(Link::try_from(value)?)),
            '#' | 'C' | 'P' | 'J' | 'W' => Ok(GfaLine::Other(value.to_string())),
            _ => Err("Invalid GFA line type"),
        }
    }
}

// ─── GFA 1.0 reader ────────────────────────────────────────────────────────────

/// Segment-level metadata kept alongside the constructed POA graph. Used by
/// downstream tools (e.g. lasagna) to translate node indices back to GFA
/// segment names + positions for GAF output.
#[derive(Debug)]
pub struct GraphSegments<Ix>
where
    Ix: IndexType,
{
    pub names: Vec<String>,
    pub start_nodes: Vec<POANodeIndex<Ix>>,
    pub end_nodes: Vec<POANodeIndex<Ix>>,
    pub segment_lengths: Vec<usize>,
}

impl<Ix> Default for GraphSegments<Ix>
where
    Ix: IndexType,
{
    fn default() -> Self {
        Self {
            names: Vec::new(),
            start_nodes: Vec::new(),
            end_nodes: Vec::new(),
            segment_lengths: Vec::new(),
        }
    }
}

/// A POA graph loaded from GFA, with the segment metadata required to write
/// GAF output back out.
pub struct POAGraphFromGFA<Ix>
where
    Ix: IndexType,
{
    pub graph: POAGraph<Ix>,
    pub graph_segments: GraphSegments<Ix>,
}

/// Parse a GFA 1.0 file into a `POAGraph`. The file may be plain text or
/// gzipped (`.gz`). `post_process()` is run before returning, so the loaded
/// graph is immediately ready to align against.
///
/// Reverse-strand links are rejected; the POA graph is undirected-strand by
/// construction and we don't attempt to fold the reverse complement here.
pub fn load_graph_from_gfa<Ix>(path: impl AsRef<Path>) -> Result<POAGraphFromGFA<Ix>, PoastaIOError>
where
    Ix: IndexType,
{
    let p = path.as_ref();

    let is_gzipped = p
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    let file = File::open(p).map_err(|source| PoastaIOError::FileReadError { source })?;
    let reader_inner: BufReader<Box<dyn Read>> = BufReader::new(if is_gzipped {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    });

    let mut graph = POAGraph::<Ix>::new();
    let mut graph_segments = GraphSegments::<Ix>::default();
    let mut name_to_ix: FxHashMap<String, usize> = FxHashMap::default();
    let mut links_to_add: Vec<Link> = Vec::new();

    for line in reader_inner.lines() {
        let line = line.map_err(|source| PoastaIOError::FileReadError { source })?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        match GfaLine::try_from(trimmed) {
            Ok(gfa_line) => match gfa_line {
                GfaLine::Segment(segment) => {
                    if let Some(seq) = segment.sequence {
                        let weights = vec![1u32; seq.len()];
                        let (start, end) = graph
                            .add_nodes_for_sequence(seq.as_bytes(), &weights, 0, seq.len())
                            .ok_or(PoastaIOError::InvalidFormat)?;

                        name_to_ix.insert(segment.sid.clone(), graph_segments.names.len());
                        graph_segments.names.push(segment.sid.clone());
                        graph_segments.start_nodes.push(start);
                        graph_segments.end_nodes.push(end);
                        graph_segments.segment_lengths.push(seq.len());
                    } else {
                        tracing::warn!(
                            sid = %segment.sid,
                            "omitting segment with no sequence"
                        );
                    }
                }
                GfaLine::Link(link) => {
                    if link.strand1 == Strand::Reverse || link.strand2 == Strand::Reverse {
                        tracing::error!(?link, "reverse-strand links are not supported");
                        return Err(PoastaIOError::InvalidFormat);
                    }

                    let from_ix = name_to_ix.get(&link.sid1).copied();
                    let to_ix = name_to_ix.get(&link.sid2).copied();

                    // Either segment may appear after the link line. Defer.
                    if from_ix.is_none() || to_ix.is_none() {
                        links_to_add.push(link);
                        continue;
                    }

                    let from = graph_segments.end_nodes[from_ix.unwrap()];
                    let to = graph_segments.start_nodes[to_ix.unwrap()];

                    graph.add_edge(from, to, 0, 1);
                }
                _ => (),
            },
            Err(e) => {
                tracing::warn!(line = %trimmed, error = %e, "failed to parse GFA line");
            }
        }
    }

    for link in links_to_add {
        let from_ix = name_to_ix.get(&link.sid1).copied();
        let to_ix = name_to_ix.get(&link.sid2).copied();

        if from_ix.is_none() || to_ix.is_none() {
            tracing::warn!(
                from = %link.sid1,
                to = %link.sid2,
                "omitting link: segment ID does not exist"
            );
            continue;
        }

        let from = graph_segments.end_nodes[from_ix.unwrap()];
        let to = graph_segments.start_nodes[to_ix.unwrap()];

        graph.add_edge(from, to, 0, 1);
    }

    graph
        .post_process()
        .map_err(|_| PoastaIOError::InvalidFormat)?;

    tracing::info!(
        nodes = graph.node_count(),
        segments = graph_segments.names.len(),
        "loaded GFA graph"
    );

    Ok(POAGraphFromGFA {
        graph,
        graph_segments,
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::traits::GraphBase;

    #[test]
    fn test_parse_header() {
        let h = Header::try_from("H\tVN:Z:1.0").unwrap();
        assert_eq!(h.version, Some("1.0".to_string()));
    }

    #[test]
    fn test_parse_segment() {
        let s = Segment::try_from("S\tseg1\tACGT").unwrap();
        assert_eq!(s.sid, "seg1");
        assert_eq!(s.sequence, Some("ACGT".to_string()));
    }

    #[test]
    fn test_parse_link() {
        let l = Link::try_from("L\ta\t+\tb\t+\t0M").unwrap();
        assert_eq!(l.sid1, "a");
        assert_eq!(l.sid2, "b");
        assert_eq!(l.strand1, Strand::Forward);
    }

    #[test]
    fn test_load_test_gfa_round_trip() {
        let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test.gfa");
        let POAGraphFromGFA {
            graph,
            graph_segments,
        } = load_graph_from_gfa::<u32>(&p).expect("load tests/test.gfa");

        assert_eq!(graph_segments.names, vec!["1", "2", "3", "4", "5", "6"]);
        // tests/test.gfa segment lengths: 7, 2, 3, 10, 2, 2
        assert_eq!(graph_segments.segment_lengths, vec![7, 2, 3, 10, 2, 2]);

        // Total real nodes = sum of segment lengths.
        let total_real: usize = graph_segments.segment_lengths.iter().sum();
        // + 2 sentinel nodes (start/end).
        assert_eq!(graph.node_count(), total_real + 2);
    }
}
