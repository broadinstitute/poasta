//! DOT (Graphviz) serialization for POA graphs.

use std::io::{self, Write};

use petgraph::graph::IndexType;
use petgraph::visit::{EdgeRef, IntoEdgeReferences, IntoNodeReferences};

use crate::graph::poa::POAGraph;

/// Write a POA graph in Graphviz DOT format.
///
/// The first line is a DOT comment containing the sequence to be aligned,
/// so it can be parsed by external scripts.
///
/// Node attributes:
/// - `label`: node symbol character
/// - `xlabel`: topological rank (rendered as an external label)
/// - `lmin`, `lmax`: shortest/longest path from start (non-visible custom attributes)
///
/// Edge attributes:
/// - `label`: space-separated sequence IDs
/// - `penwidth`: logarithmically scaled by number of sequences, clamped to [1.0, 8.0]
pub fn poa_graph_to_dot<Ix, W>(
    graph: &POAGraph<Ix>,
    mut writer: W,
    seq_to_align: &[u8],
) -> io::Result<()>
where
    Ix: IndexType,
    W: Write,
{
    let seq_str = String::from_utf8_lossy(seq_to_align);
    writeln!(writer, "// sequence: {seq_str}")?;
    writeln!(writer, "digraph poa {{")?;
    writeln!(writer, "  rankdir=LR;")?;

    // Nodes
    for (node, data) in graph.graph.node_references() {
        let symbol = data.symbol as char;
        let rank = data.rank.index();
        let lmin = data.depth_min.index();
        let lmax = data.depth_max.index();
        let id = node.index();
        writeln!(
            writer,
            r#"  N{id} [label="{symbol}" xlabel="{rank}" lmin={lmin} lmax={lmax}];"#
        )?;
    }

    // Edges
    for edge in graph.graph.edge_references() {
        let src = edge.source().index();
        let dst = edge.target().index();
        let data = edge.weight();
        let n_seqs = data.sequence_ids.len();
        let penwidth: f64 = if n_seqs == 0 {
            1.0
        } else {
            ((1.0_f64 + n_seqs as f64).log2() * 2.0).clamp(1.0, 8.0)
        };
        let seq_ids: Vec<String> = data.sequence_ids.iter().map(|id| id.to_string()).collect();
        let label = seq_ids.join(" ");
        writeln!(
            writer,
            r#"  N{src} -> N{dst} [label="{label}" penwidth={penwidth:.2}];"#
        )?;
    }

    writeln!(writer, "}}")?;
    Ok(())
}
