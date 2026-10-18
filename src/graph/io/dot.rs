//! DOT (Graphviz) serialization for POA graphs.

use std::io::{self, Write};

use crate::graph::traits::IndexType;

use crate::graph::poa::POAGraph;
use crate::graph::traits::GraphBase;

/// Write a POA graph in Graphviz DOT format.
///
/// The first line is a DOT comment containing the sequence to be aligned,
/// so it can be parsed by external scripts.
///
/// Nodes are named by rank, as in the debug TSVs of the search.
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
    for node in graph.all_nodes_iter() {
        let symbol = graph.node_symbol(node) as char;
        let (lmin, lmax) = graph.depth(node);
        // The id is the rank. `xlabel` repeats it for Graphviz and the Python
        // tools, which read it as the rank.
        let id = node.index();
        writeln!(
            writer,
            r#"  N{id} [label="{symbol}" xlabel="{id}" lmin={lmin} lmax={lmax}];"#
        )?;
    }

    // Edges
    for (src, dst, data) in graph.all_edges() {
        let (src, dst) = (src.index(), dst.index());
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
