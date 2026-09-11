//! Debug output for inspecting the A* search.
//!
//! [`subscriber::AlignStateLayer`] is a `tracing` layer that turns the A*
//! engine's trace events into per-sequence TSV files. [`write_graph_dot`]
//! dumps the graph a sequence is aligned against, so the two can be plotted
//! together.

pub mod subscriber;

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use petgraph::graph::IndexType;

use crate::graph::io::dot::poa_graph_to_dot;
use crate::graph::poa::POAGraph;

/// Write `dir/graph_for_<seq_name>.dot` for the graph `seq` is about to be
/// aligned to. The name goes through [`sanitise_filename`].
pub fn write_graph_dot<Ix: IndexType>(
    dir: &Path,
    seq_name: &str,
    graph: &POAGraph<Ix>,
    seq: &[u8],
) -> io::Result<()> {
    let path = dir.join(format!("graph_for_{}.dot", sanitise_filename(seq_name)));
    let mut w = BufWriter::new(File::create(path)?);
    poa_graph_to_dot(graph, &mut w, seq)?;
    w.flush()
}

/// Replace characters that are unsafe or inconvenient in filenames with `_`.
pub fn sanitise_filename(name: &str) -> String {
    let mut out: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '.' || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if out.is_empty() {
        out.push('_');
    }
    out
}
