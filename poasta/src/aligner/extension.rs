use tracing::debug;

use super::{
    astar::AlignableGraph,
    fr_points::{to_node_pos, Diag, DiagType, PosType},
};

/// SIMD-accelerated longest common prefix length calculator
///
/// Source: https://matklad.github.io/2023/04/09/can-you-trust-a-compiler-to-optimize-your-code.html
fn common_prefix_len(x: &[u8], y: &[u8]) -> usize {
    let chunk_size = 16;

    let off = std::iter::zip(x.chunks_exact(chunk_size), y.chunks_exact(chunk_size))
        .take_while(|(xs_chunk, ys_chunk)| xs_chunk == ys_chunk)
        .count()
        * chunk_size;

    off + std::iter::zip(&x[off..], &y[off..])
        .take_while(|(x, y)| x == y)
        .count()
}

pub fn extend<G, D, O>(
    graph: &G,
    seq: &[u8],
    node: G::NodeType,
    diag: Diag<D>,
    offset: O,
) -> O
where
    G: AlignableGraph,
    D: DiagType,
    O: PosType,
{
    let node_seq = graph.node_seq(node);
    let node_len = graph.node_len(node);
    let node_pos = to_node_pos(diag, offset.as_usize());
    
    debug!("Checking on node {node:?}, length: {node_len}, current pos: {node_pos}, query offset: {offset:?}");

    if node_pos >= node_len - 1 {
        debug!("- at node end, return {offset:?}");
        return offset;
    }

    if offset.as_usize() == seq.len() {
        debug!("- at qry end, return {offset:?}");
        return offset;
    }

    let lcp = common_prefix_len(&node_seq[node_pos+1..], &seq[offset.as_usize()..]);
    let new_offset = offset.as_usize() + lcp;
    
    O::new(new_offset)
}
