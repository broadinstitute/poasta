//! Consensus sequence generation for POA graphs.
//!
//! Implements Lee's heaviest-bundling algorithm: a forward pass over the POA DAG in
//! topological order where each node follows its heaviest incoming edge. Ties go to the
//! predecessor with the higher score, then to the smaller rank for determinism. A node's
//! score is its predecessor's score plus that edge's weight. Backtracking from the end
//! sentinel yields the consensus path.
//!
//! Following the heaviest edge rather than maximising the summed weight keeps a long
//! detour that few sequences take from outscoring the edge most sequences take.

use crate::graph::traits::IndexType;

use crate::graph::poa::{POAGraph, POANodeIndex};
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

/// Compute the heaviest-bundling consensus path through the graph.
///
/// Returns the ordered list of real (non-sentinel) nodes forming the consensus. Returns
/// an empty vector when the graph has no real nodes.
pub fn heaviest_bundle_consensus<Ix>(graph: &POAGraph<Ix>) -> Vec<POANodeIndex<Ix>>
where
    Ix: IndexType,
{
    if graph.is_empty() {
        return Vec::new();
    }

    let start = graph.start_node();
    let end = graph.end_node();
    let n = graph.node_count();

    let mut score: Vec<i64> = vec![0; n];
    let mut pred: Vec<Option<POANodeIndex<Ix>>> = vec![None; n];

    for rank in 0..n {
        let v = graph.rank_to_node(rank);
        if v == start {
            continue;
        }

        // Compared as (edge weight, predecessor score, Reverse(predecessor rank)).
        let mut best: Option<(i64, i64, usize, POANodeIndex<Ix>)> = None;
        for (u, edge) in graph.in_edges(v) {
            let u_rank = graph.node_rank(u);
            let w = edge.weight as i64;

            let replace = match best {
                None => true,
                Some((best_w, best_score, best_rank, _)) => {
                    (w, score[u_rank], std::cmp::Reverse(u_rank))
                        > (best_w, best_score, std::cmp::Reverse(best_rank))
                }
            };
            if replace {
                best = Some((w, score[u_rank], u_rank, u));
            }
        }

        if let Some((w, s, _, u)) = best {
            score[rank] = s + w;
            pred[rank] = Some(u);
        }
    }

    let mut path: Vec<POANodeIndex<Ix>> = Vec::new();
    let mut curr = end;
    while let Some(p) = pred[graph.node_rank(curr)] {
        path.push(curr);
        curr = p;
    }

    path.reverse();
    path.retain(|n| *n != start && *n != end);
    path
}

/// Drop nodes from both ends of a consensus path until each end node is covered
/// by at least `min_frac` of the sequences in the graph. A node's coverage is the
/// number of sequences through it. Interior nodes are kept, so the result stays a
/// path in the graph. Returns an empty slice when no node is covered enough.
pub fn trim_consensus_ends<'a, Ix>(
    graph: &POAGraph<Ix>,
    path: &'a [POANodeIndex<Ix>],
    min_frac: f64,
) -> &'a [POANodeIndex<Ix>]
where
    Ix: IndexType,
{
    if min_frac <= 0.0 {
        return path;
    }

    let min_coverage = min_frac * graph.sequences.len() as f64;

    // A sequence's first node has no incoming edge with its id, so count those
    // sequences separately.
    let mut starts = vec![0usize; graph.node_count()];
    for seq in &graph.sequences {
        starts[graph.node_rank(graph.node(seq.start_node()))] += 1;
    }
    let covered = |&v: &POANodeIndex<Ix>| {
        let through: usize = graph
            .in_edges(v)
            .map(|(_, edge)| edge.sequence_ids.len())
            .sum();
        (starts[graph.node_rank(v)] + through) as f64 >= min_coverage
    };

    match (
        path.iter().position(covered),
        path.iter().rposition(covered),
    ) {
        (Some(first), Some(last)) => &path[first..=last],
        _ => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::io::fasta::load_graph_from_fasta_msa;

    fn graph_from_msa(rows: &[&str]) -> POAGraph<u32> {
        let mut buf = String::new();
        for (i, r) in rows.iter().enumerate() {
            buf.push_str(&format!(">s{i}\n{r}\n"));
        }
        load_graph_from_fasta_msa::<u32, _>(buf.as_bytes()).unwrap()
    }

    #[test]
    fn consensus_of_identical_sequences_equals_input() {
        let graph = graph_from_msa(&["ACGTACGT", "ACGTACGT", "ACGTACGT"]);

        let consensus = heaviest_bundle_consensus(&graph);
        let symbols: Vec<u8> = consensus.iter().map(|&n| graph.node_symbol(n)).collect();

        assert_eq!(symbols, b"ACGTACGT");
    }

    #[test]
    fn consensus_picks_majority_branch() {
        let graph = graph_from_msa(&["ACGTACGT", "ACGTACGT", "ACGTACGT", "ACTTACGT"]);

        let consensus = heaviest_bundle_consensus(&graph);
        let symbols: Vec<u8> = consensus.iter().map(|&n| graph.node_symbol(n)).collect();

        assert_eq!(symbols, b"ACGTACGT");
    }

    #[test]
    fn consensus_skips_long_minority_insertion() {
        // One of four sequences has a 10 base insertion. Its detour sums to more
        // weight than the majority edge, but that edge is the heavier one.
        let graph = graph_from_msa(&[
            "ACGT----------ACGT",
            "ACGT----------ACGT",
            "ACGT----------ACGT",
            "ACGTTTTTTTTTTTACGT",
        ]);

        let consensus = heaviest_bundle_consensus(&graph);
        let symbols: Vec<u8> = consensus.iter().map(|&n| graph.node_symbol(n)).collect();

        assert_eq!(symbols, b"ACGTACGT");
    }

    /// The example from issue #37: one of 15 sequences has an extra base at one end.
    fn overhang_graph(prefix: &str, suffix: &str) -> POAGraph<u32> {
        let base = "CAG".repeat(20);
        let mut rows = vec![format!("-{base}-"); 14];
        rows.push(format!("{prefix}{base}{suffix}"));
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        graph_from_msa(&rows)
    }

    #[test]
    fn trim_drops_minority_overhangs() {
        for graph in [overhang_graph("-", "C"), overhang_graph("T", "-")] {
            let consensus = heaviest_bundle_consensus(&graph);
            assert_eq!(consensus.len(), 61);
            assert_eq!(trim_consensus_ends(&graph, &consensus, 0.0).len(), 61);

            let trimmed = trim_consensus_ends(&graph, &consensus, 0.35);
            let symbols: Vec<u8> = trimmed.iter().map(|&n| graph.node_symbol(n)).collect();
            assert_eq!(symbols, b"CAG".repeat(20));
        }
    }

    #[test]
    fn trim_can_drop_everything() {
        let graph = graph_from_msa(&["AC--", "--GT"]);
        let consensus = heaviest_bundle_consensus(&graph);
        assert!(trim_consensus_ends(&graph, &consensus, 1.0).is_empty());
    }

    #[test]
    fn consensus_empty_graph() {
        let graph = POAGraph::<u32>::new();
        let consensus = heaviest_bundle_consensus(&graph);
        assert!(consensus.is_empty());
    }
}
