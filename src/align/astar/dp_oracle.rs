//! Plain Gotoh dynamic programming over the full graph × query matrix.
//!
//! This is the score oracle for the A* tests: slow, allocation-heavy, and
//! written for readability over speed. It supports every cost model the
//! crate ships (linear, affine, two-piece affine) and global alignment only.
//! It computes scores, not backtraces; the A* tests check their own
//! backtraces with `test_support::score_from_alignment`.

use std::convert::Infallible;
use std::marker::PhantomData;

use crate::align::cost_models::AlignmentCostModel;
use crate::align::traits::{AlignOutput, AlignmentStats};
use crate::align::score::{INF, sat};
use crate::align::traits::{AlignableGraph, AlignmentEngine};

/// Slot layout of one DP cell: `[M, D1, I1, D2, I2]`, the same order as
/// `AlignState as usize`. `M` is the closed state: the minimum over the
/// diagonal move and every gap state, which mirrors the zero-cost
/// gap-to-match edges in the A* alignment graph.
const M: usize = 0;

pub(crate) struct DpOracle<C, G>(C, PhantomData<G>);

impl<C, G> DpOracle<C, G> {
    pub(crate) fn new(costs: C) -> Self {
        Self(costs, PhantomData)
    }
}

impl<C, G> AlignmentEngine<&[u8]> for DpOracle<C, G>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    type Graph = G;
    type Success = AlignOutput<G>;
    type Error = Infallible;

    fn align(&self, graph: &G, query: &[u8]) -> Result<AlignOutput<G>, Infallible> {
        let costs = &self.0;
        let m = query.len();
        let w = m + 1;
        let n = graph.node_count();
        let start = graph.start_node();
        let end = graph.end_node();
        let mismatch = costs.mismatch() as u32;

        // Gap pieces as (open, extend, deletion slot, insertion slot).
        // ponytail: single-piece models report 0 for piece two, so a zero
        // extend cost is the "no second piece" signal.
        let mut pieces = vec![(
            costs.gap_open() as u32,
            costs.gap_extend() as u32,
            1usize,
            2usize,
        )];
        if costs.gap_extend2() > 0 {
            pieces.push((costs.gap_open2() as u32, costs.gap_extend2() as u32, 3, 4));
        }

        let mut dp = vec![[INF; 5]; n * w];
        let cell = |rank: usize, q: usize| rank * w + q;

        for rank in 0..n {
            let v = graph.rank_to_node(rank);

            if v == start {
                dp[cell(rank, 0)][M] = 0;
                for q in 1..=m {
                    for &(go, ge, _, ins) in &pieces {
                        dp[cell(rank, q)][ins] = go + q as u32 * ge;
                    }
                    dp[cell(rank, q)][M] = pieces
                        .iter()
                        .map(|p| dp[cell(rank, q)][p.3])
                        .min()
                        .unwrap();
                }
                continue;
            }

            // Moves that come from predecessors: deletions and the diagonal.
            // The end sentinel consumes no query base, so its "diagonal" is a
            // free horizontal move, as in the A* DFA.
            for pred in graph.predecessors(v) {
                let pr = graph.node_rank(pred);
                for q in 0..=m {
                    let p = dp[cell(pr, q)];
                    let c = &mut dp[cell(rank, q)];
                    for &(go, ge, del, _) in &pieces {
                        c[del] = c[del].min(sat(p[M], go + ge)).min(sat(p[del], ge));
                    }
                    if v == end {
                        c[M] = c[M].min(p[M]);
                    } else if q >= 1 {
                        let p_diag = dp[cell(pr, q - 1)][M];
                        let step = if graph.is_symbol_equal(v, query[q - 1]) {
                            0
                        } else {
                            mismatch
                        };
                        dp[cell(rank, q)][M] = dp[cell(rank, q)][M].min(sat(p_diag, step));
                    }
                }
            }

            // Moves within the node: insertions, then close every state
            // into M. Insertions are allowed at the end sentinel too, which
            // mirrors `expand_ref_graph_end` in the A* alignment graph.
            for q in 0..=m {
                for &(go, ge, _, ins) in &pieces {
                    if q >= 1 {
                        let prev = dp[cell(rank, q - 1)];
                        dp[cell(rank, q)][ins] = sat(prev[M], go + ge).min(sat(prev[ins], ge));
                    }
                }
                let c = &mut dp[cell(rank, q)];
                for &(_, _, del, ins) in &pieces {
                    c[M] = c[M].min(c[del]).min(c[ins]);
                }
            }
        }

        let score = dp[cell(graph.node_rank(end), m)][M];
        Ok(AlignOutput {
            score,
            alignment: vec![],
            stats: AlignmentStats {
                max_bandwidth: 0,
                cells_computed: n * w * 5,
                fraction_of_full_matrix: 1.0,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cost_models::affine::Affine;
    use crate::align::cost_models::linear::Linear;
    use crate::align::cost_models::two_piece::TwoPieceAffine;
    use crate::graph::alignment::AddAlignment;
    use crate::graph::poa::POAGraph;

    fn graph(seqs: &[&[u8]]) -> POAGraph<u32> {
        let mut g = POAGraph::new();
        for (i, s) in seqs.iter().enumerate() {
            g.add_alignment(&format!("s{i}"), s, None, &vec![1; s.len()])
                .unwrap();
        }
        g
    }

    fn score<C: AlignmentCostModel>(costs: C, g: &POAGraph<u32>, q: &[u8]) -> u32 {
        DpOracle::<C, POAGraph<u32>>::new(costs)
            .align(g, q)
            .unwrap()
            .score
    }

    #[test]
    fn affine_basic_moves() {
        let c = Affine::new(1, 2, 1);
        let g = graph(&[b"ACGT"]);
        assert_eq!(score(c, &g, b"ACGT"), 0);
        assert_eq!(score(c, &g, b"ACXT"), 1);
        assert_eq!(score(c, &g, b"ACT"), 3);
        assert_eq!(score(c, &g, b"ACGGT"), 3);
        assert_eq!(score(c, &g, b"AT"), 4);
    }

    #[test]
    fn linear_gap_has_no_open_cost() {
        let c = Linear::new(1, 2);
        let g = graph(&[b"ACGT"]);
        assert_eq!(score(c, &g, b"AT"), 4);
        assert_eq!(score(c, &g, b"ACGTGG"), 4);
    }

    #[test]
    fn two_piece_long_gap_takes_second_piece() {
        // piece 1: open 2, extend 2; piece 2: open 6, extend 1.
        let c = TwoPieceAffine::new(1, 2, 2, 6, 1);
        let g = graph(&[b"ACGTACGTACGT"]);
        // One-base deletion: piece 1 (4) beats piece 2 (7).
        assert_eq!(score(c, &g, b"ACGTACGTACG"), 4);
        // Eight-base deletion: piece 2 (14) beats piece 1 (18).
        assert_eq!(score(c, &g, b"ACGT"), 14);
    }

    #[test]
    fn bubble_takes_best_branch() {
        let c = Affine::new(1, 2, 1);
        let g = graph(&[b"ACGAT", b"ACGT"]);
        assert_eq!(score(c, &g, b"ACGAT"), 0);
        assert_eq!(score(c, &g, b"ACGT"), 0);
        assert_eq!(score(c, &g, b"AXGT"), 1);
    }

}
