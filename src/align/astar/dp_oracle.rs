//! Plain Gotoh dynamic programming over the full graph × query matrix.
//!
//! This is the score oracle for the A* tests: slow, allocation-heavy, and
//! written for readability over speed. It supports every cost model the
//! crate ships (linear, affine, two-piece affine), in global mode and in the
//! CLI's semi-global mode (graph ends free, query aligned end to end). It
//! computes scores, not backtraces; the A* tests check their own backtraces
//! with `test_support::score_from_alignment`. `cost_from` gives the optimal
//! remaining cost from any alignment state, which checks heuristics for
//! admissibility without running a search.

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

pub(crate) struct DpOracle<C, G> {
    costs: C,
    /// Semi-global: every node except the end sentinel is a free start at
    /// query offset 0, and every node except the start sentinel is a valid
    /// end at the last query offset. Mirrors `AlignmentSpan::SemiGlobal`.
    ends_free: bool,
    _graph: PhantomData<G>,
}

impl<C, G> DpOracle<C, G> {
    pub(crate) fn new(costs: C) -> Self {
        Self { costs, ends_free: false, _graph: PhantomData }
    }

    pub(crate) fn new_semi_global(costs: C) -> Self {
        Self { costs, ends_free: true, _graph: PhantomData }
    }
}

impl<C, G> DpOracle<C, G>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    /// Optimal cost from alignment state `(node, offset, state)` to any valid
    /// end, or `INF` if no end is reachable.
    pub(crate) fn cost_from(&self, graph: &G, query: &[u8], node: G::Node, offset: usize, state: usize) -> u32 {
        self.fill(graph, query, Some((node, offset, state)))
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

    fn align(&mut self, graph: &G, query: &[u8]) -> Result<AlignOutput<G>, Infallible> {
        let score = self.fill(graph, query, None);
        let (n, w) = (graph.node_count(), query.len() + 1);
        Ok(AlignOutput {
            score,
            alignment: vec![],
            stats: AlignmentStats {
                cells_computed: n * w * 5,
                states_expanded: n * w * 5,
                fraction_of_full_matrix: 1.0,
                ..Default::default()
            },
        })
    }
}

impl<C, G> DpOracle<C, G>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    /// Fill the DP table. `from` = `None` starts where the mode starts;
    /// `Some((node, offset, state))` starts at that one state with cost 0.
    fn fill(&self, graph: &G, query: &[u8], from: Option<(G::Node, usize, usize)>) -> u32 {
        let costs = &self.costs;
        let m = query.len();
        let w = m + 1;
        let n = graph.node_count();
        let start = graph.start_node();
        let end = graph.end_node();
        let mismatch = costs.mismatch() as u32;

        // Gap pieces as (open, extend, deletion slot, insertion slot).
        // Single-piece models report 0 for piece two, so a zero
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

            if from.is_none() && v == start {
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

            // Where the search begins: one given state, or in semi-global
            // mode offset 0 of every real node.
            match from {
                Some((fv, fq, fs)) if fv == v => dp[cell(rank, fq)][fs] = 0,
                None if self.ends_free && v != end => dp[cell(rank, 0)][M] = 0,
                _ => {}
            }

            // Moves within the node: insertions, then close every state
            // into M. Insertions are allowed at the end sentinel too, which
            // mirrors `expand_ref_graph_end` in the A* alignment graph.
            for q in 0..=m {
                for &(go, ge, _, ins) in &pieces {
                    if q >= 1 {
                        let prev = dp[cell(rank, q - 1)];
                        let c = &mut dp[cell(rank, q)][ins];
                        *c = (*c).min(sat(prev[M], go + ge)).min(sat(prev[ins], ge));
                    }
                }
                let c = &mut dp[cell(rank, q)];
                for &(_, _, del, ins) in &pieces {
                    c[M] = c[M].min(c[del]).min(c[ins]);
                }
            }
        }

        if self.ends_free {
            (0..n)
                .filter(|&rank| graph.rank_to_node(rank) != start)
                .map(|rank| dp[cell(rank, m)][M])
                .min()
                .unwrap_or(INF)
        } else {
            dp[cell(graph.node_rank(end), m)][M]
        }
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

    #[test]
    fn semi_global_frees_both_graph_ends() {
        let c = Affine::new(4, 6, 2);
        let g = graph(&[b"TTTTACGTACCCC"]);
        let oracle = DpOracle::<Affine, POAGraph<u32>>::new_semi_global(c);
        assert_eq!(oracle.fill(&g, b"ACGTA", None), 0);
        assert_eq!(oracle.fill(&g, b"ACGGTA", None), 8);
        assert!(score(c, &g, b"ACGTA") > 0);
    }

    #[test]
    fn cost_from_start_is_the_global_score() {
        let cases: [(&[&[u8]], &[u8]); 4] = [
            (&[b"ACGT"], b"ACXT"),
            (&[b"ACGTACGTACGT"], b"ACGT"),
            (&[b"ACGAT", b"ACGT"], b"AXGT"),
            (&[b"ACGAT", b"ACCT"], b"ACGGAT"),
        ];
        for (seqs, q) in cases {
            let g = graph(seqs);
            for_costs(&g, q);
        }

        fn for_costs(g: &POAGraph<u32>, q: &[u8]) {
            check(Affine::new(1, 2, 1), g, q);
            check(Linear::new(1, 2), g, q);
            check(TwoPieceAffine::new(1, 2, 2, 6, 1), g, q);
        }

        fn check<C: AlignmentCostModel + Copy>(c: C, g: &POAGraph<u32>, q: &[u8]) {
            use crate::graph::traits::GraphWithNodeOrdering;
            let oracle = DpOracle::<C, POAGraph<u32>>::new(c);
            assert_eq!(oracle.cost_from(g, q, g.start_node(), 0, M), score(c, g, q));
        }
    }

}
