//! Inspection API for tools that show the search, such as poasta-viz.
//! Behind the `viz` feature. Nothing here changes `h`.

use super::contours::{MatchId, Point, Staircase};
use super::{Potential, SmoothHeuristic};
use crate::align::AlignState;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::offsets::OffsetType;
use crate::align::traits::AlignableGraph;
use crate::graph::traits::GraphNodeId;

/// What `h` is made of at one state. For inspecting a search, not used by
/// the engine.
#[derive(Clone, Debug, Default)]
pub struct ChainExplain {
    /// Charge per seed ahead.
    pub charge: i64,
    pub potential: Potential,
    pub t: Point,
    /// Savings of the best chain reachable from the state.
    pub credit: i64,
    /// The match whose sub-point gave the credit, then the matches its
    /// savings came through. Empty if the credit is 0.
    pub chain: Vec<MatchId>,
    pub h_chain: usize,
    pub h_min_gap: usize,
}

/// One seed match. Its id is its index in [`SmoothHeuristic::matches`].
#[derive(Clone, Debug)]
pub struct MatchInfo<N> {
    pub seed: usize,
    pub start: N,
    pub end: N,
    pub t: Point,
    /// Savings of the best chain starting here. 0 once pruned.
    pub sigma: i64,
}

/// Inspection API, for tools that show the search. `h` is unchanged by
/// any of these.
impl<N, O, C> SmoothHeuristic<N, O, C>
where
    N: GraphNodeId,
    O: OffsetType,
    C: AlignmentCostModel + Copy,
{
    /// Whether this alignment built seeds and contours.
    pub fn is_active(&self) -> bool {
        self.index.is_some()
    }

    pub fn n_seeds(&self) -> usize {
        self.index.as_ref().map_or(0, |ix| ix.seeds.n)
    }

    /// The parts of `h` at a state. `h_chain.max(h_min_gap)` equals `h`.
    pub fn explain<G>(
        &mut self,
        graph: &G,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> ChainExplain
    where
        G: AlignableGraph<Node = N>,
    {
        let h_min_gap = self.min_gap.bound(aln_node, aln_state);
        let Some(ix) = &mut self.index else {
            return ChainExplain {
                h_min_gap,
                ..Default::default()
            };
        };
        let (v, i) = (aln_node.node(), aln_node.offset().as_usize());
        let potential = ix.potential(
            graph,
            v,
            i,
            ix.potential_bound(i, aln_state),
            &mut self.scratch.walk,
        );
        let mut e = ChainExplain {
            charge: ix.costs.charge,
            potential,
            h_min_gap,
            ..Default::default()
        };
        let p = potential.value(e.charge);
        if p > 0 {
            e.t = ix.point(v, i, p);
            let best = ix.contours.best(e.t);
            e.credit = best.credit;
            e.chain = best.owner.map_or_else(Vec::new, |m| ix.contours.chain(m));
            e.h_chain = (p - best.credit).max(0) as usize;
        }
        e
    }

    /// T(x) of a state, or None once it has no potential left.
    pub fn t<G>(
        &mut self,
        graph: &G,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> Option<Point>
    where
        G: AlignableGraph<Node = N>,
    {
        let ix = self.index.as_ref()?;
        let (v, i) = (aln_node.node(), aln_node.offset().as_usize());
        let p = ix
            .potential(
                graph,
                v,
                i,
                ix.potential_bound(i, aln_state),
                &mut self.scratch.walk,
            )
            .value(ix.costs.charge);
        (p > 0).then(|| ix.point(v, i, p))
    }

    /// Every seed match, in id order, with exact savings.
    pub fn matches(&mut self) -> Vec<MatchInfo<N>> {
        let Some(ix) = &mut self.index else {
            return Vec::new();
        };
        ix.contours.repair();
        let mut out = Vec::with_capacity(ix.matches.len());
        for seed in 0..ix.seeds.n {
            for m in ix.matches.of_seed(seed) {
                let id = out.len() as MatchId;
                out.push(MatchInfo {
                    seed,
                    start: m.start,
                    end: m.end,
                    t: ix.contours.t(id),
                    sigma: ix.contours.sigma(id),
                });
            }
        }
        out
    }

    /// Delete match m as match pruning would.
    pub fn prune_match(&mut self, m: MatchId) {
        if let Some(ix) = &mut self.index {
            ix.contours.prune(m);
        }
    }

    /// Sub-points by credit, exact.
    pub fn staircases(&mut self) -> Vec<Staircase> {
        let Some(ix) = &mut self.index else {
            return Vec::new();
        };
        ix.contours.repair();
        ix.contours.staircases()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::align::AlignmentMode;
    use crate::align::astar::heuristic::AstarHeuristic;
    use crate::align::astar::heuristic::chain::SeedParams;
    use crate::align::astar::{AstarEngine, Offset};
    use crate::align::cost_models::affine::Affine;
    use crate::graph::bubbles::index::BubbleIndex;
    use crate::graph::poa::POAGraph;

    /// The inspection API must add up to `h`, before and after pruning,
    /// and name a chain exactly when there is credit.
    #[test]
    fn explain_adds_up_to_h() {
        use crate::graph::io::seq::open_sequences;
        use crate::graph::traits::GraphBase;
        let seqs: Vec<Vec<u8>> = open_sequences(std::path::Path::new("tests/dnaG/align_seq.fna"))
            .unwrap()
            .take(4)
            .map(|r| r.unwrap().1[..300].to_vec())
            .collect();
        let mut aligner = crate::align::PoastaAligner::new(
            AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(4, 6, 2)),
            POAGraph::new(),
        );
        for (i, s) in seqs[..3].iter().enumerate() {
            aligner.align_one_named(&format!("s{i}"), s).unwrap();
        }
        let g = aligner.into_graph();
        let q = &seqs[3];
        let mut h = SmoothHeuristic::<_, Offset, _>::new(
            Affine::new(4, 6, 2),
            SeedParams {
                k: 6,
                r: 0,
                prune_matches: true,
            },
        );
        h.prepare(
            &g,
            q,
            &Arc::new(BubbleIndex::new(&g)),
            AlignmentMode::Global,
        );
        let n_matches = h.matches().len();
        assert!(n_matches > 10);
        for round in 0..2 {
            for v in g.all_nodes_iter() {
                for i in (0..=q.len()).step_by(7) {
                    for st in [
                        AlignState::Match,
                        AlignState::Deletion,
                        AlignState::Insertion,
                    ] {
                        let x = AlignmentGraphNode::new(v, i as Offset);
                        let e = h.explain(&g, &x, st);
                        assert_eq!(e.h_chain.max(e.h_min_gap), h.h(&g, &x, st), "round {round}");
                        assert_eq!(e.chain.is_empty(), e.credit == 0);
                    }
                }
            }
            for m in (0..n_matches as u32).step_by(2) {
                h.prune_match(m);
            }
        }
    }
}
