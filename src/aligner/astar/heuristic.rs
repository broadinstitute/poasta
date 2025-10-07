use tracing::{debug, span, Level};

use crate::aligner::cost_models::affine::{self, Affine, AffineAstarItem, AffineAstarState};
use crate::aligner::cost_models::AlignmentCostModel;
use crate::aligner::fr_points::{Diag, DiagType, PosType};
use crate::aligner::traits::AlignableGraph;
use crate::aligner::AlignmentMode;
use crate::graph::bubbles::index::BubbleIndex;
use std::sync::Arc;

use super::runnable::{self, AstarRunnable};
use super::AlignState;

pub trait AstarHeuristic<C, G>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    fn init(
        &self,
        graph: &G,
        seq: &[u8],
        bubble_index: Arc<BubbleIndex>,
        alignment_mode: AlignmentMode,
    ) -> AstarRunnable<
        C,
        G,
        impl Fn(&C::AstarStateType<G>, &C::Item) -> usize,
        impl Fn(&C::AstarStateType<G>, &C::Item) -> bool,
    >;
}

pub struct Dijkstra<C> {
    cost_model: C,
}

impl<C> Dijkstra<C> {
    pub fn new(cost_model: C) -> Self {
        Self { cost_model }
    }
}

impl<D, O, G> AstarHeuristic<Affine<D, O>, G> for Dijkstra<Affine<D, O>>
where
    D: DiagType,
    O: PosType,
    G: AlignableGraph,
{
    fn init(
        &self,
        graph: &G,
        seq: &[u8],
        index: Arc<BubbleIndex>,
        mode: AlignmentMode,
    ) -> AstarRunnable<
        Affine<D, O>,
        G,
        impl Fn(&AffineAstarState<G, D, O>, &AffineAstarItem<D>) -> usize,
        impl Fn(&AffineAstarState<G, D, O>, &AffineAstarItem<D>) -> bool,
    > {
        let astar_state = self.cost_model.init_astar(graph, seq, index.clone(), mode);
        let index_for_prune = index.clone();

        runnable::create(
            astar_state,
            // TODO: other alignment modes (semi-global, ...)
            &[AffineAstarItem::default()],
            |_, _| 0,
            move |state, item| !affine::can_improve_bubble(graph, &index_for_prune, state, item)
        )
    }
}

pub struct MinGapCost<C> {
    cost_model: C,
}

impl<C> MinGapCost<C>
where
    C: AlignmentCostModel,
{
    pub fn new(cost_model: C) -> Self {
        MinGapCost { cost_model }
    }
}

impl<G, D, O> AstarHeuristic<Affine<D, O>, G> for MinGapCost<Affine<D, O>>
where
    G: AlignableGraph,
    D: DiagType,
    O: PosType,
{
    fn init(
        &self,
        graph: &G,
        seq: &[u8],
        index: Arc<BubbleIndex>,
        mode: AlignmentMode,
    ) -> AstarRunnable<
        Affine<D, O>,
        G,
        impl Fn(&AffineAstarState<G, D, O>, &AffineAstarItem<D>) -> usize,
        impl Fn(&AffineAstarState<G, D, O>, &AffineAstarItem<D>) -> bool,
    > {
        let astar_state = self.cost_model.init_astar(graph, seq, index.clone(), mode);

        let end_diag = Diag::new(seq.len() as isize + 1);
        let cost_model = self.cost_model.clone();
        let index_for_prune = index.clone();


        runnable::create(
            astar_state,
            // TODO: other alignment modes (semi-global, ...)
            &[AffineAstarItem::default()],
            move |_, item: &AffineAstarItem<D>| {
                let span = span!(Level::DEBUG, "min_gap_cost");
                let _enter = span.enter();

                // Subtract one to subtract the traversal to the end node
                let min_dist_to_end = index.get_min_dist_to_end(item.node_rank).saturating_sub(1);
                let max_dist_to_end = index.get_max_dist_to_end(item.node_rank).saturating_sub(1);

                let min_end_diag = item.diag + min_dist_to_end;
                let max_end_diag = item.diag + max_dist_to_end;

                let mut aln_state = item.state;
                let gap_length = if min_end_diag > end_diag {
                    // requires deletions, so if in insertion or match state, we need to open a new gap
                    if aln_state != AlignState::Deletion {
                        aln_state = AlignState::Match;
                    }

                    (min_end_diag - end_diag).as_usize()
                } else if max_end_diag < end_diag {
                    // requires insertions, so if in deletion or match state, we need to open a new gap
                    if aln_state != AlignState::Insertion {
                        aln_state = AlignState::Match;
                    }

                    (end_diag - max_end_diag).as_usize()
                } else {
                    // Can reach end diagonal without gaps
                    0usize
                };

                debug!(
                    curr_diag = ?item.diag,
                    end_diag = ?end_diag,
                    min_dist_to_end = min_dist_to_end,
                    max_dist_to_end = max_dist_to_end,
                    gap_length = gap_length
                );

                cost_model.gap_cost(aln_state, gap_length)
            },
            move |state, item| !affine::can_improve_bubble(graph, &index_for_prune, state, item)
        )
    }
}
