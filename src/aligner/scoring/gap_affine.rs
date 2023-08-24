use std::collections::VecDeque;
use std::fmt::{Display, Formatter};
use smallvec::SmallVec;

use crate::aligner::{AlignedPair, Alignment};
use crate::aligner::layers::Layer;

use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::{Diag, DiagonalPoint, DiagRange, OffsetType};
use crate::aligner::scoring::{AlignmentCosts, AlignmentCostsAffine, AlignmentCostsEdit, AlignmentCostsLinear, LayerCompute};
use crate::aligner::visited::{AlignState, Backtrace, IntervalSmallVec, VisitedIntervalOnDiag,
                              VisitedInterval, VisitedIntervalType, merge_intervals, VisitedIntervalTree};
use crate::aligner::visited::transformations::prelude::*;
use crate::aligner::visited::transformations::{TransformDeletionIterator, TransformInsertionIterator, TransformMismatchIterator};

use crate::debug::DebugOutputWriter;
use crate::debug::messages::DebugOutputMessage;


#[derive(Clone, Copy, Debug)]
pub struct GapAffine {
    cost_mismatch: u8,
    cost_gap_extend: u8,
    cost_gap_open: u8
}

impl GapAffine {
    pub fn new(cost_mismatch: u8, cost_gap_extend: u8, cost_gap_open: u8) -> Self {
        Self { cost_mismatch, cost_gap_extend, cost_gap_open }
    }
}

impl AlignmentCosts for GapAffine {
    type LayerComputeType<'a, O> = GapAffineCompute<'a, O>
        where O: OffsetType;

    fn new_layer_compute<'a, G, O>(&self, graph: &'a G) -> Self::LayerComputeType<'a, O>
    where
        G: AlignableGraph,
        O: OffsetType,
    {
        GapAffineCompute::initial(*self, graph.start_nodes())
    }
}

impl AlignmentCostsEdit for GapAffine {
    #[inline(always)]
    fn mismatch(&self) -> i64 {
        self.cost_mismatch as i64
    }
}

impl AlignmentCostsLinear for GapAffine {
    #[inline(always)]
    fn gap_extend(&self) -> i64 {
        self.cost_gap_extend as i64
    }
}

impl AlignmentCostsAffine for GapAffine {
    #[inline(always)]
    fn gap_open(&self) -> i64 {
        self.cost_gap_open as i64
    }
}

#[derive(Default, Clone, Debug)]
struct GapAffineLayer<O>
where
    O: OffsetType,
{
    layer_m: Layer<O>,
    layer_d: Layer<O>,
    layer_i: Layer<O>,
}

impl<O> GapAffineLayer<O>
where
    O: OffsetType,
{
    pub fn initial<N: NodeIndexType>(start_nodes: &[N]) -> Self {
        Self {
            layer_m: Layer::initial(start_nodes),
            layer_d: Layer::default(),
            layer_i: Layer::default()
        }
    }

    pub fn new(
        k_min: Diag, k_max: Diag,
        visited_m: VecDeque<VisitedIntervalTree<O>>,
        visited_d: VecDeque<VisitedIntervalTree<O>>,
        visited_i: VecDeque<VisitedIntervalTree<O>>,
    ) -> Self {
        Self {
            layer_m: Layer::new(k_min, k_max, visited_m),
            layer_d: Layer::new(k_min, k_max, visited_d),
            layer_i: Layer::new(k_min, k_max, visited_i)
        }
    }

    pub fn new_empty() -> Self {
        Self {
            layer_m: Layer::default(),
            layer_d: Layer::default(),
            layer_i: Layer::default()
        }
    }

    pub fn kmin(&self) -> Option<Diag> {
        [self.layer_m.kmin(), self.layer_d.kmin(), self.layer_i.kmin()]
            .into_iter()
            .flatten()
            .min()
    }

    pub fn kmax(&self) -> Option<Diag> {
        [self.layer_m.kmax(), self.layer_d.kmax(), self.layer_i.kmax()]
            .into_iter()
            .flatten()
            .max()
    }

    pub fn is_empty(&self) -> bool {
        self.layer_m.is_empty() && self.layer_d.is_empty() && self.layer_i.is_empty()
    }

    pub fn get_prev_ival(&self, state: AlignState, diag: Diag, (start, end): (O, O)) -> Option<&VisitedIntervalType<O>> {
        eprintln!("Searching for previous interval on {:?}, search range: ({:?}, {:?}), state: {:?}", diag, start, end, state);
        match state {
            AlignState::Start | AlignState::MisMatch | AlignState::Extended =>
                self.layer_m.get_visited(diag)
                    .and_then(|itree| itree.find(start, end).next()),
            AlignState::Deletion =>
                self.layer_d.get_visited(diag)
                    .and_then(|itree| itree.find(start, end).next()),
            AlignState::Insertion =>
                self.layer_i.get_visited(diag)
                    .and_then(|itree| itree.find(start, end).next()),
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state!")
        }
    }

}

pub struct GapAffineCompute<'a, O>
where
    O: OffsetType,
{
    costs: GapAffine,
    layers: Vec<GapAffineLayer<O>>,
    debug_output: Option<&'a DebugOutputWriter>,
}

impl<'a, O> GapAffineCompute<'a, O>
where
    O: OffsetType,
{
    pub fn initial<N: NodeIndexType>(costs: GapAffine, start_nodes: &[N]) -> Self {
        Self {
            costs,
            layers: vec![GapAffineLayer::initial(start_nodes)],
            debug_output: None,
        }
    }

    pub fn score_delta(&self, curr_state: AlignState, prev_state: AlignState) -> i64 {
        match curr_state {
            AlignState::Extended | AlignState::Start => 0,
            AlignState::MisMatch => match prev_state {
                AlignState::Extended => 0,
                AlignState::Start | AlignState::MisMatch => self.costs.mismatch(),
                AlignState::Insertion | AlignState::Deletion => 0,
                AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state!")
            },
            AlignState::Deletion => match prev_state {
                AlignState::Start | AlignState::Extended | AlignState::MisMatch => self.costs.gap_open() + self.costs.gap_extend(),
                AlignState::Deletion => self.costs.gap_extend(),
                AlignState::Insertion | AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state!"),
            }
            AlignState::Insertion => match prev_state {
                AlignState::Start | AlignState::Extended | AlignState::MisMatch => self.costs.gap_open() + self.costs.gap_extend(),
                AlignState::Insertion => self.costs.gap_extend(),
                AlignState::Deletion | AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state!"),
            },
            AlignState::Deletion2 | AlignState::Insertion2 => panic!("Invalid gap-affine state!")
        }
    }

    fn get_prev_layer(&self, score: i64) -> Option<&GapAffineLayer<O>> {
        if score < 0 {
            return None;
        }

        self.layers.get(score as usize)
    }

    fn get_layer_before(&self, curr_score: i64) -> Option<&GapAffineLayer<O>> {
        let min_cost_delta = [
            self.costs.mismatch(),
            self.costs.gap_open() + self.costs.gap_extend(),
            self.costs.gap_extend()
        ].into_iter().min().unwrap();

        self.get_prev_layer(curr_score - min_cost_delta)
    }

    fn next_kmin<G>(&self, graph: &G, new_score: i64) -> Diag
    where
        G: AlignableGraph,
    {
        let k_min_mismatch = self.get_prev_layer(new_score - self.costs.mismatch())
            .and_then(|layer| layer.layer_m.get_successors_kmin(graph));
        let k_min_gap_open = self.get_prev_layer(new_score - self.costs.gap_open() - self.costs.gap_extend())
            .and_then(|layer| layer.layer_m.get_successors_kmin(graph));
        let k_min_del_ext = self.get_prev_layer(new_score - self.costs.gap_extend())
            .and_then(|layer| layer.layer_d.get_successors_kmin(graph));
        let k_min_ins_ext = self.get_prev_layer(new_score - self.costs.gap_extend())
            .and_then(|layer| layer.layer_i.kmin());

        [k_min_mismatch, k_min_gap_open, k_min_del_ext, k_min_ins_ext]
            .into_iter()
            .flatten()
            .min()
            .unwrap_or(Diag::new(0)) - 1
    }

    fn next_kmax<G>(&self, graph: &G, new_score: i64) -> Diag
    where
        G: AlignableGraph,
    {
        let k_max_mismatch = self.get_prev_layer(new_score - self.costs.mismatch())
            .and_then(|layer| layer.layer_m.get_successors_kmax(graph));
        let k_max_gap_open = self.get_prev_layer(new_score - self.costs.gap_open() - self.costs.gap_extend())
            .and_then(|layer| layer.layer_m.get_successors_kmax(graph));
        let k_max_del_ext = self.get_prev_layer(new_score - self.costs.gap_extend())
            .and_then(|layer| layer.layer_d.get_successors_kmax(graph));
        let k_max_ins_ext = self.get_prev_layer(new_score - self.costs.gap_extend())
            .and_then(|layer| layer.layer_i.kmax());

        [k_max_mismatch, k_max_gap_open, k_max_del_ext, k_max_ins_ext]
            .into_iter()
            .flatten()
            .max()
            .unwrap_or(Diag::new(0)) + 1
    }
}

impl<'a, O> LayerCompute<'a, O> for GapAffineCompute<'a, O>
where
    O: OffsetType,
{
    fn with_debug(mut self, debug_writer: &'a DebugOutputWriter) -> Self {
        self.debug_output = Some(debug_writer);

        self
    }

    fn extend<G>(&mut self, graph: &G, seq: &[u8]) -> Option<(Diag, VisitedIntervalType<O>)>
    where
        G: AlignableGraph,
    {
        let last_layer = self.layers.last_mut().unwrap();
        let end = last_layer.layer_m.extend(graph, seq)
            .map(|(end_diag, end_ival)| (end_diag, end_ival.clone()));

        if let Some(debug_output) = self.debug_output {
            if !last_layer.layer_m.is_empty() {
                debug_output.log(DebugOutputMessage::new_from_layer(AlignState::Extended, &last_layer.layer_m));
            }
        }

        end
    }

    fn next_layer<G>(&mut self, graph: &G, seq: &[u8], new_score: i64)
    where
        G: AlignableGraph,
    {
        let seq_len_as_o = O::new(seq.len());

        let prev_layer_open = self.get_prev_layer(new_score - self.costs.gap_open() - self.costs.gap_extend());
        let prev_layer_ext = self.get_prev_layer(new_score - self.costs.gap_extend());
        let prev_layer_mis = self.get_prev_layer(new_score - self.costs.mismatch());

        if prev_layer_open.map_or(true, |layer| layer.is_empty())
            && prev_layer_ext.map_or(true, |layer| layer.is_empty())
            && prev_layer_mis.map_or(true, |layer| layer.is_empty()) {
            let new_layer = GapAffineLayer::new_empty();
            self.layers.push(new_layer);
            return;
        }

        // Find new k_min and k_max, by inspecting the diagonals of successors from nodes reached
        // in previous layers.
        let new_kmin = self.next_kmin(graph, new_score);
        let new_kmax = self.next_kmax(graph, new_score);
        eprintln!("New k-min: {:?}, k-max: {:?}", new_kmin, new_kmax);

        eprintln!("Building INS:");
        let new_visited_ins: VecDeque<_> = DiagRange::closed(new_kmin, new_kmax).into_iter()
            .map(|diag| {
                // Find source diagonals with intervals that will be merged into the current diagonal.
                // For insertions, there could be only two sources: either we open a new insertion, and
                // thus enter from a previous (mis)match state, or we extend an insertion, and continue
                // fro the Insertion state.
                let source_intervals: SmallVec<[IntervalSmallVec<O>; 3]> = [
                    prev_layer_open
                        .and_then(|layer| layer.layer_m.get_visited(diag - 1))
                        .map(|itree| itree.iter()
                            .cloned_with_backtrace(Backtrace::new(diag - 1, AlignState::MisMatch), true)
                            .transform_insertion()
                            .collect()
                        ),
                    prev_layer_ext
                        .and_then(|layer| layer.layer_i.get_visited(diag - 1))
                        .map(|itree| itree.iter()
                            .cloned_with_backtrace(Backtrace::new(diag - 1, AlignState::Insertion), false)
                            .transform_insertion()
                            .collect()
                        ),
                ].into_iter().flatten().collect();

                eprintln!("Curr diag: {:?}", diag);
                let source_intervals_ref: SmallVec<[&IntervalSmallVec<O>; 3]> = source_intervals.iter().collect();
                VisitedIntervalTree::from_iter(merge_intervals(&source_intervals_ref))
            })
            .collect();

        eprintln!("INS: {:?}", new_visited_ins);
        eprintln!("Building DEL:");

        let new_visited_del: VecDeque<_> = DiagRange::closed(new_kmin, new_kmax).into_iter()
            .map(|diag| {
                let source_ivals_del: SmallVec<[IntervalSmallVec<O>; 8]> = DiagRange::closed(new_kmin, new_kmax).into_iter()
                    // Iterate over all potential source diagonals
                    .filter_map(|src_diag| {
                        prev_layer_ext
                            .and_then(|layer| layer.layer_d.get_visited(src_diag))
                            .and_then(|itree| {
                                let src_ivals: IntervalSmallVec<O> = itree.iter()
                                    // Only include intervals with a successor that enters `diag`
                                    .filter(|ival| {
                                        let pred_row = (src_diag, *ival).end_point_row();
                                        let succ_row = (diag, *ival).end_point_row();

                                        graph.is_successor(pred_row, succ_row)
                                    })
                                    .cloned_with_backtrace(Backtrace::new(src_diag, AlignState::Deletion), false)
                                    .transform_deletion(diag, src_diag)
                                    .collect();

                                if !src_ivals.is_empty() { Some(src_ivals) } else { None }
                            })
                    })
                    .collect();

                let source_ivals_open: SmallVec<[IntervalSmallVec<O>; 8]> = DiagRange::closed(new_kmin, new_kmax).into_iter()
                    .filter_map(|src_diag| {
                        prev_layer_open
                            .and_then(|layer| layer.layer_m.get_visited(src_diag))
                            .and_then(|itree| {
                                let src_ivals: IntervalSmallVec<O> = itree.iter()
                                    // Only include intervals with a successor that enters `diag`
                                    .filter(|ival| {
                                        let pred_row = (src_diag, *ival).end_point_row();
                                        let succ_row = (diag, *ival).end_point_row();

                                        graph.is_successor(pred_row, succ_row)
                                    })
                                    .cloned_with_backtrace(Backtrace::new(src_diag, AlignState::MisMatch), true)
                                    .transform_deletion(diag, src_diag)
                                    .collect();

                                if !src_ivals.is_empty() { Some(src_ivals) } else { None }
                            })
                    })
                    .collect();

                let source_intervals_ref: SmallVec<[&IntervalSmallVec<O>; 16]> = source_ivals_open.iter()
                    .chain(source_ivals_del.iter())
                    .collect();

                eprintln!("Curr diag: {:?}", diag);
                VisitedIntervalTree::from_iter(merge_intervals(&source_intervals_ref))
            })
            .collect();

        eprintln!("DEL: {:?}", new_visited_del);
        eprintln!("Building MIS:");

        let new_visited_mismatch: VecDeque<_> = DiagRange::closed(new_kmin, new_kmax).into_iter().enumerate()
            .map(|(ix, diag) | {
                let mut source_ivals: SmallVec<[IntervalSmallVec<O>; 16]> = DiagRange::closed(new_kmin, new_kmax).into_iter()
                    // Iterate over all potential source diagonals
                    .filter_map(|src_diag| {
                        prev_layer_mis
                            .and_then(|layer| layer.layer_m.get_visited(src_diag))
                            .and_then(|itree| {
                                let src_ivals: IntervalSmallVec<O> = itree.iter()
                                    // Only include intervals with a successor that enters `diag`
                                    .filter(|ival| {
                                        let pred_row = (src_diag, *ival).end_point_row();
                                        let succ_row = (diag, ival.end()).row(); // Intervals are half-open, so `stop` is our actual new offset

                                        graph.is_successor(pred_row, succ_row)
                                            && ival.end.saturating_sub(O::one()) < seq_len_as_o
                                    })
                                    .cloned_with_backtrace(Backtrace::new(src_diag, AlignState::MisMatch), false)
                                    .transform_mismatch(diag, src_diag, seq_len_as_o + O::one())
                                    .collect();

                                if !src_ivals.is_empty() { Some(src_ivals) } else { None }
                            })
                    })
                    .collect();

                source_ivals.push(new_visited_del[ix].iter()
                    .cloned_with_backtrace(Backtrace::new(diag, AlignState::Deletion), false)
                    .collect()
                );
                source_ivals.push(new_visited_ins[ix].iter()
                    .cloned_with_backtrace(Backtrace::new(diag, AlignState::Insertion), false)
                    .collect()
                );

                eprintln!("Curr diag: {:?}", diag);
                let source_intervals_ref: SmallVec<[&IntervalSmallVec<O>; 8]> = source_ivals.iter().collect();
                VisitedIntervalTree::from_iter(merge_intervals(&source_intervals_ref))
            })
            .collect();

        eprintln!("MISMATCH: {:?}", new_visited_mismatch);

        let new_layer = GapAffineLayer::new(new_kmin, new_kmax,
                                            new_visited_mismatch, new_visited_del, new_visited_ins);

        if let Some(debug_output) = self.debug_output {
            debug_output.log(DebugOutputMessage::NewScore(new_score));

            debug_output.log(DebugOutputMessage::new_from_layer(AlignState::MisMatch, &new_layer.layer_m));
            debug_output.log(DebugOutputMessage::new_from_layer(AlignState::Deletion, &new_layer.layer_d));
            debug_output.log(DebugOutputMessage::new_from_layer(AlignState::Insertion, &new_layer.layer_i));
        }

        self.layers.push(new_layer)
    }

    fn backtrace<G: AlignableGraph>(&self, graph: &G, seq: &[u8], end_diag: Diag,
                                    end_interval: &VisitedIntervalType<O>) -> Alignment<G::NodeIndex> {
        let mut alignment = Alignment::with_capacity(seq.len());

        let mut curr_state = if end_interval.backtrace().prev_state == AlignState::Extended {
            AlignState::Extended
        } else {
            AlignState::MisMatch
        };

        let mut curr_diag = end_diag;
        let mut curr_score = self.layers.len() as i64 - 1;
        let mut curr_ival = end_interval;

        loop {
            let mut curr_row = (curr_diag, curr_ival).end_point_row::<G::NodeIndex>().as_usize();
            let mut curr_offset = (curr_ival.end - O::one()).as_usize();

            let bt = curr_ival.backtrace();
            let score_delta = self.score_delta(curr_state, bt.prev_state);
            eprintln!("Current interval: {:?} {:?} (state: {:?}; score: {:?}; delta: {:?})", curr_diag, curr_ival, curr_state, curr_score, score_delta);
            eprintln!("Current (row, offset): ({:?}, {:?})", curr_row, curr_offset);

            let search_range = if (curr_state == AlignState::Extended || curr_state == AlignState::MisMatch)
                    && bt.prev_state == AlignState::Extended {
                (curr_ival.start - O::one(), curr_ival.start)
            } else if curr_state == AlignState::Insertion {
                (curr_ival.start - O::one(), curr_ival.end - O::one())
            } else {
                (curr_ival.start, curr_ival.end)
            };

            let prev_ival = bt.prev_diag.and_then(|p_diag| {
                self.get_prev_layer(curr_score - score_delta)
                    .and_then(|layer|
                        layer.get_prev_ival(bt.prev_state, p_diag, search_range))
            });

            let offset_diff = if let Some(prev) = prev_ival {
                eprintln!("Previous interval: {:?}", prev);
                (curr_ival.end - prev.end).as_usize()
            } else {
                eprintln!("No previous interval (i.e. interval from start).");
                (curr_ival.end - curr_ival.start).as_usize() - 1
            };

            match curr_state {
                AlignState::Start | AlignState::MisMatch | AlignState::Extended => {
                    for _ in 0..offset_diff {
                        let node = graph.row_to_node_ix(curr_row);
                        alignment.push(AlignedPair { rpos: Some(node), qpos: Some(curr_offset - 1) });

                        eprintln!("{:?}, {:?} - {:?} (row: {:?}) ", alignment.last().unwrap(),
                            graph.get_symbol(node), char::from(seq[curr_offset - 1]), curr_row);

                        curr_row = curr_row.saturating_sub(1);
                        curr_offset = curr_offset.saturating_sub(1);
                    }
                },
                AlignState::Deletion => {
                    let node = graph.row_to_node_ix(curr_row);
                    alignment.push(AlignedPair { rpos: Some(node), qpos: None });

                    eprintln!("{:?}, {:?} - | (row: {:?}) ", alignment.last().unwrap(),
                              graph.get_symbol(node),  curr_row);
                }
                AlignState::Insertion => {
                    alignment.push(AlignedPair { rpos: None, qpos: Some(curr_offset - 1) });

                    eprintln!("{:?}, {:?} - | (row: {:?}) ", alignment.last().unwrap(),
                              char::from(seq[curr_offset - 1]),  curr_row);
                },
                AlignState::Deletion2 | AlignState::Insertion2 => panic!("Invalid gap-affine state!"),
            }

            if let Some(prev) = prev_ival {
                curr_ival = prev;
                curr_score -= score_delta;
                curr_state = bt.prev_state;
                curr_diag = bt.prev_diag.unwrap();
            } else {
                break;
            }

            eprintln!();
        }

        alignment.reverse();
        alignment
    }
}

impl<'a, O> Display for GapAffineCompute<'a, O>
where
    O: OffsetType,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        todo!()
    }

}
