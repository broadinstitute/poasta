use std::fmt::Debug;
use crate::debug::messages::DebugOutputMessage;
use crate::alignment::{AlignedPair, Alignment};
use crate::debug::DebugOutputWriter;

use crate::graph::{NodeByRank, POAGraph};
use crate::wavefront::GraphWavefront;
use crate::wavefront::offsets::{OffsetPrimitive, OffsetContainer, PrevCandidate, BacktraceState, OffsetCell, Backtrace};
use super::WFCompute;

/// A struct holding the query offsets when using gap-affine costs, thus holding offsets for
/// the (mis)match, insertion, and deletion states.
#[derive(Clone, Debug, Default)]
struct GraphWavefrontGapAffine<Offset: OffsetPrimitive> {
    wavefront_m: GraphWavefront<Offset>,
    wavefront_i: GraphWavefront<Offset>,
    wavefront_d: GraphWavefront<Offset>,
}

impl<Offset: OffsetPrimitive> GraphWavefrontGapAffine<Offset> {
    pub fn initial() -> Self {
        Self {
            wavefront_m: GraphWavefront::initial(),
            wavefront_i: GraphWavefront::new(),
            wavefront_d: GraphWavefront::new()
        }
    }

    pub fn new_with_offsets(offsets_m: OffsetContainer<Offset>, offsets_i: OffsetContainer<Offset>,
                            offsets_d: OffsetContainer<Offset>) -> Self {
        Self {
            wavefront_m: GraphWavefront::new_with_offsets(offsets_m),
            wavefront_i: GraphWavefront::new_with_offsets(offsets_i),
            wavefront_d: GraphWavefront::new_with_offsets(offsets_d),
        }
    }

    pub fn reached(&self, node_rank: usize, offset: Offset) -> bool {
        if let Some(node_offset) = self.wavefront_m.get(node_rank) {
            offset >= node_offset
        } else {
            false
        }
    }

    pub fn extend_candidates(&self) -> impl DoubleEndedIterator<Item=(usize, Offset)> + '_ {
        self.wavefront_m.iter_endpoints()
    }

    pub fn update_extended_path(&mut self, start_node: usize, path: Vec<(usize, usize)>) {
        let Some(max_rank) = path.last().map(|v| v.0) else {
            return
        };

        self.wavefront_m.resize(max_rank);

        let mut prev_node = start_node;
        let path_len = path.len();
        for (i, (node, offset)) in path.into_iter().enumerate() {
            let offset_conv = Offset::from(offset).unwrap();

            // Mark last node in path as an alignment end-point
            let new_cell = if i + 1 == path_len {
                OffsetCell::new_endpoint(offset_conv, Backtrace::new(prev_node, BacktraceState::Extend))
            } else {
                OffsetCell::new_extended(offset_conv, Backtrace::new(prev_node, BacktraceState::Extend))
            };
            self.wavefront_m.set(node, new_cell);

            prev_node = node
        }
    }

    pub fn get_prev_point(&self, node: usize, prev_state: &BacktraceState) -> Option<&Backtrace> {
        match prev_state {
            BacktraceState::End | BacktraceState::Mismatch | BacktraceState::Extend | BacktraceState::InsOpen | BacktraceState::DelOpen => {
                eprintln!("Getting previous point for node {:?} with prev state = M ({:?})", node, prev_state);
                eprintln!("M: {:?}", self.wavefront_m);
                self.wavefront_m.get_backtrace(node)
            }
            BacktraceState::InsExt | BacktraceState::InsClose => {
                eprintln!("Getting previous point for node {:?} with prev state = I ({:?})", node, prev_state);
                eprintln!("I: {:?}", self.wavefront_i);
                self.wavefront_i.get_backtrace(node)
            }
             BacktraceState::DelExt | BacktraceState::DelClose => {
                 eprintln!("Getting previous point for node {:?} with prev state = D ({:?})", node, prev_state);
                 eprintln!("D: {:?}", self.wavefront_d);
                 self.wavefront_d.get_backtrace(node)
             }
            _ => None
        }
    }
}

pub struct WFComputeGapAffine<Offset: OffsetPrimitive> {
    pub cost_mismatch: i64,
    pub cost_gap_open: i64,
    pub cost_gap_extend: i64,

    wavefronts: Vec<GraphWavefrontGapAffine<Offset>>
}


impl<Offset: OffsetPrimitive> WFComputeGapAffine<Offset> {
    fn get_prev_wf(&self, score: i64) -> Option<&GraphWavefrontGapAffine<Offset>> {
        if score >= 0 && score < self.wavefronts.len() as i64 {
            Some(&self.wavefronts[score as usize])
        } else {
            None
        }
    }

    /// For all nodes with a defined offset in a wavefront, find their successors, and identify the node with highest rank
    ///
    /// This is used to know which nodes we have assess in a new wavefront that depends on the given wavefront.
    ///
    /// TODO: can potentially preprocess the graph and store highest successor rank per node
    fn get_highest_successor_rank(&self, graph: &POAGraph, wavefront: &GraphWavefront<Offset>) -> Option<usize> {
        wavefront.iter_nodes()
            .filter_map(|node| {
                graph.successors(node).max()
            }).max()
    }

    fn score_delta(&self, state: &BacktraceState) -> i64 {
        match *state {
            BacktraceState::Start => 0,
            BacktraceState::End => 0,
            BacktraceState::Extend => 0,
            BacktraceState::Mismatch => self.cost_mismatch,
            BacktraceState::InsOpen => self.cost_gap_open + self.cost_gap_extend,
            BacktraceState::InsExt => self.cost_gap_extend,
            BacktraceState::InsClose => 0,
            BacktraceState::DelOpen => self.cost_gap_open + self.cost_gap_extend,
            BacktraceState::DelExt => self.cost_gap_extend,
            BacktraceState::DelClose => 0,
            _ => panic!("Invalid backtrace state for gap-affine scoring!")
        }
    }
}

impl<Offset: OffsetPrimitive> Default for WFComputeGapAffine<Offset> {
    fn default() -> Self {
        Self {
            cost_mismatch: 4,
            cost_gap_open: 6,
            cost_gap_extend: 2,

            wavefronts: vec![GraphWavefrontGapAffine::initial()]
        }
    }
}

impl<Offset: OffsetPrimitive> WFCompute for WFComputeGapAffine<Offset> {
    type OffsetType = Offset;

    fn furthest_offsets(&self) -> Option<&OffsetContainer<Self::OffsetType>> {
        self.wavefronts.last()
            .and_then(|wf| wf.wavefront_m.get_offsets())
    }

    fn reached_end(&self, graph: &POAGraph, seq_length: usize) -> Option<usize> {
        let seq_len_as_offset = Offset::from(seq_length).unwrap();

        self.wavefronts.last()
            .and_then(|wf| {
                eprintln!("Last wavefront: {:?}", wf);
                graph.end_nodes().iter()
                    .map(|v| {
                        eprintln!("Looking for end point {:?}", (graph.get_node_rank(*v), seq_length));
                        graph.get_node_rank(*v)
                    })
                    .find(|node_rank|
                        wf.reached(*node_rank, seq_len_as_offset))
            })
    }

    fn reset(&mut self) {
        self.wavefronts.clear();
        self.wavefronts.push(GraphWavefrontGapAffine::initial());
    }

    fn extend_candidates(&self) -> Vec<(usize, Offset)> {
        match self.wavefronts.last() {
            // Reverse to extend candidates with higher node rank first
            Some(v) => v.extend_candidates().rev().collect(),
            None => vec![]
        }
    }

    fn is_further(&self, node: usize, offset: usize) -> bool {
        let offset_conv = Offset::from(offset).unwrap();
        match self.wavefronts.last() {
            Some(wf) => wf.wavefront_m.get(node).map_or(true, |furthest| furthest <= offset_conv),
            None => true
        }
    }

    fn update_extended_path(&mut self, start_node: usize, path: Vec<(usize, usize)>) {
        if let Some(wf) = self.wavefronts.last_mut() {
            wf.update_extended_path(start_node, path)
        }
    }

    fn next(&mut self, graph: &POAGraph, seq_len: usize, new_score: i64) {
        let seqlen_as_offset = Offset::from(seq_len).unwrap();

        // For which nodes can we open or extend an insertion?
        let prev_wf_ext = self.get_prev_wf(new_score - self.cost_gap_extend);
        let prev_wf_open = self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend);
        let max_ins_rank = vec![
            prev_wf_ext.and_then(|prev_wf| prev_wf.wavefront_i.max_rank()),
            prev_wf_open.and_then(|prev_wf| prev_wf.wavefront_m.max_rank())
        ].into_iter().flatten().max();

        let new_offsets_i: OffsetContainer<Offset> = if let Some(max_node) = max_ins_rank {
            (0..=max_node).map(|node| {
                vec![
                    prev_wf_open
                        .and_then(|prev_wf| prev_wf.wavefront_m.get_if_endpoint(node)
                            .map(|v| PrevCandidate::new(v, Backtrace::new(node, BacktraceState::InsOpen)))),
                    prev_wf_ext
                        .and_then(|prev_wf| prev_wf.wavefront_i.get(node)
                            .map(|v| PrevCandidate::new(v, Backtrace::new(node, BacktraceState::InsExt)))),
                ].into_iter()
                    .flatten()
                    .max()
                    .map(|candidate| {
                        let offset = candidate.offset() + Offset::one();
                        candidate.into_offset_with_bt(offset)
                    })
            })
            .collect()
        } else {
            OffsetContainer::default()
        };

        eprintln!("New I (max rank: {:?}): {:?}", max_ins_rank, new_offsets_i);

        // Find highest successor rank of source wavefronts which could end up in a deletion state
        let max_del_rank = vec![
            prev_wf_ext.and_then(|prev_wf| self.get_highest_successor_rank(graph, &prev_wf.wavefront_d)),
            prev_wf_open.and_then(|prev_wf| self.get_highest_successor_rank(graph, &prev_wf.wavefront_m)),
        ].into_iter().flatten().max();

        let new_offsets_d: OffsetContainer<Offset> = if let Some(max_node) = max_del_rank {
            (0..=max_node).map(|node| {
                graph.predecessors(node)
                    .filter_map(|pred| {
                        vec![
                            prev_wf_open.and_then(|prev_wf| prev_wf.wavefront_m.get_if_endpoint(pred)
                                .map(|offset| PrevCandidate::new(offset, Backtrace::new(pred, BacktraceState::DelOpen)))),
                            prev_wf_ext.and_then(|prev_wf| prev_wf.wavefront_d.get(pred)
                                .map(|offset| PrevCandidate::new(offset, Backtrace::new(pred, BacktraceState::DelExt)))),
                        ].into_iter().flatten().max()
                    })
                    .max()
                    .map(|max_prev| {
                        let offset = max_prev.offset();
                        max_prev.into_offset_with_bt(offset)
                    })
            })
            .collect()
        } else {
            OffsetContainer::default()
        };

        eprintln!("New D (max rank: {:?}): {:?}", max_del_rank, new_offsets_d);

        // Find highest successor rank that could end up in a (mis)match state
        let prev_wf_mis = self.get_prev_wf(new_score - self.cost_mismatch);
        let max_mm_rank = vec![
            max_ins_rank,
            max_del_rank,
            prev_wf_mis.and_then(|prev_wf| self.get_highest_successor_rank(graph, &prev_wf.wavefront_m))
        ].into_iter().flatten().max();

        let new_offsets_m: OffsetContainer<Offset> = if let Some(max_node) = max_mm_rank {
            (0..=max_node).map(|node| {
                // Candidates where the previous state was insertion or deletion
                vec![
                    new_offsets_i.get(node).and_then(|v| v.as_ref()).map(|v|
                        PrevCandidate::new(v.offset(), Backtrace::new(node, BacktraceState::InsClose))),
                    new_offsets_d.get(node).and_then(|v| v.as_ref()).map(|v|
                        PrevCandidate::new(v.offset(), Backtrace::new(node, BacktraceState::DelClose)))
                ].into_iter()
                    // Combine with predecessors in (mis)match state
                    .chain(graph.predecessors(node)
                        .map(|pred| {
                            prev_wf_mis.and_then(|prev_wf| prev_wf.wavefront_m.get_if_endpoint(pred)
                                .map(|offset|
                                    PrevCandidate::new(offset + Offset::one(), Backtrace::new(pred, BacktraceState::Mismatch))))
                        })
                    )
                    .flatten()
                    .max()
                    .and_then(|max_prev| {
                        // Prevent from going beyond the query length
                        if max_prev.offset() >= seqlen_as_offset {
                            None
                        } else {
                            let offset = max_prev.offset();
                            Some(max_prev.into_offset_with_bt(offset))
                        }
                    })
            })
            .collect()
        } else {
            OffsetContainer::default()
        };

        eprintln!("New M (max rank: {:?}): {:?}", max_mm_rank, new_offsets_m);

        self.wavefronts.push(GraphWavefrontGapAffine::new_with_offsets(
            new_offsets_m, new_offsets_i, new_offsets_d));
    }

    fn backtrace(&self, graph: &POAGraph, seq: &[u8], end_node: usize) -> Alignment {
        let mut alignment: Alignment = vec![];

        let mut curr_score = (self.wavefronts.len() - 1) as i64;
        let mut curr_state = Backtrace::new(end_node, BacktraceState::End);
        let mut curr_offset = seq.len();
        let mut curr_node = end_node;

        loop {
            eprintln!("Curr state: {:?}, score: {:?}, node: {:?}, qry offset: {:?}", curr_state.state(), curr_score, curr_node, curr_offset);

            match curr_state.state() {
                BacktraceState::End => curr_offset -= 1,
                BacktraceState::Mismatch | BacktraceState::Extend | BacktraceState::InsOpen | BacktraceState::DelOpen => {
                    if let NodeByRank::Node(node_ix) = graph.get_node_by_rank(curr_node) {
                        alignment.push(AlignedPair::new(Some(node_ix), Some(curr_offset)));
                        eprintln!("{:?} (rank: {:?}, symbols: {}-{})", alignment.last().unwrap(), curr_node,
                                  char::from(graph.graph[node_ix].symbol), char::from(seq[curr_offset]));
                    } else {
                        eprintln!("Ignoring start node.")
                    }

                    curr_offset -= 1;
                },
                BacktraceState::DelExt | BacktraceState::DelClose => {
                    if let NodeByRank::Node(node_ix) = graph.get_node_by_rank(curr_node) {
                        alignment.push(AlignedPair::new(Some(node_ix), None));
                        eprintln!("{:?} (rank: {:?}, symbol: {})", alignment.last().unwrap(), curr_node, char::from(graph.graph[node_ix].symbol));
                    } else {
                        eprintln!("Ignoring start node.")
                    }
                }
                 BacktraceState::InsExt | BacktraceState::InsClose => {
                    alignment.push(AlignedPair::new(None, Some(curr_offset)));
                    eprintln!("{:?} (rank: {:?}, insertion: {})", alignment.last().unwrap(), curr_node, char::from(seq[curr_offset]));
                    curr_offset -= 1;
                },
                _ => ()
            }

            let backtrace = self.get_prev_wf(curr_score)
                .and_then(|wf| wf.get_prev_point(curr_node, curr_state.state()))
                .unwrap();

            if *backtrace.state() == BacktraceState::Start {
                break;
            }

            eprintln!("{:?}, thus score reduced with {}-{}={}", curr_state.state(), curr_score,
                      self.score_delta(curr_state.state()),
                      curr_score - self.score_delta(curr_state.state()));

            curr_score -= self.score_delta(curr_state.state());
            curr_node = curr_state.prev_node();
            curr_state = backtrace.clone();

            eprintln!()
        }

        alignment.into_iter().rev().collect()
    }

    fn log_debug_data(&self, debug: &DebugOutputWriter) {
        if let Some(wf) = self.wavefronts.last() {
            let score = self.wavefronts.len() as i64 - 1;
            debug.log(DebugOutputMessage::new_from_wavefront(score, "match", &wf.wavefront_m));
            debug.log(DebugOutputMessage::new_from_wavefront(score, "insertion", &wf.wavefront_i));
            debug.log(DebugOutputMessage::new_from_wavefront(score, "deletion", &wf.wavefront_d));
        }
    }
}
