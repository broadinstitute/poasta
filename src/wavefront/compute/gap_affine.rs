use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::alignment::{AlignedPair, Alignment};

use crate::graph::{NodeByRank, POAGraph};
use crate::wavefront::Wavefront;
use crate::wavefront::fr_points::{OffsetPrimitive, DiagIx, WFPoint, FRPoint, FRPointContainer, ExtendCandidate, State, PrevCandidate};
use super::WFCompute;



#[derive(Debug, Clone, Default)]
struct WavefrontSetGapAffine<Offset: OffsetPrimitive> {
    wavefront_m: Wavefront<Offset>,
    wavefront_i: Wavefront<Offset>,
    wavefront_d: Wavefront<Offset>,
}

impl<Offset: OffsetPrimitive> WavefrontSetGapAffine<Offset> {
    pub fn new() -> Self {
        Self {
            wavefront_m: Wavefront::initial(),
            wavefront_i: Wavefront::new(),
            wavefront_d: Wavefront::new()
        }
    }

    pub fn new_with_fr_points(k_lo: DiagIx, k_hi: DiagIx,
                              fr_points_m: FRPointContainer<Offset>,
                              fr_points_i: FRPointContainer<Offset>,
                              fr_points_d: FRPointContainer<Offset>
    ) -> Self {
        Self {
            wavefront_m: Wavefront::new_with_fr_points(k_lo, k_hi, fr_points_m),
            wavefront_i: Wavefront::new_with_fr_points(k_lo, k_hi, fr_points_i),
            wavefront_d: Wavefront::new_with_fr_points(k_lo, k_hi, fr_points_d)
        }
    }

    pub fn extend_candidates(&self) -> impl Iterator<Item=FRPoint<Offset>> + '_ {
        eprintln!("{:?}", self.wavefront_m);

        self.wavefront_m.iter()
    }

    pub fn extend_point(&mut self, curr_score: i64, candidate: &ExtendCandidate<Offset>) -> bool {
        self.wavefront_m.extend_candidate(curr_score, candidate)
    }

    pub fn reached_point(&self, point: &FRPoint<Offset>) -> bool {
        if let Some(offset) = self.wavefront_m.get(point.diag()) {
            offset >= point.offset()
        } else {
            false
        }
    }

    pub fn get_prev_point(&self, k: DiagIx, prev_state: &State) -> Option<&PrevCandidate<Offset>> {
        match prev_state {
            State::Start | State::Match(_) | State::Mismatch(_) =>
                self.wavefront_m.get_pointer(k),
            State::Deletion(_) =>
                self.wavefront_d.get_pointer(k),
            State::Insertion(_) =>
                self.wavefront_i.get_pointer(k)
        }
    }
}

pub struct WFComputeGapAffine<Offset: OffsetPrimitive> {
    pub cost_mismatch: i64,
    pub cost_gap_open: i64,
    pub cost_gap_extend: i64,

    wavefronts: Vec<WavefrontSetGapAffine<Offset>>
}

impl<Offset: OffsetPrimitive> WFComputeGapAffine<Offset> {
    fn get_prev_wf(&self, score: i64) -> Option<&WavefrontSetGapAffine<Offset>> {
        if score >= 0 && score < self.wavefronts.len() as i64 {
            Some(&self.wavefronts[score as usize])
        } else {
            None
        }
    }

    fn get_lowest_successor_diag(&self, graph: &POAGraph,
                                          wavefront: &Wavefront<Offset>) -> Option<DiagIx> {
        wavefront.iter()
            // For each FR point, identify the corresponding node in the graph, and find successor
            // with the highest rank difference.
            .filter_map(|point| {
                let max_rank_diff = graph.neighbors_for_rank(point.rank())
                    .map(|n| graph.get_node_rank(n))
                    .map(|succ_rank| succ_rank - point.rank())  // To rank difference
                    .max();

                // Difference in rank is the same as difference in diagonal index
                max_rank_diff.map(|v| point.diag() - v as DiagIx)
            })
            .min()
    }

    fn new_k_lo(&self, graph: &POAGraph, new_score: i64) -> DiagIx {
        // Take into account that the new k can be quite a bit lower because of a successor
        // with lower rank
        let k_lo_mis = self.get_prev_wf(new_score - self.cost_mismatch)
            .and_then(|prev_wf| self.get_lowest_successor_diag(graph, &prev_wf.wavefront_m));

        let k_lo_gap_open = self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend)
            .and_then(|prev_wf| self.get_lowest_successor_diag(graph, &prev_wf.wavefront_m));

        let k_lo_del_ext = self.get_prev_wf(new_score - self.cost_gap_extend)
            .and_then(|prev_wf| self.get_lowest_successor_diag(graph, &prev_wf.wavefront_d));

        let k_lo_ins_ext = self.get_prev_wf(new_score - self.cost_gap_extend)
            .map(|prev_wf| prev_wf.wavefront_i.k_lo - 1);

        let options = vec![k_lo_mis, k_lo_gap_open, k_lo_del_ext, k_lo_ins_ext];

        options.into_iter()
            .flatten()
            .min().unwrap_or(-1)
    }

    fn new_k_hi(&self, new_score: i64) -> DiagIx {
        // For the new k_hi, we don't need to take into account successors, because a successor node
        // will always have a lower rank, not higher, because of topological order.
        let k_hi_mis = self.get_prev_wf(new_score - self.cost_mismatch)
            .map(|prev_wf| prev_wf.wavefront_m.k_hi);

        let k_hi_gap_open = self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend)
            .map(|prev_wf| prev_wf.wavefront_m.k_hi);

        let k_hi_del_ext = self.get_prev_wf(new_score - self.cost_gap_extend)
            .map(|prev_wf| prev_wf.wavefront_d.k_hi);

        let k_hi_ins_ext = self.get_prev_wf(new_score - self.cost_gap_extend)
            .map(|prev_wf| prev_wf.wavefront_i.k_hi);

        let options = vec![k_hi_mis, k_hi_gap_open, k_hi_del_ext, k_hi_ins_ext];

        options.into_iter()
            .flatten()
            .max().unwrap_or(0) + 1
    }
}

impl<Offset: OffsetPrimitive> Default for WFComputeGapAffine<Offset> {
    fn default() -> Self {
        Self {
            cost_mismatch: 4,
            cost_gap_open: 6,
            cost_gap_extend: 2,

            wavefronts: vec![WavefrontSetGapAffine::new()]
        }
    }
}

impl<Offset: OffsetPrimitive> WFCompute for WFComputeGapAffine<Offset> {
    type OffsetType = Offset;

    fn reached_end(&self, graph: &POAGraph, seq_length: usize) -> Option<FRPoint<Offset>> {
        self.wavefronts.last()
            .and_then(|wf| {
                eprintln!("Last wavefront: {:?}", wf);
                graph.end_nodes().into_iter()
                    .map(|node| {
                        let rank = graph.get_node_rank(*node);
                        let diag = seq_length as DiagIx - rank as DiagIx;
                        eprintln!("endpoint: rank: {:?}, {:?}", rank, (diag, Offset::from(seq_length).unwrap()));

                        (diag, Offset::from(seq_length).unwrap())
                    })
                    .find(|p| {
                        wf.reached_point(p)
                    })
            })
    }

    fn reset(&mut self) {
        self.wavefronts.clear();
        self.wavefronts.push(WavefrontSetGapAffine::new());
    }

    fn extend_candidates(&self) -> Vec<FRPoint<Offset>> {
        match self.wavefronts.last() {
            Some(v) => v.extend_candidates().collect(),
            None => vec![]
        }
    }

    fn extend(&mut self, candidate: &ExtendCandidate<Offset>) -> bool {
        let score = self.wavefronts.len() as i64 - 1;
        if let Some(v) = self.wavefronts.last_mut() {
            v.extend_point(score, candidate)
        } else {
            false
        }
    }

    fn next(&mut self, graph: &POAGraph, seq_len: usize, new_score: i64) {
        let seqlen_as_offset = Offset::from(seq_len).unwrap();

        let k_lo = self.new_k_lo(graph, new_score);
        let k_hi = self.new_k_hi(new_score);
        eprintln!("new k_lo: {:?}, k_hi: {:?}", k_lo, k_hi);

        let new_fr_i: FRPointContainer<Offset> = (k_lo..=k_hi)
            .map(|k| {
                let values: Vec<PrevCandidate<Offset>> = vec![
                    self.get_prev_wf(new_score - self.cost_gap_extend)
                        .and_then(|prev_wf| prev_wf.wavefront_i.get(k-1))
                        .and_then(|offset|
                            if offset < seqlen_as_offset {
                                Some(PrevCandidate((k-1, offset), State::Insertion(new_score - self.cost_gap_extend)))
                            } else {
                                None
                            }
                        ),
                    self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend)
                        .and_then(|prev_wf| prev_wf.wavefront_m.get(k-1))
                        .and_then(|offset|
                            if offset < seqlen_as_offset {
                                Some(PrevCandidate((k-1, offset), State::Mismatch(new_score - self.cost_gap_open - self.cost_gap_extend)))
                            } else {
                                None
                            }
                        ),
                ].into_iter().flatten().collect();
                eprintln!("-I k: {:?}, prev options ({:?}, {:?}): {:?}", k,
                          new_score - self.cost_gap_open - self.cost_gap_extend,
                          new_score - self.cost_gap_extend, values);

                values.into_iter()
                    .max()
                    .map(|prev_max| {
                        let new_offset = prev_max.offset() + Offset::one();
                        prev_max.into_offset_with_bt(new_offset)
                    })
            }).collect();

        eprintln!("I_{}: {:?}", new_score, new_fr_i);
        eprintln!();

        let new_fr_d: FRPointContainer<Offset> = (k_lo..=k_hi)
            .map(|k| {
                let pred_candidates_d = (k+1..=k_hi)
                    // Get previous furthest reaching points from source wavefront
                    .filter_map(|k_prev| {
                        self.get_prev_wf(new_score - self.cost_gap_extend)
                            .and_then(|prev_wf| prev_wf.wavefront_d.get(k_prev))
                            .map(|offset|
                                PrevCandidate((k_prev, offset), State::Deletion(new_score - self.cost_gap_extend))
                            )
                    })
                    // Only include those that involve a valid edge
                    .filter(|prev_fr_point| {
                        let prev_rank = prev_fr_point.rank();
                        let new_rank = (k, prev_fr_point.offset()).rank();

                        if new_rank >= graph.max_rank() {
                            return false
                        }

                        graph.is_neighbor_rank(prev_rank, new_rank)
                    });

                let pred_candidates_m = (k+1..=k_hi)
                    // Get previous furthest reaching points from source wavefront
                    .filter_map(|k_prev| {
                        self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend)
                            .and_then(|prev_wf| prev_wf.wavefront_m.get(k_prev))
                            .map(|offset|
                                PrevCandidate((k_prev, offset), State::Mismatch(new_score - self.cost_gap_open - self.cost_gap_extend))
                            )
                    })
                    // Only include those that involve a valid edge
                    .filter(|prev_fr_point| {
                        let prev_rank = prev_fr_point.rank();
                        let new_rank = (k, prev_fr_point.offset()).rank();

                        if new_rank >= graph.max_rank() {
                            return false
                        }

                        graph.is_neighbor_rank(prev_rank, new_rank)
                    });

                pred_candidates_d.chain(pred_candidates_m)
                    .max()
                    .map(|prev_max| {
                        let new_offset = prev_max.offset();
                        prev_max.into_offset_with_bt(new_offset)
                    })
            }).collect();

        eprintln!("D_{}: {:?}", new_score, new_fr_d);
        eprintln!();

        let new_fr_m: FRPointContainer<Offset> = (k_lo..=k_hi)
            .map(|k| {
                // Previous was insertion/deletion state
                let prev_indel: Vec<PrevCandidate<Offset>> = vec![
                    new_fr_i[(k - k_lo) as usize].as_ref()
                        .map(|offset|
                                PrevCandidate((k, offset.offset), State::Insertion(new_score))),
                    new_fr_d[(k - k_lo) as usize].as_ref()
                        .map(|offset|
                            PrevCandidate((k, offset.offset), State::Deletion(new_score)))
                ].into_iter().flatten().collect();

                // Previous was (mis)match state, check all possible predecessors
                let prev_mis: Vec<PrevCandidate<Offset>> = (k..=k_hi)
                    // Get previous furthest reaching points from source wavefront
                    .filter_map(|k_prev| {
                        self.get_prev_wf(new_score - self.cost_mismatch)
                            .and_then(|prev_wf| prev_wf.wavefront_m.get(k_prev))
                            .map(|offset|
                                PrevCandidate((k_prev, offset), State::Mismatch(new_score - self.cost_mismatch)))
                    })
                    // Only include those that involve a valid edge
                    .filter(|prev_fr_point| {
                        let prev_rank = prev_fr_point.rank();
                        // Offset plus one because we are performing a mismatch
                        let new_point = (k, prev_fr_point.offset() + Offset::one());
                        let new_rank = new_point.rank();

                        if new_rank >= graph.max_rank() || new_point.offset() > seqlen_as_offset {
                            false
                        } else {
                            graph.is_neighbor_rank(prev_rank, new_rank)
                        }
                    })
                    .collect();


                eprintln!("- prev options (k: {:?}-{:?}): {:?}, indel: {:?}",
                    k, k_hi, prev_mis, prev_indel);

                // In case of ties, max() returns the last element of an iterator.
                // By listing previous candidates in order of insertion, deletion and (mis)match,
                // we ensure that (mis)match state gets priority over deletions, which in turn gets
                // priority over insertions
                prev_indel.into_iter()
                    .chain(prev_mis.into_iter())
                    .max()
                    .map(|prev_max| {
                        let new_offset = match prev_max.1 {
                            State::Match(_) | State::Mismatch(_) => prev_max.offset() + Offset::one(),
                            _ => prev_max.offset()
                        };
                        prev_max.into_offset_with_bt(new_offset)
                    })
            }).collect();
        eprintln!("M_{}: {:?}", new_score, new_fr_m);

        self.wavefronts.push(WavefrontSetGapAffine::new_with_fr_points(k_lo, k_hi, new_fr_m, new_fr_i, new_fr_d));
    }

    fn backtrace(&self, graph: &POAGraph, end_point: FRPoint<Offset>) -> Alignment {
        let mut alignment: Alignment = vec![];

        let mut curr_state = State::Match(self.wavefronts.len() as i64 - 1);
        let mut curr_offset = end_point.offset();
        let mut curr_diag = end_point.diag();

        loop {
            let curr_score = curr_state.score();
            let mut curr_rank = (curr_diag, curr_offset).rank();

            let curr_pointer = self.get_prev_wf(curr_score)
                .and_then(|wf| wf.get_prev_point(curr_diag, &curr_state))
                .unwrap();

            let PrevCandidate(prev_point, prev_state) = curr_pointer;

            eprintln!("Current (rank, offset): ({:?}, {:?}), diag: {:?}, state: {:?}", curr_rank, curr_offset, curr_diag, curr_state);
            eprintln!("Prev (rank, offset): ({:?}, {:?}), diag: {:?}, state: {:?})", prev_point.rank(), prev_point.offset(), prev_point.diag(), prev_state);

            let offset_diff: i64 = (curr_offset - prev_point.offset())
                .try_into()
                .unwrap_or_else(|_| panic!("Could not convert offset difference"));

            match curr_state {
                State::Start | State::Match(_) | State::Mismatch(_) => {
                    for _ in 0..offset_diff {
                        let offset: i64 = curr_offset.try_into()
                            .unwrap_or_else(|_| panic!("Could not convert current offset"));

                        if let NodeByRank::Node(node_ix) = graph.get_node_by_rank(curr_rank) {
                            alignment.push(AlignedPair::new(Some(node_ix), Some(offset as usize - 1)));
                            eprintln!("{:?} (rank: {:?})", alignment.last().unwrap(), curr_rank);
                        } else {
                            eprintln!("Ignoring start node.")
                        }

                        curr_offset = curr_offset - Offset::one();
                        curr_rank = (curr_diag, curr_offset).rank();
                    }
                },
                State::Deletion(_) => {
                    if let NodeByRank::Node(node_ix) = graph.get_node_by_rank(curr_rank) {
                        alignment.push(AlignedPair::new(Some(node_ix), None));
                        eprintln!("{:?} (rank: {:?})", alignment.last().unwrap(), curr_rank);
                    } else {
                        eprintln!("Ignoring start node.")
                    }
                }
                State::Insertion(_) => {
                    let offset: i64 = curr_offset.try_into()
                        .unwrap_or_else(|_| panic!("Could not convert current offset"));

                    alignment.push(AlignedPair::new(None, Some(offset as usize - 1)));
                    eprintln!("{:?} (rank: {:?})", alignment.last().unwrap(), curr_rank);
                    curr_offset = curr_offset - Offset::one();
                }
            }

            if *prev_state == State::Start {
                break;
            }

            eprintln!("- Setting new node rank to {:?} (old: {:?})", prev_point.rank(), curr_rank);
            curr_diag = prev_point.diag();
            curr_state = prev_state.clone();
            eprintln!()
        }

        alignment.into_iter().rev().collect()
    }

    fn write_csv(&self, writer: &mut BufWriter<File>) -> Result<(), Box<dyn Error>> {
        writeln!(writer, "score\tstate\tk\toffset\trank\tprev")?;

        for (score, wf) in self.wavefronts.iter().enumerate() {
            for (k, p) in (wf.wavefront_m.k_lo..=wf.wavefront_m.k_hi)
                .zip(wf.wavefront_m.iter())
            {
                writeln!(writer, "{}\t{}\t{}\t{:?}\t{}\t{:?}", score, "match", k, p.offset(), p.rank(),
                    wf.wavefront_m.get_pointer(k))?;
            }

            for (k, p) in (wf.wavefront_d.k_lo..=wf.wavefront_d.k_hi)
                .zip(wf.wavefront_d.iter())
            {
                writeln!(writer, "{}\t{}\t{}\t{:?}\t{}\t{:?}", score, "deletion", k, p.offset(), p.rank(),
                         wf.wavefront_d.get_pointer(k))?;
            }

            for (k, p) in (wf.wavefront_i.k_lo..=wf.wavefront_i.k_hi)
                .zip(wf.wavefront_i.iter())
            {
                writeln!(writer, "{}\t{}\t{}\t{:?}\t{}\t{:?}", score, "insertion", k, p.offset(), p.rank(),
                         wf.wavefront_i.get_pointer(k))?;
            }
        }

        Ok(())
    }
}
