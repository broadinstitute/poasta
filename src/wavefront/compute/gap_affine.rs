use std::convert::identity;
use std::fmt::Debug;

use crate::graph::POAGraph;
use crate::wavefront::{OffsetPrimitive, DiagIx, WFPoint, FRPoint, Wavefront, FRPointContainer};
use super::WFCompute;

#[derive(Debug, Clone, Default)]
struct WavefrontSetGapAffine<Offset: OffsetPrimitive> {
    wavefront_m: Wavefront<Offset>,
    wavefront_i: Wavefront<Offset>,
    wavefront_d: Wavefront<Offset>,
}

impl<Offset: OffsetPrimitive> WavefrontSetGapAffine<Offset> {
    pub fn new(k_lo: DiagIx, k_hi: DiagIx,
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

    pub fn initial() -> Self {
        Self {
            wavefront_m: Wavefront::new_with_fr_points(0, 0, vec![Some(Offset::zero())].into()),
            wavefront_i: Wavefront::new(),
            wavefront_d: Wavefront::new()
        }
    }

    pub fn k_lo(&self) -> i64 {
        let ks = vec![self.wavefront_m.k_lo, self.wavefront_i.k_lo, self.wavefront_d.k_lo];

        ks.into_iter().min().unwrap_or(0)
    }

    pub fn k_hi(&self) -> i64 {
        let ks = vec![self.wavefront_m.k_hi, self.wavefront_i.k_hi, self.wavefront_d.k_hi];

        ks.into_iter().max().unwrap_or(0)
    }

    pub fn extend_candidates(&self) -> impl Iterator<Item=FRPoint<Offset>> + '_ {
        self.wavefront_m.iter()
    }

    pub fn extend_point(&mut self, point: &FRPoint<Offset>) {
        eprintln!("Extending point {:?} to increase offset with one", point);
        self.wavefront_m.update_fr_point((point.diag(), point.offset() + Offset::one()));
    }

    pub fn reached_point(&self, point: &FRPoint<Offset>) -> bool {
        if let Some(p) = self.wavefront_m.get_fr_point(point.diag()) {
            p.offset() == point.offset()
        } else {
            false
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
            eprintln!("- Prev score {:?}, wf {:?}", score, &self.wavefronts[score as usize]);
            Some(&self.wavefronts[score as usize])
        } else {
            None
        }
    }

    fn get_lowest_successor_diag<Seq: Eq + Clone>(&self, graph: &POAGraph<Seq>,
                                          wavefront: &Wavefront<Offset>) -> Option<DiagIx> {
        wavefront.iter()
            // For each FR point, identify the corresponding node in the graph, and find successor
            // with the highest rank difference.
            .map(|point| {
                let node = graph.get_node_by_rank(point.rank());
                let max_rank_diff = graph.graph.neighbors(node)
                    .map(|n| graph.get_node_rank(n))
                    .map(|succ_rank| succ_rank - point.rank())  // To rank difference
                    .max();

                // Difference in rank is the same as difference in diagonal index
                max_rank_diff.map(|v| point.diag() - v as DiagIx)
            })
            .filter_map(identity)
            .min()
    }

    fn new_k_lo<Seq: Eq + Clone>(&self, graph: &POAGraph<Seq>, new_score: i64) -> DiagIx {
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

        options.into_iter().
            filter_map(identity).
            min().unwrap_or(-1)
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

        options.into_iter().
            filter_map(identity).
            max().unwrap_or(0) + 1
    }
}

impl<Offset: OffsetPrimitive> Default for WFComputeGapAffine<Offset> {
    fn default() -> Self {
        Self {
            cost_mismatch: 4,
            cost_gap_open: 6,
            cost_gap_extend: 2,

            wavefronts: vec![WavefrontSetGapAffine::initial()]
        }
    }
}

impl<Offset: OffsetPrimitive> WFCompute<Offset> for WFComputeGapAffine<Offset> {
    fn reached_point(&self, point: &FRPoint<Offset>) -> bool {
        match self.wavefronts.last() {
            Some(v) => v.reached_point(point),
            None => false
        }
    }

    fn extend_candidates(&self) -> Vec<FRPoint<Offset>> {
        match self.wavefronts.last() {
            Some(v) => v.extend_candidates().collect(),
            None => vec![]
        }
    }

    fn extend(&mut self, point: &FRPoint<Offset>) {
        if let Some(v) = self.wavefronts.last_mut() {
            v.extend_point(point);
        }
    }

    fn next<Seq: Eq + Clone>(&mut self, graph: &POAGraph<Seq>, new_score: i64) {
        let k_lo = self.new_k_lo(graph, new_score);
        let k_hi = self.new_k_hi(new_score);
        eprintln!("new k_lo: {:?}, k_hi: {:?}", k_lo, k_hi);

        let new_fr_i: FRPointContainer<Offset> = (k_lo..=k_hi)
            .map(|k| {
                let values: Vec<Offset> = vec![
                    self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend)
                        .and_then(|prev_wf| prev_wf.wavefront_m.get(k-1)),
                    self.get_prev_wf(new_score - self.cost_gap_extend)
                        .and_then(|prev_wf| prev_wf.wavefront_i.get(k-1)),
                ].into_iter().filter_map(identity).collect();
                eprintln!("-I k: {:?}, prev options ({:?}, {:?}): {:?}", k,
                          new_score - self.cost_gap_open - self.cost_gap_extend,
                          new_score - self.cost_gap_extend, values);

                values.into_iter()
                    .max()
                    .map(|max_offset| max_offset + Offset::one())
            }).collect();

        eprintln!("I_{}: {:?}", new_score, new_fr_i);
        eprintln!();

        let new_fr_d: FRPointContainer<Offset> = (k_lo..=k_hi)
            .map(|k| {
                let pred_candidates_m = (k+1..=k_hi)
                    // Get previous furthest reaching points from source wavefront
                    .filter_map(|k_prev| {
                        self.get_prev_wf(new_score - self.cost_gap_open - self.cost_gap_extend)
                            .and_then(|prev_wf| prev_wf.wavefront_m.get(k_prev))
                            .map(|offset| (k_prev, offset.clone()))
                    })
                    // Only include those that involve a valid edge
                    .filter(|prev_fr_point| {
                        let prev_rank = prev_fr_point.rank();
                        let new_rank = (k, prev_fr_point.offset()).rank();

                        if new_rank >= graph.graph.node_count() {
                            return false
                        }

                        let node1 = graph.get_node_by_rank(prev_rank);
                        let node2 = graph.get_node_by_rank(new_rank);

                        let is_neighbor = graph.is_neighbor(node1, node2);
                        eprintln!(" - {:?} -> {:?} -- is true edge: {:?}",
                                  prev_rank, new_rank, is_neighbor);

                        is_neighbor
                    })
                    .map(|fr_point| fr_point.offset());

                let pred_candidates_d = (k+1..=k_hi)
                    // Get previous furthest reaching points from source wavefront
                    .filter_map(|k_prev| {
                        self.get_prev_wf(new_score - self.cost_gap_extend)
                            .and_then(|prev_wf| prev_wf.wavefront_d.get(k_prev))
                            .map(|offset| (k_prev, offset.clone()))
                    })
                    // Only include those that involve a valid edge
                    .filter(|prev_fr_point| {
                        let prev_rank = prev_fr_point.rank();
                        let new_rank = (k, prev_fr_point.offset()).rank();

                        if new_rank >= graph.graph.node_count() {
                            return false
                        }

                        let node1 = graph.get_node_by_rank(prev_rank);
                        let node2 = graph.get_node_by_rank(new_rank);

                        let is_neighbor = graph.is_neighbor(node1, node2);
                        eprintln!("- {:?} -> {:?} -- is true edge: {:?}",
                                  prev_rank, new_rank, is_neighbor);

                        is_neighbor
                    })
                    .map(|fr_point| fr_point.offset());

                pred_candidates_m.chain(pred_candidates_d).max()
            }).collect();

        eprintln!("D_{}: {:?}", new_score, new_fr_d);
        eprintln!();

        let new_fr_m: FRPointContainer<Offset> = (k_lo..=k_hi)
            .map(|k| {
                let prev_mis: Vec<Offset> = (k..=k_hi)
                    // Get previous furthest reaching points from source wavefront
                    .filter_map(|k_prev| {
                        eprintln!("k': {:?}", k_prev);
                        self.get_prev_wf(new_score - self.cost_mismatch)
                            .and_then(|prev_wf| prev_wf.wavefront_m.get(k_prev))
                            .map(|offset| (k_prev, offset.clone()))
                    })
                    // Only include those that involve a valid edge
                    .filter(|prev_fr_point| {
                        let prev_rank = prev_fr_point.rank();
                        // Offset plus one because we are performing a (mis)match
                        let new_rank = (k, prev_fr_point.offset() + Offset::one()).rank();

                        if new_rank >= graph.graph.node_count() {
                            return false
                        }

                        let node1 = graph.get_node_by_rank(prev_rank);
                        let node2 = graph.get_node_by_rank(new_rank);

                        let is_neighbor = graph.is_neighbor(node1, node2);
                        eprintln!("- {:?} -> {:?} -- is true edge: {:?}",
                                 prev_rank, new_rank, is_neighbor);

                        is_neighbor
                    })
                    .map(|fr_point| fr_point.offset())
                    .collect();

                // Previous insertion/deletion state
                let prev_indel: Vec<Offset> = vec![
                    new_fr_i[(k - k_lo) as usize],
                    new_fr_d[(k - k_lo) as usize]
                ].into_iter().filter_map(identity).collect();

                eprintln!("- prev options (k: {:?}-{:?}): {:?}, indel: {:?}",
                    k, k_hi, prev_mis, prev_indel);

                prev_mis.into_iter()
                    .chain(prev_indel.into_iter())
                    .max().map(|v| v + Offset::one())
            }).collect();
        eprintln!("M_{}: {:?}", new_score, new_fr_m);

        self.wavefronts.push(WavefrontSetGapAffine::new(k_lo, k_hi, new_fr_m, new_fr_i, new_fr_d));
    }
}
