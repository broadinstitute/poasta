extern crate num;
use num::{FromPrimitive, PrimInt};

use std::collections::VecDeque;
use std::cmp;
use std::convert::identity;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

use crate::graph;
use graph::POAGraph;

pub mod compute;


const SCORE_MISMATCH: u32 = 2;
const SCORE_GAP: u32 = 4;

pub type DiagIx = i64;

trait OffsetPrimitive: FromPrimitive + PrimInt
    + Into<DiagIx>
    + Debug + Default + Copy { }

pub type FRPoint<Offset: OffsetPrimitive> = (DiagIx, Offset);
pub type FRPointContainer<Offset: OffsetPrimitive> = VecDeque<Option<Offset>>;

trait WFPoint<Offset: OffsetPrimitive>: Debug {
    fn rank(&self) -> usize;
    fn offset(&self) -> Offset;
    fn diag(&self) -> DiagIx;
}

impl<Offset: OffsetPrimitive> WFPoint<Offset> for FRPoint<Offset> {
    fn rank(&self) -> usize {
        (self.1.into() - self.0) as usize
    }

    fn offset(&self) -> Offset {
        self.1
    }

    fn diag(&self) -> DiagIx {
        self.0
    }
}

#[derive(Debug, Clone)]
struct Wavefront<Offset: OffsetPrimitive> {
    pub k_lo: DiagIx,
    pub k_hi: DiagIx,

    furthest_points: FRPointContainer<Offset>,
}

impl<Offset: OffsetPrimitive> Wavefront<Offset> {
    pub fn new() -> Self {
        Wavefront::default()
    }

    pub fn new_with_size(k_lo: DiagIx, k_hi: DiagIx) -> Self {
        let mut wf = Self::default();
        let len = ((k_hi - k_lo) + 1) as usize;
        wf.furthest_points.resize(len, None);

        wf
    }

    fn diag_to_ix(&self, diag: DiagIx) -> Option<usize> {
        if diag >= self.k_lo && diag <= self.k_hi {
            Some((diag - self.k_lo) as usize)
        } else {
            None
        }
    }

    fn ensure_size(&mut self, diag: DiagIx) {
        let diff = if diag < self.k_lo {
            (self.k_lo + diag) as usize
        } else if diag > self.k_hi {
            (diag - self.k_hi) as usize
        } else {
            0
        };

        if diff != 0 {
            self.furthest_points.reserve(self.furthest_points.len() + diff);
            if diag < self.k_lo {
                for i in 0..diff {
                    self.furthest_points.push_front(None);
                }
            } else {
                for i in 0..diff {
                    self.furthest_points.push_back(None);
                }
            }
        }
    }

    pub fn update_fr_point(&mut self, point: &FRPoint<Offset>) {
        self.ensure_size(point.diag());
        let ix = self.diag_to_ix(point.diag()).unwrap();

        match &self.furthest_points[ix] {
            Some(ref mut v) => *v = cmp::max(*v, point.offset()),
            None => self.furthest_points[ix as usize] = Some(point.offset())
        };

        self.k_lo = cmp::min(point.diag(), self.k_lo);
        self.k_hi = cmp::max(point.diag(), self.k_hi);
    }

    pub fn get_fr_point(&self, k: DiagIx) -> Option<FRPoint<Offset>> {
        let ix = self.diag_to_ix(k)?;

        return if let Some(offset) = &self.furthest_points[ix] {
            Some((k, *offset))
        } else {
            None
        }
    }

    pub fn iter(&self) -> impl Iterator<Item=FRPoint<Offset>> + '_ {
        self.furthest_points.iter()
            .zip(self.k_lo..=self.k_hi)
            .filter_map(|v| match v {
                (Some(offset), diag) => Some((diag, offset.clone())),
                _ => None
            })
    }
}

impl<Offset: OffsetPrimitive> Default for Wavefront<Offset> {
    fn default() -> Self {
        Wavefront {
            k_lo: 0,
            k_hi: 0,
            furthest_points: vec![Some(Offset::default())].into(),
        }
    }
}

impl<Offset: OffsetPrimitive> Index<DiagIx> for Wavefront<Offset> {
    type Output = Option<Offset>;
    fn index(&self, index: DiagIx) -> &Self::Output {
        assert!(index >= self.k_lo && index <= self.k_hi);

        &self.furthest_points[(index - self.k_lo) as usize]
    }
}

impl<Offset: OffsetPrimitive> IndexMut<DiagIx> for Wavefront<Offset> {
    fn index_mut(&mut self, index: DiagIx) -> &mut Self::Output {
        assert!(index >= self.k_lo && index <= self.k_hi);

        &mut self.furthest_points[(index - self.k_lo) as usize]
    }
}


#[derive(Debug, Clone)]
struct WFGapAffine {
    pub penalty_mismatch: usize,
    pub penalty_gap_open: usize,
    pub penalty_gap_extend: usize
}

impl Default for WFGapAffine {
    fn default() -> Self {
        WFGapAffine { penalty_mismatch: 4, penalty_gap_open: 8, penalty_gap_extend: 6 }
    }
}

impl WavefrontPenalties for WFGapAffine {
    fn smallest_score_step(&self) -> usize {
        let values = vec![self.penalty_mismatch, self.penalty_gap_open, self.penalty_gap_extend];

        values.into_iter().min().unwrap()
    }
}

trait WFCompute<Offset: OffsetPrimitive> {
    fn reached_point(&self, point: &FRPoint<Offset>) -> bool;

    fn extend_candidates(&self) -> Vec<FRPoint<Offset>>;
    fn extend(&mut self, point: &FRPoint<Offset>);
    fn next<Seq: Eq>(&mut self, graph: POAGraph<Seq>, new_score: usize);
}

struct WavefrontPOAligner<'a, A, Offset, P>
where
    Offset: OffsetPrimitive,
    P: WavefrontPenalties
{
    graph: &'a POAGraph<A>,
    method: P,
    wavefronts: Vec<WavefrontSetGapAffine<Offset>>
}

impl<'a, A, Offset, P> WavefrontPOAligner<'a, A, Offset, P>
where
    A: Eq,
    Offset: OffsetPrimitive,
    P: WavefrontPenalties,
{
    pub fn new(graph: &'a POAGraph<A>) -> Self {
        WavefrontPOAligner {
            graph,
            compute: P::default(),
            wavefronts: vec![WavefrontSetGapAffine::<Offset>::default()]
        }
    }

    pub fn align(&mut self, seq: &[A]) {
        let k_end = seq.len() as DiagIx - self.graph.graph.node_count() as i64;
        let offset_end = Offset::from(seq.len() - 1).unwrap();
        let fr_end = (k_end, offset_end) as FRPoint<Offset>;

        let mut score = 0;
        loop {
            self.wf_extend(self.wavefronts.last_mut().unwrap(), seq);

            if self.wavefronts[score].reached_point(&fr_end) {
                println!("Reached the end!");
                break;
            }

            score += self.penalties.smallest_score_step();

            let new_wavefront = self.wavefronts[score-1].wf_next(self, seq, score);

        }
    }

    pub fn get_source_wf(&self, score_diff: usize) -> Option<&WavefrontSetGapAffine<Offset>> {
        let steps = score_diff / self.penalties.smallest_score_step();
        debug_assert!(score_diff % self.penalties.smallest_score_step() == 0);

        let ix = self.wavefronts.len() - steps - 1;

        if ix >= 0 {
            Some(&self.wavefronts[ix])
        } else {
            None
        }
    }


    fn wf_extend(&mut self, wf: &mut WavefrontSetGapAffine<Offset>, seq: &[A]) {
        // Check if we can extend without mismatches along our current diagonals. Follow paths
        // in the graph in a depth first manner, and update furthest reaching points accordingly.
        let mut stack: Vec<FRPoint<Offset>> = Vec::new();

        for point in wf.extend_candidates() {
            stack.push(point)
        }

        println!("Stack at start: {:?}", stack);

        while let Some(point) = stack.pop() {
            let rank = point.rank();
            let node = self.graph.get_node_by_rank(rank);

            println!("Popped item {:?}, node_rank: {}", point, rank);
            println!("Remaining stack: {:?}", stack);

            // Sequence matches, add successors to stack
            if self.graph.graph[node].code == seq[point.offset().into() as usize] {
                wf.extend_point(&point);

                let mut num_neighbors = 0;
                for neighbor in self.graph.graph.neighbors(node) {
                    // `succ_rank` is always greater than `rank` because of topological order
                    let succ_rank = self.graph.get_node_rank(neighbor);
                    let new_k = point.diag() - (succ_rank - rank - 1) as DiagIx;

                    // Add neighbor with updated diagonal and one step further along `seq`
                    stack.push((new_k, point.offset() + Offset::one()));
                    num_neighbors += 1;
                }
            }

            println!("{:?}", wf);
        }
    }

    fn wf_next(&mut self, current_score: usize) {
        // First, determine the new k_lo and k_hi
        let mut curr_wf = self.wavefronts.last().unwrap();
        let k_lo = self.new_k_lo(curr_wf);
        let k_hi = self.new_k_hi(curr_wf);



    }

    fn new_k_lo(&self, wf: &WavefrontSetGapAffine<Offset>) -> DiagIx {
        wf.all_fr_points()
            .map(|point| {
                let node = self.graph.get_node_by_rank(point.rank());
                let max_rank_diff_opt = self.graph.graph.neighbors(node)
                    .map(|n| self.graph.get_node_rank(n))
                    .map(|succ_rank| succ_rank - point.rank())  // To rank difference
                    .max();

                // Difference in rank is the same as difference in diagonal index
                return if let Some(max_rank_diff) = max_rank_diff_opt {
                    Some(point.diag() - max_rank_diff as DiagIx)
                } else {
                    None
                }
            })
            .filter_map(identity)
            .min()
            .unwrap_or(wf.k_lo() - 1)
    }

    fn new_k_hi(&self, wf: &WavefrontSetGapAffine<Offset>) -> DiagIx {
        wf.k_hi() + 1
    }

}
