use crate::aligner::aln_graph::AlignState;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentCosts;
use crate::graphs::NodeIndexType;

#[derive(Copy, Clone)]
pub struct GapLinear {
    cost_mismatch: u8,
    cost_gap: u8
}

impl AlignmentCosts for GapLinear {
    type AlignmentGraphType<'a, N, O> = LinearAlignmentGraph<'a, N, O>
    where
        N: NodeIndexType,
        O: OffsetType;

    #[inline]
    fn mismatch(&self) -> u8 {
        self.cost_mismatch
    }

    #[inline]
    fn gap_open(&self) -> u8 {
        0
    }

    #[inline]
    fn gap_extend(&self) -> u8 {
        self.cost_gap
    }

    #[inline]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline]
    fn gap_extend2(&self) -> u8 {
        0
    }

    #[inline]
    fn gap_cost(&self, _: AlignState, length: usize) -> usize {
        length * self.cost_gap as usize
    }
}

pub struct LinearAlignmentGraph<'a, N, O>;