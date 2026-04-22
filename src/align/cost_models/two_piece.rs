use super::AlignmentCostModel;
use crate::align::kernels::two_piece::TwoPieceAffineKernel;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TwoPieceAffine {
    equal: u8,
    mismatch: u8,
    gap_open: u8,
    gap_extend: u8,
    gap_open2: u8,
    gap_extend2: u8,
}

impl TwoPieceAffine {
    pub fn new(
        equal: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
        gap_open2: u8,
        gap_extend2: u8,
    ) -> Self {
        Self { equal, mismatch, gap_open, gap_extend, gap_open2, gap_extend2 }
    }
}

impl AlignmentCostModel for TwoPieceAffine {
    type Kernel = TwoPieceAffineKernel;
    #[inline(always)] fn equal(&self) -> u8 { self.equal }
    #[inline(always)] fn mismatch(&self) -> u8 { self.mismatch }
    #[inline(always)] fn gap_open(&self) -> u8 { self.gap_open }
    #[inline(always)] fn gap_extend(&self) -> u8 { self.gap_extend }
    #[inline(always)] fn gap_open2(&self) -> u8 { self.gap_open2 }
    #[inline(always)] fn gap_extend2(&self) -> u8 { self.gap_extend2 }
    #[inline(always)] fn min_gap_extend(&self) -> u8 { self.gap_extend.min(self.gap_extend2) }
}
