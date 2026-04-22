use super::AlignmentCostModel;
use crate::align::kernels::affine::AffineKernel;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Affine {
    equal: u8,
    mismatch: u8,
    gap_open: u8,
    gap_extend: u8,
}

impl Affine {
    pub fn new(match_: u8, mismatch: u8, gap_open: u8, gap_extend: u8) -> Self {
        Self {
            equal: match_,
            mismatch,
            gap_open,
            gap_extend,
        }
    }
}

impl AlignmentCostModel for Affine {
    type Kernel = AffineKernel;

    #[inline(always)]
    fn equal(&self) -> u8 {
        self.equal
    }

    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.mismatch
    }

    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.gap_open
    }

    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.gap_extend
    }

    #[inline(always)]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn min_gap_extend(&self) -> u8 {
        self.gap_extend
    }
}
