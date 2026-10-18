use super::AlignmentCostModel;

/// Linear gap model: gap cost = gap_extend × length (no open penalty).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Linear {
    mismatch: u8,
    gap_extend: u8,
}

impl Linear {
    pub fn new(mismatch: u8, gap_extend: u8) -> Self {
        Self {
            mismatch,
            gap_extend,
        }
    }
}

impl AlignmentCostModel for Linear {
    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.mismatch
    }

    #[inline(always)]
    fn gap_open(&self) -> u8 {
        0
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
}
