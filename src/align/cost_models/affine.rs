use super::AlignmentCostModel;
use crate::align::kernels::affine::AffineKernel;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Affine {
    mismatch: u8,
    gap_open: u8,
    gap_extend: u8,
}

impl Affine {
    pub fn new(mismatch: u8, gap_open: u8, gap_extend: u8) -> Self {
        Self {
            mismatch,
            gap_open,
            gap_extend,
        }
    }
}

impl AlignmentCostModel for Affine {
    type Kernel = AffineKernel;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::AlignState;

    fn costs() -> Affine {
        Affine::new(4, 2, 1) // mismatch=4, gap_open=2, gap_extend=1
    }

    #[test]
    fn gap_cost_zero_length_is_zero() {
        assert_eq!(costs().gap_cost(AlignState::Match, 0), 0);
    }

    #[test]
    fn gap_cost_from_match_charges_open_and_extend() {
        assert_eq!(costs().gap_cost(AlignState::Match, 3), 2 + 3);
    }

    #[test]
    fn gap_cost_already_in_gap_charges_extend_only() {
        assert_eq!(costs().gap_cost(AlignState::Insertion, 3), 3);
        assert_eq!(costs().gap_cost(AlignState::Deletion, 3), 3);
    }
}
