use std::fmt::Debug;

use super::AlignmentCostModel;

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
    fn gap_cost(&self, current_state: AlignState, gap_length: usize) -> usize {
        if gap_length == 0 {
            return 0;
        }

        let gap_open = match current_state {
            AlignState::Insertion | AlignState::Deletion => 0,
            AlignState::Match => self.gap_open,
            _ => panic!("Invalid current state {current_state:?} for gap affine scoring model!"),
        };

        gap_open as usize + (gap_length * self.gap_extend as usize)
    }
}

