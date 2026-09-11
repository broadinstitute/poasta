use super::AlignmentCostModel;
use crate::align::AlignState;
use crate::align::kernels::two_piece::TwoPieceAffineKernel;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TwoPieceAffine {
    mismatch: u8,
    gap_open: u8,
    gap_extend: u8,
    gap_open2: u8,
    gap_extend2: u8,
}

impl TwoPieceAffine {
    pub fn new(mismatch: u8, gap_open: u8, gap_extend: u8, gap_open2: u8, gap_extend2: u8) -> Self {
        Self {
            mismatch,
            gap_open,
            gap_extend,
            gap_open2,
            gap_extend2,
        }
    }
}

impl AlignmentCostModel for TwoPieceAffine {
    type Kernel = TwoPieceAffineKernel;
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
        self.gap_open2
    }
    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        self.gap_extend2
    }
    #[inline(always)]
    fn min_gap_extend(&self) -> u8 {
        self.gap_extend.min(self.gap_extend2)
    }

    /// Override: pick the cheaper of the two pieces when opening from `Match`;
    /// no new open cost when already inside a gap of the matching piece —
    /// same "already open, extend only" shape as the single-piece default,
    /// applied per-piece. (Both tracks open independently from `Match`, per
    /// this crate's `TwoPieceAffineKernel` recurrence — see
    /// `engine::astar::aln_graph::TwoPieceAlignmentGraph`'s doc comment for
    /// why this differs from `stable`'s transition-based model.)
    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize {
        if length == 0 {
            return 0;
        }

        match current_state {
            AlignState::Insertion | AlignState::Deletion => length * self.gap_extend as usize,
            AlignState::Insertion2 | AlignState::Deletion2 => length * self.gap_extend2 as usize,
            AlignState::Match => {
                let cost1 = self.gap_open as usize + (length * self.gap_extend as usize);
                let cost2 = self.gap_open2 as usize + (length * self.gap_extend2 as usize);
                cost1.min(cost2)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn costs() -> TwoPieceAffine {
        TwoPieceAffine::new(4, 10, 2, 8, 1)
    }

    #[test]
    fn gap_cost_zero_length_is_zero() {
        assert_eq!(costs().gap_cost(AlignState::Match, 0), 0);
    }

    #[test]
    fn gap_cost_from_match_takes_cheaper_piece() {
        let c = costs();
        assert_eq!(c.gap_cost(AlignState::Match, 1), 9); // min(10+2*1, 8+1*1) = min(12, 9) = 9
        assert_eq!(c.gap_cost(AlignState::Match, 2), 10); // min(14, 10) = 10
        assert_eq!(c.gap_cost(AlignState::Match, 3), 11); // min(16, 11) = 11
    }

    #[test]
    fn gap_cost_within_first_piece_has_no_reopen_cost() {
        let c = costs();
        assert_eq!(c.gap_cost(AlignState::Insertion, 1), 2); // extend only: 1*2
        assert_eq!(c.gap_cost(AlignState::Deletion, 2), 4); // 2*2
    }

    #[test]
    fn gap_cost_within_second_piece_has_no_reopen_cost() {
        let c = costs();
        assert_eq!(c.gap_cost(AlignState::Insertion2, 1), 1); // 1*1
        assert_eq!(c.gap_cost(AlignState::Deletion2, 2), 2); // 2*1
    }
}
