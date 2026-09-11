pub mod affine;
pub mod linear;
pub mod two_piece;

use crate::align::AlignState;
use crate::align::kernels::DPKernel;

pub trait AlignmentCostModel {
    type Kernel: DPKernel<Costs = Self>;

    fn mismatch(&self) -> u8;

    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;

    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;

    /// Smallest per-step gap cost across all gap modes the model supports.
    /// Used as the denominator of the Ukkonen termination bound
    /// `indels ≤ score / min_gap_extend`.
    fn min_gap_extend(&self) -> u8;

    /// Cost of a gap of `length` given the current alignment state.
    ///
    /// `length == 0` always costs `0`. If already inside a gap (`Insertion` /
    /// `Deletion`), no new gap-open penalty is charged — the gap is merely
    /// extended. Two-piece models charge the matching piece's open+extend
    /// when already inside that piece's gap state, and take the cheaper of
    /// the two pieces when opening a fresh gap from `Match`.
    ///
    /// Used by the A* heuristic (minimum remaining gap cost) and by the
    /// bubble-exit pruning check, both of which only ever query this in the
    /// `Match` state; the per-state branches below mirror `stable`'s scoring
    /// model for parity with its `AstarHeuristic`/bubble-pruning code.
    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize {
        if length == 0 {
            return 0;
        }

        let gap_open = match current_state {
            AlignState::Insertion | AlignState::Deletion => 0,
            AlignState::Match => self.gap_open() as usize,
            AlignState::Insertion2 | AlignState::Deletion2 => {
                panic!("Invalid state {current_state:?} for a single-piece gap cost model!")
            }
        };

        gap_open + (length * self.gap_extend() as usize)
    }
}
