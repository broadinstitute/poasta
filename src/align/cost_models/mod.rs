//! Cost models: what a mismatch and a gap cost.
//!
//! [`AlignmentCostModel`] exposes costs for up to two gap pieces; the
//! [`linear`], [`affine`] and [`two_piece`] models implement it. The A*
//! engine attaches its per-model machinery through `astar::AstarCostModel`.

pub mod affine;
pub mod linear;
pub mod two_piece;

use crate::align::AlignState;

pub trait AlignmentCostModel {
    fn mismatch(&self) -> u8;

    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;

    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;

    /// Cost of a gap of `length` given the current alignment state.
    ///
    /// `length == 0` always costs `0`. Already inside a gap (`Insertion` or
    /// `Deletion`) no new gap-open penalty is charged, the gap only extends.
    /// Two-piece models extend within the piece they are in and take the
    /// cheaper piece when opening from `Match`.
    ///
    /// Used by the A* heuristic (minimum remaining gap cost) and by the
    /// bubble-exit pruning check.
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
