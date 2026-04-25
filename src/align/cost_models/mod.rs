pub mod affine;
pub mod linear;
pub mod two_piece;

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
}
