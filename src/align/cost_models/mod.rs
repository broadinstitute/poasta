use std::fmt::Debug;
use std::sync::Arc;


pub mod affine;

pub trait AlignmentCostModel {
    fn mismatch(&self) -> u8;

    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;

    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;

    fn gap_cost(&self, current_state: AlignState, gap_length: usize) -> usize;
}
