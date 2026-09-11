//! Query-offset type used by the A* engine.
//!
//! Query offsets are stored in every visited cell and every queued state, so
//! their width dominates memory traffic. `OffsetType` names the handful of
//! operations the engine needs. `u32` is the only implementation.

use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, Sub};

pub trait OffsetType:
    Copy + Ord + Hash + Debug + Default + Add<Output = Self> + Sub<Output = Self>
{
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn increase_one(&self) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
}

impl OffsetType for u32 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        self + 1
    }

    #[inline(always)]
    fn zero() -> Self {
        0
    }

    #[inline(always)]
    fn one() -> Self {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u32_offset_round_trips() {
        let o = u32::new(42);
        assert_eq!(o.as_usize(), 42);
        assert_eq!(o.increase_one(), u32::new(43));
        assert_eq!(u32::zero(), 0);
        assert_eq!(u32::one(), 1);
        assert_eq!(u32::zero() + u32::one(), u32::one());
        assert_eq!(o - u32::one(), u32::new(41));
    }
}
