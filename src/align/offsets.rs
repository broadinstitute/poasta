//! Query-offset type used by the A* engine.
//!
//! Query offsets are stored in every visited cell and every queued state, so
//! their width dominates memory traffic. `OffsetType` names the handful of
//! operations the engine needs. It is implemented for `u16`, `u32` and
//! `u64`. The CLI uses `u32` (`astar::Offset`).

use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, Sub};

pub trait OffsetType:
    Copy + Ord + Hash + Debug + Default + Add<Output = Self> + Sub<Output = Self>
{
    /// Largest offset the type holds. The engine rejects longer queries.
    const MAX: usize;

    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn increase_one(&self) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
}

macro_rules! impl_offset_type {
    ($($t:ty),*) => {$(
        impl OffsetType for $t {
            // `u64::MAX as usize` truncates on 32-bit targets, where no
            // query is that long anyway.
            const MAX: usize = <$t>::MAX as usize;

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
    )*};
}

impl_offset_type!(u16, u32, u64);

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trips<O: OffsetType>() {
        let o = O::new(42);
        assert_eq!(o.as_usize(), 42);
        assert_eq!(o.increase_one(), O::new(43));
        assert_eq!(O::zero().as_usize(), 0);
        assert_eq!(O::one().as_usize(), 1);
        assert_eq!(O::zero() + O::one(), O::one());
        assert_eq!(o - O::one(), O::new(41));
        assert_eq!(O::new(O::MAX).as_usize(), O::MAX);
    }

    #[test]
    fn offsets_round_trip() {
        round_trips::<u16>();
        round_trips::<u32>();
        round_trips::<u64>();
    }
}
