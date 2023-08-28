use std::fmt::Debug;
use std::hash::Hash;

use num::{FromPrimitive, Unsigned, One, Bounded};
use num::traits::{SaturatingAdd, SaturatingSub};

pub trait OffsetType: FromPrimitive + Unsigned + PartialEq + Eq +
    PartialOrd + Ord + Default + Copy + Hash + Debug + Bounded + SaturatingSub + SaturatingAdd
{
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn as_isize(&self) -> isize;
    fn increase_one(&self) -> Self;
}

impl OffsetType for u8 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u16 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
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
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u64 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }
    
    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}
