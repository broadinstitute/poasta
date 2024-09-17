use std::collections::VecDeque;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{BitAnd, Not, Shr, Add, Sub};

use num::{FromPrimitive, Unsigned, One, Bounded};
use num::traits::{SaturatingAdd, SaturatingSub};

/// Index types are integer types to represent the diagonal or query position 
pub trait IndexType: FromPrimitive + Unsigned + PartialEq + Eq
    + PartialOrd + Ord + Default + Copy + Hash + Debug + Bounded 
    + Add + Sub + SaturatingSub + SaturatingAdd
    + Not<Output=Self> + BitAnd<Output=Self> + Shr<Output=Self> + 'static
{
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn as_isize(&self) -> isize;
    fn increase_one(&self) -> Self;
}

impl IndexType for u8 {
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

impl IndexType for u16 {
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

impl IndexType for u32 {
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

impl IndexType for u64 {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Diag<D>(D);

impl<D> Diag<D>
    where D: IndexType
{
    #[inline(always)]
    pub fn as_usize(&self) -> usize {
        self.0.as_usize()
    }
}

impl<D> Add for Diag<D>
    where D: IndexType
{
    type Output = Self;
    
    fn add(self, rhs: Self) -> Self::Output {
        Diag(self.0 + rhs.0)
    }
}

impl<D> Sub for Diag<D>
    where D: IndexType
{
    type Output = Self;
    
    fn sub(self, rhs: Self) -> Self::Output {
        Diag(self.0 - rhs.0)
    }
}

impl<D> Add<D> for Diag<D>
    where D: IndexType
{
    type Output = Self;
    
    fn add(self, rhs: D) -> Self::Output {
        Diag(self.0 + rhs)
    }
}

impl<D> Sub<D> for Diag<D>
    where D: IndexType
{
    type Output = Self;
    
    fn sub(self, rhs: D) -> Self::Output {
        Diag(self.0 - rhs)
    }
}


/// Holds the furthest reaching points for a node in the POA graph
#[derive(Debug, Clone)]
pub struct NodeFrPoints<D> {
    fr_points: VecDeque<D>,
    kmin: Diag<D>,
    kmax: Diag<D>,
}

impl<D> NodeFrPoints<D>
    where D: IndexType
{
    pub fn new() -> Self {
        Self {
            fr_points: VecDeque::new(),
            kmin: Diag(D::default()),
            kmax: Diag(D::default()),
        }
    }
    
    pub fn len(&self) -> usize {
        self.fr_points.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.fr_points.is_empty()
    }
    
    pub fn get(&self, diag: Diag<D>) -> Option<D> {
        let ix = diag - self.kmin;
        
        self.fr_points.get(ix.as_usize()).copied()
    }
    
    pub fn is_further(&self, diag: Diag<D>, offset: D) -> bool {
        if self.is_empty() {
            return true;
        }
        
        if self.kmin > diag || self.kmax < diag {
            return true;
        }
        
        self.get(diag).unwrap() < offset
    }

    pub fn kmin(&self) -> Diag<D> {
        self.kmin
    }

    pub fn kmax(&self) -> Diag<D> {
        self.kmax
    }

    pub fn fr_points(&self) -> &VecDeque<D> {
        &self.fr_points
    }
}

impl<D> Default for NodeFrPoints<D>
    where D: Default
{
    fn default() -> Self {
        Self {
            fr_points: VecDeque::default(),
            kmin: Diag(D::default()),
            kmax: Diag(D::default()),
        }
    }
}
