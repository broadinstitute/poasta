use std::collections::VecDeque;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, AddAssign, BitAnd, Not, Shr, Sub, SubAssign};

use num::traits::{SaturatingAdd, SaturatingSub};
use num::{Bounded, FromPrimitive, One, Signed, Unsigned};

pub trait NumOperations:
    FromPrimitive
    + PartialEq
    + Eq
    + PartialOrd
    + Ord
    + Default
    + Copy
    + Hash
    + Debug
    + Bounded
    + Add
    + Sub
    + SaturatingSub
    + SaturatingAdd
    + Not<Output = Self>
    + BitAnd<Output = Self>
    + Shr<Output = Self>
    + Send
    + Sync
    + 'static
{
}

impl NumOperations for u8 {}
impl NumOperations for u16 {}
impl NumOperations for u32 {}
impl NumOperations for u64 {}

impl NumOperations for i8 {}
impl NumOperations for i16 {}
impl NumOperations for i32 {}
impl NumOperations for i64 {}

/// Position types are integer types to represent a query position
///
/// The most significant bit is used as a `visited` flag. 
pub trait PosType: NumOperations + Unsigned {
    /// Create a new offset from a usize, saturating to the maximum value if the input is too large
    fn new(value: usize) -> Self;
    
    /// Return the value without the `visited` bit
    fn value(&self) -> Self;
    
    /// Get the value as usize
    fn as_usize(&self) -> usize;
    
    /// Get the value as isize
    fn as_isize(&self) -> isize;
    
    /// Increase the offset by one, resets the `visited` bit
    fn increase_one(&self) -> Self;
    
    /// Check if the `visited` bit is set
    fn is_visited(&self) -> bool; 
    
    /// Set the `visited` bit
    fn set_visited(&mut self, visited: bool);
    
    fn max() -> Self;
}

impl PosType for u8 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        (value as Self) & Self::MAX >> 1
    }
    
    #[inline(always)]
    fn value(&self) -> Self {
        *self & (Self::MAX >> 1)
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        let value = *self & (Self::MAX >> 1);
        value as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        let value = *self & (Self::MAX >> 1);
        value as isize
    }

    /// Increase the offset by one, resets the `visited` bit
    #[inline(always)]
    fn increase_one(&self) -> Self {
        let value = *self & (Self::MAX >> 1);
        if value > <Self as PosType>::max() {
            panic!("Overflow in increase_one");
        }
        
        value + Self::one()
    }
    
    #[inline(always)]
    fn is_visited(&self) -> bool {
        *self & (1 << (Self::BITS - 1)) != 0
    }
    
    #[inline(always)]
    fn set_visited(&mut self, visited: bool) {
        if visited {
            *self |= 1 << (Self::BITS - 1);
        } else {
            *self &= !(1 << (Self::BITS - 1));
        }
    }
    
    #[inline(always)]
    fn max() -> Self {
        Self::MAX >> 1
    }
}

impl PosType for u16 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        (value as Self) & Self::MAX >> 1
    }
    
    #[inline(always)]
    fn value(&self) -> Self {
        *self & (Self::MAX >> 1)
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        let value = *self & (Self::MAX >> 1);
        value as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        let value = *self & (Self::MAX >> 1);
        value as isize
    }

    /// Increase the offset by one, resets the `visited` bit
    #[inline(always)]
    fn increase_one(&self) -> Self {
        let value = *self & (Self::MAX >> 1);
        if value > <Self as PosType>::max() {
            panic!("Overflow in increase_one");
        }
        
        value + Self::one()
    }
    
    #[inline(always)]
    fn is_visited(&self) -> bool {
        *self & (1 << (Self::BITS - 1)) != 0
    }
    
    #[inline(always)]
    fn set_visited(&mut self, visited: bool) {
        if visited {
            *self |= 1 << (Self::BITS - 1);
        } else {
            *self &= !(1 << (Self::BITS - 1));
        }
    }
    
    #[inline(always)]
    fn max() -> Self {
        Self::MAX >> 1
    }
}

impl PosType for u32 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        (value as Self) & Self::MAX >> 1
    }
    
    #[inline(always)]
    fn value(&self) -> Self {
        *self & (Self::MAX >> 1)
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        let value = *self & (Self::MAX >> 1);
        value as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        let value = *self & (Self::MAX >> 1);
        value as isize
    }

    /// Increase the offset by one, resets the `visited` bit
    #[inline(always)]
    fn increase_one(&self) -> Self {
        let value = *self & (Self::MAX >> 1);
        if value > <Self as PosType>::max() {
            panic!("Overflow in increase_one");
        }
        
        value + Self::one()
    }
    
    #[inline(always)]
    fn is_visited(&self) -> bool {
        *self & (1 << (Self::BITS - 1)) != 0
    }
    
    #[inline(always)]
    fn set_visited(&mut self, visited: bool) {
        if visited {
            *self |= 1 << (Self::BITS - 1);
        } else {
            *self &= !(1 << (Self::BITS - 1));
        }
    }
    
    #[inline(always)]
    fn max() -> Self {
        Self::MAX >> 1
    }
}

impl PosType for u64 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        (value as Self) & Self::MAX >> 1
    }
    
    #[inline(always)]
    fn value(&self) -> Self {
        *self & (Self::MAX >> 1)
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        let value = *self & (Self::MAX >> 1);
        value as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        let value = *self & (Self::MAX >> 1);
        value as isize
    }

    /// Increase the offset by one, resets the `visited` bit
    #[inline(always)]
    fn increase_one(&self) -> Self {
        let value = *self & (Self::MAX >> 1);
        if value > <Self as PosType>::max() {
            panic!("Overflow in increase_one");
        }
        
        value + Self::one()
    }
    
    #[inline(always)]
    fn is_visited(&self) -> bool {
        *self & (1 << (Self::BITS - 1)) != 0
    }
    
    #[inline(always)]
    fn set_visited(&mut self, visited: bool) {
        if visited {
            *self |= 1 << (Self::BITS - 1);
        } else {
            *self &= !(1 << (Self::BITS - 1));
        }
    }
    
    #[inline(always)]
    fn max() -> Self {
        Self::MAX >> 1
    }
}

pub trait DiagType: NumOperations + Signed {
    fn new(value: isize) -> Self;
    fn as_isize(&self) -> isize;
    fn as_usize(&self) -> usize;

    fn increase_one(&self) -> Self;
    fn decrease_one(&self) -> Self;
}

impl DiagType for i8 {
    #[inline(always)]
    fn new(value: isize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }

    #[inline(always)]
    fn decrease_one(&self) -> Self {
        *self - Self::one()
    }
}

impl DiagType for i16 {
    #[inline(always)]
    fn new(value: isize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }

    #[inline(always)]
    fn decrease_one(&self) -> Self {
        *self - Self::one()
    }
}

impl DiagType for i32 {
    #[inline(always)]
    fn new(value: isize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }

    #[inline(always)]
    fn decrease_one(&self) -> Self {
        *self - Self::one()
    }
}

impl DiagType for i64 {
    #[inline(always)]
    fn new(value: isize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }

    #[inline(always)]
    fn decrease_one(&self) -> Self {
        *self - Self::one()
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Diag<D>(D);

impl<D> Diag<D>
where
    D: DiagType,
{
    #[inline(always)]
    pub fn new(v: isize) -> Self {
        Diag(D::new(v))
    }

    #[inline(always)]
    pub fn as_usize(&self) -> usize {
        self.0.as_usize()
    }

    #[inline(always)]
    pub fn as_isize(&self) -> isize {
        self.0.as_isize()
    }
}

impl<D> Add<usize> for Diag<D>
where
    D: DiagType,
{
    type Output = Self;

    fn add(self, rhs: usize) -> Self::Output {
        Diag(D::new(self.0.as_isize() + rhs as isize))
    }
}

impl<D> Add<isize> for Diag<D>
where
    D: DiagType,
{
    type Output = Self;

    fn add(self, rhs: isize) -> Self::Output {
        Diag(D::new(self.0.as_isize() + rhs))
    }
}

impl<D> Add<Diag<D>> for Diag<D>
where
    D: DiagType,
{
    type Output = Self;

    fn add(self, rhs: Diag<D>) -> Self {
        Diag(self.0 + rhs.0)
    }
}

impl<D> AddAssign<usize> for Diag<D>
where
    D: DiagType,
{
    fn add_assign(&mut self, rhs: usize) {
        self.0 = self.0 + D::new(rhs as isize);
    }
}

impl<D> AddAssign<isize> for Diag<D>
where
    D: DiagType,
{
    fn add_assign(&mut self, rhs: isize) {
        self.0 = self.0 + D::new(rhs);
    }
}

impl<D> Sub<usize> for Diag<D>
where
    D: DiagType,
{
    type Output = Self;

    fn sub(self, rhs: usize) -> Self::Output {
        Diag(D::new(self.0.as_isize() - rhs as isize))
    }
}

impl<D> Sub<isize> for Diag<D>
where
    D: DiagType,
{
    type Output = Self;

    fn sub(self, rhs: isize) -> Self::Output {
        Diag(D::new(self.0.as_isize() - rhs))
    }
}

impl<D> Sub<Diag<D>> for Diag<D>
where
    D: DiagType,
{
    type Output = Self;

    fn sub(self, rhs: Diag<D>) -> Self::Output {
        Diag(self.0 - rhs.0)
    }
}

impl<D> SubAssign<usize> for Diag<D>
where
    D: DiagType,
{
    fn sub_assign(&mut self, rhs: usize) {
        self.0 = self.0 - D::new(rhs as isize);
    }
}

impl<D> SubAssign<isize> for Diag<D>
where
    D: DiagType,
{
    fn sub_assign(&mut self, rhs: isize) {
        self.0 = self.0 - D::new(rhs);
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Score(u32);

impl Score {
    #[inline(always)]
    pub fn new(value: u32) -> Self {
        Score(value)
    }

    #[inline(always)]
    pub fn score(&self) -> u32 {
        self.0
    }

    #[inline(always)]
    pub fn as_usize(&self) -> usize {
        self.0 as usize
    }
}

impl Add<usize> for Score {
    type Output = Self;

    fn add(self, rhs: usize) -> Self {
        Score(self.0 + rhs as u32)
    }
}

impl Add<u8> for Score {
    type Output = Self;

    fn add(self, rhs: u8) -> Self {
        Score(self.0 + rhs as u32)
    }
}

impl Sub<usize> for Score {
    type Output = Self;

    fn sub(self, rhs: usize) -> Self {
        Score(self.0 - rhs as u32)
    }
}

impl Sub<u8> for Score {
    type Output = Self;

    fn sub(self, rhs: u8) -> Self {
        Score(self.0 - rhs as u32)
    }
}

impl Add<Score> for Score {
    type Output = Self;

    fn add(self, rhs: Score) -> Self {
        Score(self.0 + rhs.0)
    }
}

impl Sub<Score> for Score {
    type Output = Self;

    fn sub(self, rhs: Score) -> Self {
        Score(self.0 - rhs.0)
    }
}

impl From<u8> for Score {
    fn from(value: u8) -> Self {
        Score(value as u32)
    }
}

#[derive(Debug, Clone, Default)]
pub struct Diagonals<D, O> {
    /// The furthest reached query position for each diagonal.
    diagonals: VecDeque<O>,
    
    /// The smallest reached diagonal. `diagonals[0]` represents this diagonal.
    kmin: Diag<D>,
}

impl<D, O> Diagonals<D, O>
where
    D: DiagType,
    O: PosType,
{
    
    pub fn len(&self) -> usize {
        self.diagonals.len()
    }

    pub fn is_empty(&self) -> bool {
        self.diagonals.is_empty()
    }

    fn ensure_space(&mut self, diag: Diag<D>) {
        if self.is_empty() {
            self.diagonals.resize(8, O::default());
            self.kmin = diag;
            return;
        }

        let kmax = self.kmin + self.len() - 1usize;
        if diag < self.kmin {
            let extra = (self.kmin - diag).as_usize();
            self.diagonals.reserve(extra);
            for _ in 0..extra {
                self.diagonals.push_front(O::default());
            }
            self.kmin = diag;
        } else if diag > kmax {
            let extra = (diag - kmax).as_usize();
            self.diagonals.reserve(extra);
            for _ in 0..extra {
                self.diagonals.push_back(O::default());
            }
        }
    }

    pub fn get_furthest(&self, diag: Diag<D>) -> Option<O> {
        let ix = (diag - self.kmin).as_usize();

        self.diagonals.get(ix).copied()
    }
    
    pub fn set_furthest(&mut self, diag: Diag<D>, offset: O) {
        self.ensure_space(diag);

        let ix = (diag - self.kmin).as_usize();
        self.diagonals[ix] = offset;
    }

    pub fn is_further(&self, diag: Diag<D>, offset: O) -> bool {
        if self.is_empty() {
            return true;
        }

        let ix = (diag - self.kmin).as_usize();
        self.diagonals.get(ix).map(|&o| o < offset).unwrap_or(true)
    }

    pub fn update_if_further(&mut self, diag: Diag<D>, offset: O) -> bool {
        self.ensure_space(diag);

        let ix = (diag - self.kmin).as_usize();

        if self.diagonals[ix] < offset {
            self.diagonals[ix] = offset;
            true
        } else {
            false
        }
    }
    
    pub fn is_visited(&self, diag: Diag<D>) -> bool {
        if diag < self.kmin || diag >= self.kmin + self.len() {
            return false;
        }
        
        let ix = (diag - self.kmin).as_usize();
        self.diagonals[ix].is_visited()
    }
    
    pub fn set_visited(&mut self, diag: Diag<D>, visited: bool) {
        self.ensure_space(diag);
        
        let ix = (diag - self.kmin).as_usize();
        self.diagonals[ix].set_visited(visited);
    }
}

/// Holds reached query positions for node diagonals.
#[derive(Debug, Clone)]
pub struct NodeFrPoints<D, O> {
    /// Per score, the furthest reached query position for each diagonal.
    fr_points: VecDeque<Diagonals<D, O>>,

    /// The lowest score reached in this node. fr_points[0] refers to furthest reached points at this score.
    score_min: Score,
}

impl<D, O> NodeFrPoints<D, O>
where
    D: DiagType,
    O: PosType,
{
    pub fn is_empty(&self) -> bool {
        self.fr_points.is_empty()
    }

    pub fn get_furthest(&self, score: Score, diag: Diag<D>) -> Option<O> {
        if score < self.score_min {
            return None;
        }
        
        let ix = (score - self.score_min).as_usize();

        self.fr_points.get(ix).and_then(|v| v.get_furthest(diag))
    }
    
    pub fn set_furthest(&mut self, score: Score, diag: Diag<D>, offset: O) {
        self.ensure_space(score);
        
        let ix = (score - self.score_min).as_usize();
        self.fr_points[ix].set_furthest(diag, offset)
    }

    fn ensure_space(&mut self, score: Score) {
        if self.is_empty() {
            self.fr_points.resize(8, Diagonals::default());
            self.score_min = score;
            return;
        }

        let score_max = self.score_min + self.fr_points.len() - 1usize;
        if score < self.score_min {
            let extra = (self.score_min - score).as_usize();
            self.fr_points.reserve(extra);
            for _ in 0..extra {
                self.fr_points.push_front(Diagonals::default());
            }
            self.score_min = score;
        } else if score > score_max {
            let extra = (score - score_max).as_usize();
            self.fr_points.reserve(extra);
            for _ in 0..extra {
                self.fr_points.push_back(Diagonals::default());
            }
        }
    }

    pub fn update_if_further(&mut self, score: Score, diag: Diag<D>, offset: O) -> bool {
        self.ensure_space(score);
        let ix = (score - self.score_min).as_usize();

        self.fr_points[ix].update_if_further(diag, offset)
    }

    /// Is a given query position further than any currently reached points for a specific diagonal?
    pub fn is_further(&self, score: Score, diag: Diag<D>, offset: O) -> bool {
        if self.is_empty() {
            return true;
        }
        
        if score < self.score_min {
            return true;
        }

        let ix = (score - self.score_min).as_usize();
        self.fr_points
            .get(ix)
            .map(|v| v.is_further(diag, offset))
            .unwrap_or(true)
    }
    
    pub fn is_visited(&self, score: Score, diag: Diag<D>) -> bool {
        if self.is_empty() {
            return false;
        }
        
        if score < self.score_min {
            return false;
        }
        
        let ix = (score - self.score_min).as_usize();
        self.fr_points
            .get(ix)
            .map(|v| v.is_visited(diag))
            .unwrap_or(false)
    }
    
    pub fn set_visited(&mut self, score: Score, diag: Diag<D>, visited: bool) {
        self.ensure_space(score);
        
        let ix = (score - self.score_min).as_usize();
        self.fr_points[ix].set_visited(diag, visited);
    }
}

impl<D, O> Default for NodeFrPoints<D, O>
where
    D: Default,
{
    fn default() -> Self {
        Self {
            fr_points: VecDeque::default(),
            score_min: Score::default(),
        }
    }
}

#[inline(always)]
pub fn to_node_pos<D>(diag: Diag<D>, query_offset: usize) -> usize
where
    D: DiagType,
{
    query_offset
        .checked_add_signed(-diag.as_isize())
        .unwrap()
}
