use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Bound, Range, Sub, SubAssign};

use num::{FromPrimitive, Unsigned, One, Bounded, PrimInt};
use serde::{Deserialize, Serialize};
use crate::graphs::NodeIndexType;

type DiagType = i64;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Diag(pub DiagType);

impl Diag {
    pub fn new(diag: DiagType) -> Self {
        Diag(diag)
    }

    pub fn new_from_coords(row: DiagType, col: DiagType) -> Self {
        Diag(col - row)
    }

    pub fn diff(&self, rhs: &Self) -> usize {
        (self.0 - rhs.0).abs() as usize
    }
}

impl Add for Diag {
    type Output = Diag;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self::Output {
        Diag(self.0 + rhs.0)
    }
}

impl Sub for Diag {
    type Output = Diag;

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self::Output {
        Diag(self.0 - rhs.0)
    }
}

impl<T: Into<DiagType>> Add<T> for Diag {
    type Output = Diag;

    #[inline(always)]
    fn add(self, rhs: T) -> Self::Output {
        Diag(self.0 + rhs.into())
    }
}

impl<T: Into<DiagType>> Sub<T> for Diag {
    type Output = Diag;

    #[inline(always)]
    fn sub(self, rhs: T) -> Self::Output {
        Diag(self.0 - rhs.into())
    }
}

impl<T: Into<DiagType>> From<T> for Diag {
    fn from(value: T) -> Self {
        Diag(value.into())
    }
}

impl<T: Into<DiagType>> AddAssign<T> for Diag {
    fn add_assign(&mut self, rhs: T) {
        self.0 += rhs.into();
    }
}

impl<T: Into<DiagType>> SubAssign<T> for Diag {
    fn sub_assign(&mut self, rhs: T) {
        self.0 -= rhs.into();
    }
}

pub struct DiagRange {
    k_min: Bound<Diag>,
    k_max: Bound<Diag>
}

/// A struct that represents a range of diagonals for easy iteration
///
/// # Example
///
/// ```
/// use poasta::aligner::offsets::{Diag, DiagRange};
///
/// let range: Vec<_> = DiagRange::closed(Diag(3), Diag(7)).into_iter().collect();
/// assert_eq!(range, vec![Diag(3), Diag(4), Diag(5), Diag(6), Diag(7)]);
///
/// let range2: Vec<_> = DiagRange::half_open(Diag(3), Diag(7)).into_iter().collect();
/// assert_eq!(range2, vec![Diag(3), Diag(4), Diag(5), Diag(6)]);
/// ```
impl DiagRange {
    pub fn closed(k_min: Diag, k_max: Diag) -> Self {
        Self { k_min: Bound::Included(k_min), k_max: Bound::Included(k_max) }
    }

    pub fn half_open(k_min: Diag, k_max: Diag) -> Self {
        Self { k_min: Bound::Included(k_min), k_max: Bound::Excluded(k_max) }
    }
}

impl IntoIterator for DiagRange {
    type Item = Diag;
    type IntoIter = DiagRangeIterator;

    fn into_iter(self) -> Self::IntoIter {
        let start = match self.k_min {
            Bound::Included(d) => d.0,
            Bound::Excluded(d) => d.0 + 1,
            Bound::Unbounded => panic!("Unbounded range not supported!")
        };
        let end = match self.k_max {
            Bound::Included(d) => d.0 + 1,
            Bound::Excluded(d) => d.0,
            Bound::Unbounded => panic!("Unbounded range not supported!")
        };

        Self::IntoIter { iterator: Range { start, end }}
    }
}

pub struct DiagRangeIterator {
    iterator: Range<i64>
}

impl Iterator for DiagRangeIterator {
    type Item = Diag;

    fn next(&mut self) -> Option<Self::Item> {
        self.iterator.next().map(|v| Diag(v.into()))
    }
}

pub trait OffsetType: PrimInt + FromPrimitive + Unsigned + PartialEq + Eq +
    PartialOrd + Ord + Default + Copy + Clone + Hash + Debug + Bounded + Sync + Send +
    AddAssign + SubAssign + Serialize
{
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
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
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for usize {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}


/// A trait that enables easy conversion from (row, offset) tuple to a diagonal
///
/// # Example
/// ```
/// use poasta::aligner::offsets::{Diag, MatrixPoint};
///
/// let point1 = (0 as u32, 0 as u32);
/// let point2 = (1 as u32, 1 as u32);
/// let point3 = (0 as u32, 1 as u32);
///
/// assert_eq!(point1.diag(), point2.diag());
/// assert_eq!(point3.diag(), Diag(1))
/// ```
pub trait MatrixPoint<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    fn row(&self) -> N;
    fn offset(&self) -> O;

    #[inline(always)]
    fn diag(&self) -> Diag {
        Diag::new_from_coords(self.row().as_usize() as i64, self.offset().as_usize() as i64)
    }
}

impl<N, O> MatrixPoint<N, O> for (N, O)
where
    N: NodeIndexType,
    O: OffsetType,
{
    #[inline(always)]
    fn row(&self) -> N {
        self.0
    }

    #[inline(always)]
    fn offset(&self) -> O {
        self.1
    }
}

/// A trait that allows easy conversion from a point on a diagonal to a point in the matrix
///
/// # Example
/// ```
/// use poasta::aligner::offsets::{Diag, DiagonalPoint};
///
/// let point = (Diag(1), 4 as u32);
/// assert_eq!(point.row::<u32>(), 3);
///
/// let point2 = (Diag(-3), 2 as u32);
/// assert_eq!(point2.row::<u32>(), 5);
/// ```
pub trait DiagonalPoint<O>
where
    O: OffsetType,
{
    #[inline(always)]
    fn row<N: NodeIndexType>(&self) -> N {
        N::new((self.offset().as_usize() as i64 - self.diag().0) as usize)
    }

    fn offset(&self) -> O;
    fn diag(&self) -> Diag;
}

impl<O> DiagonalPoint<O> for (Diag, O)
where
    O: OffsetType
{
    #[inline(always)]
    fn offset(&self) -> O {
        self.1
    }

    #[inline(always)]
    fn diag(&self) -> Diag {
        self.0
    }
}

