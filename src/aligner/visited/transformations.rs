use std::marker::PhantomData;

use crate::aligner::offsets::{Diag, OffsetType};
use crate::aligner::visited::AlignState;
use super::{Backtrace, VisitedInterval};

pub mod prelude {
    #[allow(unused_imports)]
    pub use super::{CloneWithBacktraceIterator, IntervalMoveRightIterator, IntervalExtendLengthIterator};
}

pub struct ClonedWithBacktrace<I, O> {
    iter: I,
    backtrace: Backtrace,
    keep_extended: bool,
    dummy: PhantomData<O>,
}

impl<'a, I, V, O> Iterator for ClonedWithBacktrace<I, O>
where
    I: Iterator<Item=&'a V>,
    V: VisitedInterval<O> + 'a,
    O: OffsetType + 'a,
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|v| {
            let mut new = v.clone();
            let mut new_bt = self.backtrace.clone();
            if self.keep_extended && new.backtrace().prev_state == AlignState::Extended {
                new_bt.prev_state = AlignState::Extended;
            }

            *new.backtrace_mut() = new_bt;

            new
        })
    }
}

pub trait CloneWithBacktraceIterator<'a, V, O> : Iterator<Item=&'a V> + Sized
where
    V: VisitedInterval<O> + 'a,
    O: OffsetType + 'a,
{
    fn cloned_with_backtrace(self, backtrace: Backtrace, keep_extended: bool) -> ClonedWithBacktrace<Self, O> {
        ClonedWithBacktrace { iter: self, backtrace, keep_extended, dummy: PhantomData }
    }
}

impl<'a, T, V, O> CloneWithBacktraceIterator<'a, V, O> for T
where
    T: Iterator<Item=&'a V>,
    V: VisitedInterval<O> + 'a,
    O: OffsetType + 'a,
{ }


pub struct TransformMismatch<I, V, O>
where
    I: Iterator<Item=V>,
    V: VisitedInterval<O>,
    O: OffsetType,
{
    input_iter: I,
    new_diag: Diag,
    source_diag: Diag,
    max_offset: O,
}

impl<I, V, O> Iterator for TransformMismatch<I, V, O>
where
    I: Iterator<Item=V>,
    V: VisitedInterval<O>,
    O: OffsetType,
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.input_iter.next().map(|mut ival| {
            // Increase length by one, representing the new mismatched (node, query)
            // Only increase if we don't increase the maximum (i.e., going beyond the query length)
            *ival.end_mut() = std::cmp::min(self.max_offset, ival.end() + O::one());

            // If changing diagonals, reduce length to 1, because otherwise the interval would mark
            // cells as visited incorrectly
            *ival.start_mut() = if self.new_diag == self.source_diag { ival.start() } else { ival.end() - O::one() };

            ival
        })
    }
}

pub trait TransformMismatchIterator<V, O> : Iterator<Item=V> + Sized
where
    V: VisitedInterval<O>,
    O: OffsetType,
{
    fn transform_mismatch(self, new_diag: Diag, source_diag: Diag, max_offset: O) -> TransformMismatch<Self, V, O> {
        TransformMismatch {
            input_iter: self,
            new_diag,
            source_diag,
            max_offset
        }
    }
}

impl<I, V, O> TransformMismatchIterator<V, O> for I
where
    I: Iterator<Item=V>,
    V: VisitedInterval<O>,
    O: OffsetType,
{ }

pub struct TransformInsertion<I, V, O>
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{
    input_iter: I,
    dummy: PhantomData<O>
}

impl<I, V, O> Iterator for TransformInsertion<I, V, O>
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.input_iter.next().map(|mut ival| {
            // For insertions, we simply move the interval one position to the right
            *ival.start_mut() = ival.start().increase_one();
            *ival.end_mut() = ival.end().increase_one();

            ival
        })
    }
}

pub trait TransformInsertionIterator<V, O> : Iterator<Item=V> + Sized
where
    V: VisitedInterval<O>,
    O: OffsetType,
{
    fn transform_insertion(self) -> TransformInsertion<Self, V, O> {
        TransformInsertion {
            input_iter: self,
            dummy: PhantomData
        }
    }
}

impl<I, V, O> TransformInsertionIterator<V, O> for I
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{ }


pub struct TransformDeletion<I, V, O>
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{
    input_iter: I,
    new_diag: Diag,
    source_diag: Diag,
    dummy: PhantomData<O>
}

impl<I, V, O> Iterator for TransformDeletion<I, V, O>
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.input_iter.next().map(|mut ival| {
            // Increase the length if the source interval was clipped.
            if ival.clipped() > O::zero() && ival.start() > O::zero() {
                *ival.start_mut() -= O::one();
                *ival.clipped_mut() -= O::one();
            }

            // For a deletion, we simply move the interval a diagonal down.
            if self.new_diag != self.source_diag - 1 {
                // Clip the interval in case of diagonal switch, other than the one directly below
                let curr_clipped_len = ival.clipped();
                let new_clipped_len = ival.len() - O::one();
                *ival.clipped_mut() += std::cmp::max(curr_clipped_len, new_clipped_len);
                *ival.start_mut() = ival.end() - O::one();
            }

            ival
        })
    }
}

pub trait TransformDeletionIterator<V, O> : Iterator<Item=V> + Sized
where
    V: VisitedInterval<O>,
    O: OffsetType,
{
    fn transform_deletion(self, new_diag: Diag, source_diag: Diag) -> TransformDeletion<Self, V, O> {
        TransformDeletion {
            input_iter: self,
            new_diag,
            source_diag,
            dummy: PhantomData
        }
    }
}

impl<I, V, O> TransformDeletionIterator<V, O> for I
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{ }

/// Iterator adaptor that moves every interval in the input iterator a number of positions
/// to the right.
///
/// For more documentation, see [`IntervalMoveRightIterator::move_right_by'].
pub struct IntervalMoveRight<I, O> {
    iter: I,
    offset: O
}

impl<I, V, O> Iterator for IntervalMoveRight<I, O>
    where
        I: Iterator<Item=V>,
        V: VisitedInterval<O>,
        O: OffsetType,
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|mut v| {
            *v.start_mut() += self.offset;
            *v.end_mut() += self.offset;

            v
        })
    }
}
pub trait IntervalMoveRightIterator<V, O> : Iterator<Item=V> + Sized
where
    V: VisitedInterval<O>,
    O: OffsetType
{
    /// Move every interval in an input iterator to the right
    ///
    /// # Example
    /// ```
    /// use poasta::aligner::offsets::OffsetType;
    /// use poasta::aligner::visited::{VisitedInterval, VisitedIntervalData, VisitedIntervalType};
    /// use poasta::aligner::visited::transformations::prelude::*;
    ///
    /// type Ival = VisitedIntervalType<u32>;
    ///
    /// let intervals = vec![
    ///     Ival::new(2..4, VisitedIntervalData::initial()),
    ///     Ival::new(3..7, VisitedIntervalData::initial()),
    /// ];
    ///
    /// let moved: Vec<_> = intervals.into_iter().move_right_by(3).collect();
    /// assert_eq!(moved, vec![
    ///     Ival::new(5..7, VisitedIntervalData::initial()),
    ///     Ival::new(6..10, VisitedIntervalData::initial()),
    /// ])
    /// ```
    fn move_right_by(self, offset: O) -> IntervalMoveRight<Self, O> {
        IntervalMoveRight { iter: self, offset }
    }
}

impl<T, V, O> IntervalMoveRightIterator<V, O> for T
where
    T: Iterator<Item=V>,
    V: VisitedInterval<O>,
    O: OffsetType,
{ }


/// Iterator adapter that increases the length of each input interval with a set amount.
///
/// For more information, see [`IntervalExtendLengthIterator::extend_len_by'].
pub struct IntervalExtendLength<I, O> {
    iter: I,
    length: O,
    maximum: Option<O>,
}

impl<I, V, O> Iterator for IntervalExtendLength<I, O>
where
    I: Iterator<Item=V>,
    V: VisitedInterval<O>,
    O: OffsetType
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|mut v| {
            if let Some(max) = self.maximum {
                let end = v.end();
                *v.end_mut() = std::cmp::min(end + self.length, max);
            } else {
                *v.end_mut() += self.length;
            }

            v
        })
    }
}

pub trait IntervalExtendLengthIterator<V, O> : Iterator<Item=V> + Sized
where
    V: VisitedInterval<O>,
    O: OffsetType,
{
    /// Extend the length of all intervals in an input iterator.
    ///
    /// # Example
    /// ```
    /// use poasta::aligner::offsets::OffsetType;
    /// use poasta::aligner::visited::{VisitedInterval, VisitedIntervalData, VisitedIntervalType};
    /// use poasta::aligner::visited::transformations::prelude::*;
    ///
    /// type Ival = VisitedIntervalType<u32>;
    ///
    /// let intervals = vec![
    ///     Ival::new(2..4, VisitedIntervalData::initial()),
    ///     Ival::new(3..7, VisitedIntervalData::initial()),
    /// ];
    ///
    /// let increased_len: Vec<_> = intervals.into_iter().extend_len_by(3).collect();
    /// assert_eq!(increased_len, vec![
    ///     Ival::new(2..7, VisitedIntervalData::initial()),
    ///     Ival::new(3..10, VisitedIntervalData::initial()),
    /// ])
    /// ```
    fn extend_len_by(self, length: O) -> IntervalExtendLength<Self, O> {
        IntervalExtendLength { iter: self, length, maximum: None }
    }

    fn extend_len_by_with_max(self, length: O, maximum: O) -> IntervalExtendLength<Self, O> {
        IntervalExtendLength { iter: self, length, maximum: Some(maximum) }
    }
}

impl<T, V, O> IntervalExtendLengthIterator<V, O> for T
where
    T: Iterator<Item=V>,
    V: VisitedInterval<O>,
    O: OffsetType,
{ }
