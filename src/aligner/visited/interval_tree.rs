//! Interval tree, a data structure for efficiently storing and searching intervals.
//!
//! This implementation is based on the implementation in rust-bio, which in turn is based on
//! the sorted array version as described/given in
//! https://github.com/lh3/cgranges / https://github.com/lh3/cgranges/blob/master/cpp/IITree.h
//!
//! # Example
//! ```
//! use poasta::aligner::visited::interval_tree::{Interval, IntervalTree};
//! use poasta::aligner::offsets::OffsetType;
//! use std::iter::FromIterator;
//!
//! let mut tree = IntervalTree::<usize, ()>::from_iter(vec![(12..15, ()), (25..32, ()), (1..5, ())]);
//! assert!(tree.contains(12));
//! assert!(tree.contains(27));
//! assert!(!tree.contains(15));
//! assert!(!tree.contains(20));
//!
//! let mut tree2 = IntervalTree::<usize, ()>::from_sorted(vec![(1..5, ()), (12..15, ()), (25..32, ())]);
//! assert!(tree2.contains(12));
//! assert!(tree2.contains(27));
//! assert!(!tree2.contains(15));
//! assert!(!tree2.contains(20));
//! ```

use std::hash::{Hash, Hasher};
use std::cmp::{min, Ordering};
use std::iter::{FromIterator, Map};
use std::ops::{Deref, Range};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use crate::aligner::offsets::OffsetType;


/// A type representing an half open interval, ranging from [start, end).
///
/// Wraps a Rust `std::ops::Range` object, and adds some convenience methods.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Interval<O, D>(Range<O>, D);

impl<O, D> Interval<O, D>
where
    O: OffsetType
{
    pub fn new(range: Range<O>, data: D) -> Self {
        assert!(range.start < range.end);

        Self(range, data)
    }

    pub fn new_size_one(pos: O, data: D) -> Self {
        Self::new(pos..(pos + O::one()), data)
    }

    #[inline(always)]
    pub fn start_mut(&mut self) -> &mut O {
        &mut self.0.start
    }

    #[inline(always)]
    pub fn end_mut(&mut self) -> &mut O {
        &mut self.0.end
    }

    #[inline(always)]
    pub fn data(&self) -> &D {
        &self.1
    }

    #[inline(always)]
    pub fn data_mut(&mut self) -> &mut D {
        &mut self.1
    }

    #[inline]
    pub fn overlaps(&self, other: &Self) -> bool {
        (self.start <= other.start && self.end > other.start)
            || (self.start > other.start && other.end > self.start)
    }

    #[inline]
    pub fn overlaps_or_adjacent(&self, other: &Self) -> bool {
        (self.start <= other.start && self.end >= other.start)
            || (self.start > other.start && other.end >= self.start)
    }

    #[inline]
    pub fn contains(&self, pos: O) -> bool {
        pos >= self.start && pos < self.end
    }
}

impl<O, D> Default for Interval<O, D>
    where
        O: OffsetType,
        D: Default,
{
    fn default() -> Self {
        Self(O::zero()..O::one(), D::default())
    }
}

impl<O, D> PartialOrd for Interval<O, D>
    where
        O: OffsetType,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<O, D> Ord for Interval<O, D>
    where
        O: OffsetType
{
    fn cmp(&self, other: &Self) -> Ordering {
        if self.start < other.start {
            Ordering::Less
        } else if self.start == other.start {
            // Flip comparison, make sure larger intervals come up first
            if self.end > other.end {
                Ordering::Less
            } else if self.end == other.end {
                Ordering::Equal
            } else {
                Ordering::Greater
            }
        } else {
            Ordering::Greater
        }
    }
}

impl<O, D> PartialEq for Interval<O, D>
    where
        O: OffsetType
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<O, D> Hash for Interval<O, D>
where
    O: OffsetType
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state)
    }
}

impl<O: OffsetType, D> Eq for Interval<O, D> { }

/// Convert a `Range` into an `Interval`. This conversion will panic if the `Range` has end < start
impl<O: OffsetType, D> From<(Range<O>, D)> for Interval<O, D> {
    fn from((r, d): (Range<O>, D)) -> Self {
        Interval::new(r, d)
    }
}

/// Use the `Deref` operator to get a reference to `Range` wrapped by the `Interval` newtype.
impl<O: OffsetType, D> Deref for Interval<O, D> {
    type Target = Range<O>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}


#[derive(Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
pub struct InternalEntry<O: OffsetType, D> {
    interval: Interval<O, D>,
    max: O,
}

impl<O: OffsetType, D, const S: usize> Default for IntervalTree<O, D, S> {
    fn default() -> Self {
        IntervalTree {
            entries: SmallVec::default(),
            max_level: 0,
            indexed: true,
        }
    }
}

#[derive(Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
pub struct IntervalTree<O: OffsetType, D, const S: usize = 4> {
    entries: SmallVec<[InternalEntry<O, D>; S]>,
    max_level: usize,
    indexed: bool,
}

impl<O, D, V> FromIterator<V> for IntervalTree<O, D>
    where
        V: Into<Interval<O, D>>,
        O: OffsetType,
{
    fn from_iter<T: IntoIterator<Item=V>>(iter: T) -> Self {
        let mut tree = Self::new();
        iter.into_iter()
            .for_each(|ival| tree.insert(ival.into()));
        tree.sort_and_index();
        tree
    }
}

impl<O: OffsetType, D, const S: usize> IntervalTree<O, D, S> {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn from_sorted<T, I>(container: T) -> Self
    where
        T: IntoIterator<Item=I>,
        I: Into<Interval<O, D>>,
        O: OffsetType
    {
        let entries: SmallVec<[InternalEntry<O, D>; S]> = container.into_iter()
            .map(|ival| {
                let interval = ival.into();
                let max = interval.end;
                InternalEntry {
                    interval,
                    max
                }
            })
            .collect();

        let mut last_start = O::min_value();
        for entry in entries.iter() {
            assert!(entry.interval.start >= last_start);
            last_start = entry.interval.start;
        }

        let mut tree = Self { entries, max_level: 0, indexed: true };
        tree.index_core();

        tree
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item=&Interval<O, D>> + '_ {
        self.entries.iter()
            .map(|entry| &entry.interval)
    }

    pub fn insert<I: Into<Interval<O, D>>>(&mut self, interval: I) {
        let interval = interval.into();
        let max = interval.end;
        self.entries.push(InternalEntry {
            interval,
            max,
        });
        self.indexed = false;
    }

    pub fn insert_and_index<I: Into<Interval<O, D>>>(&mut self, interval: I) {
        let interval = interval.into();
        let max = interval.end;
        self.entries.push(InternalEntry {
            interval,
            max,
        });
        self.sort_and_index();
    }

    pub fn sort_and_index(&mut self) {
        if !self.indexed {
            self.entries.sort_unstable_by(|a, b| a.interval.cmp(&b.interval));
            self.index_core();
            self.indexed = true;
        }
    }

    fn index_core(&mut self) {
        let a = &mut self.entries;
        if a.is_empty() {
            return;
        }

        let n = a.len();
        let mut last_i = 0;
        let mut last_value = a[0].max;
        (0..n).step_by(2).for_each(|i| {
            last_i = i;
            a[i].max = a[i].interval.end;
            last_value = a[i].max;
        });
        let mut k = 1;
        while (1 << k) <= n {
            // process internal nodes in the bottom-up order
            let x = 1 << (k - 1);
            let i0 = (x << 1) - 1; // i0 is the first node
            let step = x << 2;
            for i in (i0..n).step_by(step) {
                // traverse all nodes at level k
                let end_left = a[i - x].max; // max value of the left child
                let end_right = if i + x < n { a[i + x].max } else { last_value }; // max value of the right child
                let end = max3(a[i].interval.end, end_left, end_right);
                a[i].max = end;
            }
            last_i = if (last_i >> k & 1) > 0 {
                last_i - x
            } else {
                last_i + x
            };
            if last_i < n && a[last_i].max > last_value {
                last_value = a[last_i].max
            }
            k += 1;
        }
        self.max_level = k - 1;
    }

    pub fn find(&self, start: O, end: O) -> IntervalOverlapIter<'_, O, D, S> {
        if !self.indexed {
            panic!("This interval tree has not been indexed yet. Call `index()` first.")
        }

        IntervalOverlapIter::new(self, start, end)
    }

    /// Check if a position is contained within one of the intervals
    pub fn contains(&self, pos: O) -> bool {
        self.find(pos, pos + O::one()).next().is_some()
    }

    #[inline]
    pub fn get(&self, index: usize) -> Option<&Interval<O, D>> {
        self.entries.get(index).map(|entry| &entry.interval)
    }
}

impl<O, D, const S: usize> IntoIterator for IntervalTree<O, D, S>
where
    O: OffsetType,
{
    type Item = Interval<O, D>;
    type IntoIter = Map<
        <SmallVec<[InternalEntry<O, D>; S]> as IntoIterator>::IntoIter,
        fn(InternalEntry<O, D>) -> Self::Item
    >;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter().map(|v| v.interval)
    }
}


impl<'a, O, D, const S: usize> IntoIterator for &'a IntervalTree<O, D, S>
where
    O: OffsetType + 'a,
    D: 'a
{
    type Item = &'a Interval<O, D>;
    type IntoIter = Map<
        <&'a SmallVec<[InternalEntry<O, D>; S]> as IntoIterator>::IntoIter,
        fn(&'a InternalEntry<O, D>) -> Self::Item
    >;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.iter().map(|v| &v.interval)
    }
}


fn max3<T: Ord>(a: T, b: T, c: T) -> T {
    a.max(b.max(c))
}

#[derive(Clone, Copy)]
struct StackCell {
    // node
    x: usize,
    // level
    k: usize,
    // false if left child hasn't been processed
    w: bool,
}

impl StackCell {
    fn empty() -> Self {
        Self {
            x: 0,
            k: 0,
            w: false,
        }
    }
}

pub struct IntervalOverlapIter<'a, O, D, const S: usize>
where
    O: OffsetType,
{
    tree: &'a IntervalTree<O, D, S>,
    start: O,
    end: O,
    stack: [StackCell; 64],
    t: usize,
    subtree: SmallVec<[&'a Interval<O, D>; 8]>,
    subtree_ix: usize
}

impl<'a, O, D, const S: usize> IntervalOverlapIter<'a, O, D, S>
where
    O: OffsetType,
{
    fn new(tree: &'a IntervalTree<O, D, S>, start: O, end: O) -> Self {
        let mut new_iter = Self {
            tree, start, end, stack: [StackCell::empty(); 64], t: 1,
            subtree: SmallVec::default(), subtree_ix: 0
        };

        // push the root; this is a top down traversal
        new_iter.stack[0].k = tree.max_level;
        new_iter.stack[0].x = (1 << tree.max_level) - 1;
        new_iter.stack[0].w = false;

        new_iter
    }
}

impl<'a, O, D, const S: usize> Iterator for IntervalOverlapIter<'a, O, D, S>
where
    O: OffsetType,
{
    type Item = &'a Interval<O, D>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.subtree.is_empty() && self.subtree_ix < self.subtree.len() {
            self.subtree_ix += 1;
            return Some(self.subtree[self.subtree_ix - 1]);
        } else if self.subtree_ix >= self.subtree.len() {
            self.subtree.clear();
            self.subtree_ix = 0;
        }

        let n = self.tree.len();
        let a = &self.tree.entries;
        while self.t > 0 {
            self.t -= 1;
            let StackCell { k, x, w } = self.stack[self.t];
            if k <= 3 {
                // we are in a small subtree; traverse every node in this subtree
                let i0 = x >> k << k;
                let i1 = min(i0 + (1 << (k + 1)) - 1, n);
                for node in a.iter().take(i1).skip(i0) {
                    if node.interval.start >= self.end {
                        break;
                    }

                    if self.start < node.interval.end {
                        self.subtree.push(&node.interval);
                    }
                }

                if let Some(ival) = self.subtree.first() {
                    self.subtree_ix = 1;
                    return Some(*ival);
                }
            } else if !w {
                // if left child not processed
                let y = x - (1 << (k - 1)); // the left child of x; NB: y may be out of range (i.e. y>=n)
                self.stack[self.t].k = k;
                self.stack[self.t].x = x;
                self.stack[self.t].w = true; // re-add node x, but mark the left child having been processed
                self.t += 1;
                if y >= n || a[y].max > self.start {
                    // push the left child if y is out of range or may overlap with the query
                    self.stack[self.t].k = k - 1;
                    self.stack[self.t].x = y;
                    self.stack[self.t].w = false;
                    self.t += 1;
                }
            } else if x < n && a[x].interval.start < self.end {
                // need to push the right child
                self.stack[self.t].k = k - 1;
                self.stack[self.t].x = x + (1 << (k - 1));
                self.stack[self.t].w = false;
                self.t += 1;

                if self.start < a[x].interval.end {
                    return Some(&a[x].interval);
                }
            }
        }

        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.tree.len()))
    }

}


#[cfg(test)]
mod tests {
    use super::{Interval, IntervalTree};

    type IntervalTest = Interval<u32, ()>;
    type IntervalTreeTest = IntervalTree<u32, ()>;

    #[test]
    fn test_interval_overlaps() {
        let ival1 = IntervalTest::new(3..5, ());
        let ival2 = IntervalTest::new( 7..10, ());
        let ival3 = IntervalTest::new(5..9, ());

        assert!(!ival1.overlaps(&ival2));
        assert!(!ival2.overlaps(&ival1));
        assert!(!ival1.overlaps(&ival3));
        assert!(!ival3.overlaps(&ival1));
        assert!(ival1.overlaps_or_adjacent(&ival3));
        assert!(ival3.overlaps_or_adjacent(&ival1));
        assert!(ival3.overlaps(&ival2));
        assert!(ival2.overlaps(&ival3));
    }

    #[test]
    fn test_contains() {
        let ival_tree = IntervalTreeTest::from_iter(vec![
            (3..5, ()),
            (10..11, ()),
            (11..13, ()),
            (20..30, ()),
        ]);

        assert!(ival_tree.contains(3));
        assert!(ival_tree.contains(4));
        assert!(ival_tree.contains(10));
        assert!(ival_tree.contains(11));
        assert!(ival_tree.contains(25));

        assert!(!ival_tree.contains(1));
        assert!(!ival_tree.contains(5));
        assert!(!ival_tree.contains(8));
        assert!(!ival_tree.contains(50));
    }

    #[test]
    fn test_find() {
        let ival_tree = IntervalTreeTest::from_iter(vec![
            (3..5, ()),
            (10..11, ()),
            (11..13, ()),
            (20..30, ()),
        ]);

        let found: Vec<_> = ival_tree.find(10, 25).collect();
        assert_eq!(found, vec![
            &Interval::from((10..11, ())),
            &Interval::from((11..13, ())),
            &Interval::from((20..30, ())),
        ]);

        let ival_tree2 = IntervalTreeTest::from_iter(vec![
            ((3..5), ()),
            ((10..11), ())
        ]);

        let found2: Vec<_> = ival_tree2.find(4, 5).collect();
        assert_eq!(found2, vec![&Interval::from((3..5, ()))]);

        let found3: Vec<_> = ival_tree2.find(10, 11).collect();
        assert_eq!(found3, vec![&Interval::from((10..11, ()))]);

        let found4: Vec<_> = ival_tree2.find(11, 12).collect();
        assert!(found4.is_empty());
    }
}
