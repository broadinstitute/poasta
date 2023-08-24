use std::cmp::Ordering;
use std::collections::binary_heap::{BinaryHeap, PeekMut};

use crate::aligner::OffsetType;
use super::intervals::VisitedInterval;



/// An iterator that returns elements from one or more, pre-sorted, input iterators, and yields
/// elements of those iterators in merged, sorted order.
///
/// This iterator does not check if the input iterators are sorted, and so it is an error
/// if they are not.
///
/// Algorithm adapted from
/// https://dev.to/creativcoder/merge-k-sorted-arrays-in-rust-1b2f
pub struct SortedContainerMerger<I, V> {
    candidates: BinaryHeap<MergeCandidate<I, V>>,
    known_size: Option<usize>
}

impl<I, V> SortedContainerMerger<I, V> {
    pub fn new<C>(containers: &[C]) -> Self
        where
            I: Iterator<Item=V>,
            V: Ord + Eq,
            C: IntoIterator<IntoIter=I, Item=V> + Copy,
    {
        let mut candidates = BinaryHeap::with_capacity(containers.len());
        for container in containers.into_iter() {
            let mut c_iter = container.into_iter();

            if let Some(curr) = c_iter.next() {
                candidates.push(MergeCandidate { curr, iter: c_iter });
            }
        }

        Self { candidates, known_size: None }
    }

    pub fn new_exact_size<C>(containers: &[C]) -> Self
        where
            I: ExactSizeIterator<Item=V>,
            V: Ord + Eq,
            C: IntoIterator<IntoIter=I, Item=V> + Copy,
    {
        let mut candidates = BinaryHeap::with_capacity(containers.len());
        let mut total_length = 0;
        for container in containers {
            let mut c_iter = container.into_iter();
            total_length += c_iter.len();

            if let Some(curr) = c_iter.next() {
                candidates.push(MergeCandidate { curr, iter: c_iter });
            }
        }

        Self { candidates, known_size: Some(total_length) }
    }
}

impl<I, V> Iterator for SortedContainerMerger<I, V>
    where
        I: Iterator<Item=V>,
        V: Copy + Ord + Eq
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(mut candidate) = self.candidates.peek_mut() {
            let to_return = candidate.curr;

            if let Some(next) = candidate.iter.next() {
                candidate.curr = next;
            } else {
                PeekMut::pop(candidate);
            }

            return Some(to_return);
        }

        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, self.known_size)
    }
}

struct MergeCandidate<I, V> {
    curr: V,
    iter: I,
}

impl<I, V> PartialEq for MergeCandidate<I, V>
    where
        V: PartialEq
{
    fn eq(&self, other: &Self) -> bool {
        self.curr == other.curr
    }
}

impl<I, V> Eq for MergeCandidate<I, V>
    where
        V: Eq
{ }

impl<I, V> PartialOrd for MergeCandidate<I, V>
    where
        V: PartialOrd + PartialEq,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.curr.partial_cmp(&self.curr)
    }
}

impl<I, V> Ord for MergeCandidate<I, V>
    where
        V: Ord + Eq
{
    fn cmp(&self, other: &Self) -> Ordering {
        other.curr.cmp(&self.curr)
    }
}

pub fn merge_intervals<'a, C, I, V, O>(source_ivals: &'a [C]) -> Vec<V>
    where
        C: IntoIterator<IntoIter=I, Item=&'a V> + Copy,
        I: ExactSizeIterator<Item=&'a V>,
        V: VisitedInterval<O> + 'a,
        O: OffsetType + 'a,
{
    let mut merge_iter = SortedContainerMerger::new_exact_size(source_ivals);
    let num = merge_iter.size_hint().1.unwrap();
    let mut new_intervals: Vec<V> = Vec::with_capacity(num);
    let mut ival_to_add: Option<V> = merge_iter.next().cloned();

    for incoming_ival in merge_iter {
        if let Some(ref mut to_add) = ival_to_add {
            eprintln!("      Curr to add: {}{} ({}, {:?})",
                      " ".repeat(to_add.start().as_usize()),
                      "-".repeat(to_add.len().as_usize()),
                      to_add.extendable(),
                      to_add.backtrace().prev_state);
            eprintln!("Incoming interval: {}{} ({}, {:?})",
                      " ".repeat(incoming_ival.start().as_usize()),
                      "-".repeat(incoming_ival.len().as_usize()),
                      incoming_ival.extendable(),
                      incoming_ival.backtrace().prev_state);

            if incoming_ival.end() <= to_add.end() {
                // New interval is contained within previous, ignore
                eprintln!("Contained.");
                continue;
            }

            if to_add.overlaps(incoming_ival) {
                eprintln!("Overlaps.");

                // We already checked if the incoming interval was contained, so there should be at
                // least some new suffix.
                let mut curr_prefix = std::mem::replace(to_add, incoming_ival.clone());

                // `to_add` is now the new incoming interval
                *to_add.start_mut() = std::cmp::max(curr_prefix.start(), incoming_ival.start());

                if curr_prefix.start() < incoming_ival.start() {
                    *curr_prefix.end_mut() = incoming_ival.start();
                    curr_prefix.set_extendable(false);

                    let new_clipped = incoming_ival.clipped() + curr_prefix.len();
                    // Only set new clipped if the clipped part extends beyond the previous interval
                    if new_clipped > curr_prefix.len() + curr_prefix.clipped() {
                        *to_add.clipped_mut() = new_clipped;
                    }

                    new_intervals.push(curr_prefix);
                }
            } else {
                eprintln!("Non-overlap.");
                let new = std::mem::replace(to_add, incoming_ival.clone());
                new_intervals.push(new);
            }
        } else {
            ival_to_add = Some(incoming_ival.clone());
        }
    }

    if let Some(to_add) = ival_to_add {
        new_intervals.push(to_add);
    }

    eprintln!("All merged intervals:");
    for ival in new_intervals.iter() {
        eprintln!("{:?}", ival);
        eprintln!("                   {}{} ({}, {:?})",
                  " ".repeat(ival.start().as_usize()),
                  "-".repeat(ival.len().as_usize()),
                  ival.extendable(),
                  ival.backtrace().prev_state);
    }

    new_intervals
}

#[cfg(test)]
mod tests {
    use crate::aligner::visited::{IntervalSmallVec, merge_intervals, VisitedInterval, VisitedIntervalData, VisitedIntervalType};
    use super::SortedContainerMerger;

    type Ival = VisitedIntervalType<u32>;
    type IvalSmallVec = IntervalSmallVec<u32>;

    #[test]
    fn test_merge() {
        let c1 = vec![2, 8, 10];
        let c2 = vec![1, 4, 9];
        let c3 = vec![3, 4, 11];

        let containers = [&c1, &c2, &c3];
        let merge_iter = SortedContainerMerger::new_exact_size(&containers);

        let merged: Vec<_> = merge_iter.copied().collect();
        assert_eq!(merged, vec![1, 2, 3, 4, 4, 8, 9, 10, 11]);
    }

    #[test]
    fn test_sort_intervals() {
        let c1 = vec![
            Ival::new(3..10, VisitedIntervalData::initial()),
            Ival::new(40..45, VisitedIntervalData::initial()),
        ];

        let c2 = vec![
            Ival::new(1..25, VisitedIntervalData::initial()),
            Ival::new(25..30, VisitedIntervalData::initial()),
            Ival::new(50..60, VisitedIntervalData::initial()),
        ];

        let c3 = vec![
            Ival::new(11..18, VisitedIntervalData::initial()),
            Ival::new(20..24, VisitedIntervalData::initial()),
            Ival::new(26..29, VisitedIntervalData::initial()),
            Ival::new(50..70, VisitedIntervalData::initial()),
        ];

        let containers = [&c1, &c2, &c3];
        let merge_iter = SortedContainerMerger::new_exact_size(&containers);
        let merged: Vec<_> = merge_iter.collect();

        assert_eq!(merged, vec![&c2[0], &c1[0], &c3[0], &c3[1], &c2[1], &c3[2], &c1[1], &c3[3], &c2[2]])
    }

    #[test]
    fn test_merge_intervals() {
        let c1 = vec![
            Ival::new(3..10, VisitedIntervalData::initial()),
            Ival::new(40..45, VisitedIntervalData::initial()),
        ];

        let c2 = vec![
            Ival::new(1..25, VisitedIntervalData::initial()),
            Ival::new(25..30, VisitedIntervalData::initial()),
            Ival::new(42..48, VisitedIntervalData::initial()),
            Ival::new(50..60, VisitedIntervalData::initial()),
        ];

        let c3 = vec![
            Ival::new(11..18, VisitedIntervalData::initial()),
            Ival::new(20..24, VisitedIntervalData::initial()),
            Ival::new(26..29, VisitedIntervalData::initial()),
            Ival::new(50..70, VisitedIntervalData::initial()),
        ];

        let containers = [&c1, &c2, &c3];
        let merged_tree = merge_intervals(&containers);

        let truth_data = vec![
            (1..25, true),
            (25..30, true),
            (40..42, false),
            (42..48, true),
            (50..70, true)
        ];

        for (ival, (true_range, extendable)) in merged_tree.iter().zip(truth_data) {
            eprintln!("{:?}", ival);
            assert_eq!(ival.start, true_range.start);
            assert_eq!(ival.end, true_range.end);
            assert_eq!(ival.extendable(), extendable);
        }
    }
}
