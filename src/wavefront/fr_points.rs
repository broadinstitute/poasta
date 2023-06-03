use std::cmp::{Ordering, Ord};
use std::collections::VecDeque;
use std::fmt::Debug;

use num::{FromPrimitive, PrimInt};

pub type DiagIx = i64;

pub trait OffsetPrimitive: FromPrimitive + PrimInt
    + TryInto<i64> + Default + Copy + Debug
{
}


impl OffsetPrimitive for u8 { }
impl OffsetPrimitive for u16 { }
impl OffsetPrimitive for u32 { }
impl OffsetPrimitive for u64 { }

pub type FRPoint<Offset> = (DiagIx, Offset);

pub trait WFPoint<Offset: OffsetPrimitive>: Debug {
    fn new(k: DiagIx, offset: Offset) -> Self;

    fn rank(&self) -> usize;
    fn offset(&self) -> Offset;
    fn diag(&self) -> DiagIx;

    fn move_one(&self) -> Self;
}

impl<Offset: OffsetPrimitive> WFPoint<Offset> for FRPoint<Offset> {
    fn new(k: DiagIx, offset: Offset) -> Self {
        (k, offset)
    }

    fn rank(&self) -> usize {
        match self.1.try_into() {
            Ok(v) => {
                (v - self.0) as usize
            },
            Err(_) => panic!("Could not convert offset to usize!")
        }
    }

    fn offset(&self) -> Offset {
        self.1
    }

    fn diag(&self) -> DiagIx {
        self.0
    }

    fn move_one(&self) -> Self {
        (self.0, self.1 + Offset::one())
    }
}

#[derive(Clone, Debug)]
pub enum PrevState {
    Start,
    Match(i64),
    Mismatch(i64),
    Deletion(i64),
    Insertion(i64)
}

impl PrevState {
    pub fn score(&self) -> i64 {
        match self {
            PrevState::Start => 0,
            PrevState::Match(score) => *score,
            PrevState::Mismatch(score) => *score,
            PrevState::Deletion(score) => *score,
            PrevState::Insertion(score) => *score
        }
    }
}

#[derive(Clone, Debug)]
pub struct PrevCandidate<Offset: OffsetPrimitive>(pub FRPoint<Offset>, pub PrevState);

impl<Offset: OffsetPrimitive> PrevCandidate<Offset> {
    pub fn rank(&self) -> usize {
        self.0.rank()
    }

    pub fn offset(&self) -> Offset {
        self.0.offset()
    }

    pub fn diag(&self) -> DiagIx {
        self.0.diag()
    }

    pub fn score(&self) -> i64 {
        self.1.score()
    }

    pub fn into_offset_with_bt(self, new_offset: Offset) -> OffsetWithBacktrace<Offset> {
        OffsetWithBacktrace {
            offset: new_offset,
            prev: Some(self)
        }
    }
}

impl<Offset: OffsetPrimitive> Ord for PrevCandidate<Offset> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.offset().cmp(&other.0.offset())
    }
}

impl<Offset: OffsetPrimitive> PartialOrd for PrevCandidate<Offset> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<Offset: OffsetPrimitive> PartialEq for PrevCandidate<Offset> {
    fn eq(&self, other: &Self) -> bool {
        self.0.offset() == other.0.offset()
    }
}

impl<Offset: OffsetPrimitive> Eq for PrevCandidate<Offset> {

}


#[derive(Clone, Debug)]
pub struct OffsetWithBacktrace<Offset: OffsetPrimitive> {
    pub offset: Offset,
    pub prev: Option<PrevCandidate<Offset>>
}

impl<Offset: OffsetPrimitive> Default for OffsetWithBacktrace<Offset> {
    fn default() -> Self {
        Self {
            offset: Offset::zero(),
            prev: Some(PrevCandidate((0, Offset::zero()), PrevState::Start))
        }
    }
}

impl<Offset: OffsetPrimitive> From<FRPoint<Offset>> for OffsetWithBacktrace<Offset> {
    fn from(value: FRPoint<Offset>) -> Self {
        Self {
            offset: value.offset(),
            prev: None
        }
    }
}

impl<Offset: OffsetPrimitive> Ord for OffsetWithBacktrace<Offset> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.offset.cmp(&other.offset)
    }
}

impl<Offset: OffsetPrimitive> PartialOrd for OffsetWithBacktrace<Offset> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<Offset: OffsetPrimitive> PartialEq for OffsetWithBacktrace<Offset> {
    fn eq(&self, other: &Self) -> bool {
        self.offset == other.offset
    }
}

impl<Offset: OffsetPrimitive> Eq for OffsetWithBacktrace<Offset> { }


#[derive(Clone, Debug)]
pub struct ExtendCandidate<Offset: OffsetPrimitive>(pub FRPoint<Offset>, pub Option<FRPoint<Offset>>);

impl<Offset: OffsetPrimitive> ExtendCandidate<Offset> {
    pub fn new_for_successor(&self, succ: FRPoint<Offset>) -> Self {
        Self(succ, Some(self.curr()))
    }

    pub fn curr(&self) -> FRPoint<Offset> {
        self.0
    }

    pub fn pred(&self) -> Option<FRPoint<Offset>> {
        self.1
    }

    pub fn make_offset_with_bt(&self, curr_score: i64) -> OffsetWithBacktrace<Offset> {
        OffsetWithBacktrace {
            offset: self.0.move_one().offset(),
            prev: self.1.and_then(
                // Only update backtrace if we switch diagonals
                |prev| if self.curr().diag() != prev.diag() {
                    eprintln!("{:?} -> {:?}, changed diagonals", self.curr(), self.pred());
                    Some(PrevCandidate(prev, PrevState::Match(curr_score)))
                } else {
                    None
                }
            )
        }
    }
}

pub type FRPointContainer<Offset> = VecDeque<Option<OffsetWithBacktrace<Offset>>>;
