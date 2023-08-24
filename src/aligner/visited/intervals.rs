use std::cmp::Ordering;
use std::fmt::Debug;

use crate::aligner::offsets::{Diag, DiagonalPoint, OffsetType};
use crate::graphs::NodeIndexType;

use smallvec::SmallVec;
use serde::{Deserialize, Serialize};
use super::interval_tree::{Interval, IntervalTree};

pub type VisitedIntervalType<O> = Interval<O, VisitedIntervalData<O>>;
pub type VisitedIntervalTree<O> = IntervalTree<O, VisitedIntervalData<O>>;

const IVALS_ON_STACK: usize = 4;
pub type IntervalSmallVec<O> = SmallVec<[VisitedIntervalType<O>; IVALS_ON_STACK]>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AlignState {
    Insertion,
    Insertion2,
    Deletion,
    Deletion2,
    MisMatch,
    Extended,
    Start,
}

impl PartialOrd for AlignState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(match self {
            Self::Insertion | Self::Insertion2 => match other {
                Self::Insertion | Self::Insertion2 => Ordering::Equal,
                _ => Ordering::Less
            },
            Self::Deletion | Self::Deletion2 => match other {
                Self::Insertion | Self::Insertion2 => Ordering::Greater,
                Self::Deletion | Self::Deletion2 => Ordering::Equal,
                _ => Ordering::Less
            },
            Self::MisMatch | Self::Extended => match other {
                Self::MisMatch | Self::Extended => Ordering::Equal,
                Self::Start => Ordering::Less,
                _ => Ordering::Greater
            },
            Self::Start => match other {
                Self::Start => Ordering::Equal,
                _ => Ordering::Greater
            }
        })
    }
}

impl Ord for AlignState {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}


#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Backtrace {
    pub(crate) prev_diag: Option<Diag>,
    pub(crate) prev_state: AlignState,
}

impl Backtrace {
    pub fn new(diag: Diag, prev_state: AlignState) -> Self {
        Self {
            prev_diag: Some(diag),
            prev_state
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct VisitedIntervalData<O> {
    pub extendable: bool,
    pub clipped: O,
    pub(crate) bt: Backtrace
}

impl<O> VisitedIntervalData<O>
where
    O: OffsetType
{
    pub fn new(extendable: bool, clipped: O, bt: Backtrace) -> Self {
        Self {
            extendable,
            clipped,
            bt
        }
    }

    pub fn initial() -> Self {
        Self {
            extendable: true,
            clipped: O::zero(),
            bt: Backtrace { prev_diag: None, prev_state: AlignState::Start }
        }
    }
}


pub trait VisitedInterval<O> : Clone + Debug + Ord + Eq
where
    O: OffsetType,
{

    fn start(&self) -> O;
    fn end(&self) -> O;
    fn len(&self) -> O {
        self.end() - self.start()
    }

    fn is_empty(&self) -> bool {
        self.start() == self.end()
    }

    fn start_mut(&mut self) -> &mut O;
    fn end_mut(&mut self) -> &mut O;

    fn extendable(&self) -> bool;
    fn set_extendable(&mut self, value: bool);

    fn clipped(&self) -> O;
    fn clipped_mut(&mut self) -> &mut O;

    fn backtrace(&self) -> &Backtrace;
    fn backtrace_mut(&mut self) -> &mut Backtrace;

    fn overlaps(&self, other: &Self) -> bool;
}

impl<O> VisitedInterval<O> for VisitedIntervalType<O>
where
    O: OffsetType,
{
    #[inline(always)]
    fn start(&self) -> O {
        self.start
    }

    #[inline(always)]
    fn end(&self) -> O {
        self.end
    }

    #[inline(always)]
    fn start_mut(&mut self) -> &mut O {
        self.start_mut()
    }

    #[inline(always)]
    fn end_mut(&mut self) -> &mut O {
        self.end_mut()
    }

    #[inline(always)]
    fn extendable(&self) -> bool {
        self.data().extendable
    }

    #[inline(always)]
    fn set_extendable(&mut self, value: bool) {
        self.data_mut().extendable = value;
    }

    #[inline(always)]
    fn clipped(&self) -> O {
        self.data().clipped
    }

    #[inline(always)]
    fn clipped_mut(&mut self) -> &mut O {
        &mut self.data_mut().clipped
    }

    #[inline(always)]
    fn backtrace(&self) -> &Backtrace {
        &self.data().bt
    }

    #[inline(always)]
    fn backtrace_mut(&mut self) -> &mut Backtrace {
        &mut self.data_mut().bt
    }

    #[inline]
    fn overlaps(&self, other: &Self) -> bool {
        self.overlaps(other)
    }
}

pub trait VisitedIntervalOnDiag<O> {
    fn end_point_row<N: NodeIndexType>(&self) -> N;
}

impl<O, T> VisitedIntervalOnDiag<O> for (Diag, &T)
    where
        O: OffsetType,
        T: VisitedInterval<O>,
{
    #[inline(always)]
    fn end_point_row<N: NodeIndexType>(&self) -> N {
        (self.0, self.1.end() - O::one()).row()
    }
}
