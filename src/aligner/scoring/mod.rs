// pub mod gap_linear;
pub mod gap_affine;

use crate::aligner::aln_graph::{AlignmentGraph, AlignState};

pub use gap_affine::GapAffine;
use std::default;
use std::ops::{Add, AddAssign, Sub, SubAssign, Bound};
use std::fmt::{Display, Formatter};
use std::cmp::Ordering;
use nonmax::NonMaxU32;
use crate::aligner::astar::AstarQueue;
use crate::aligner::offsets::OffsetType;
use crate::graphs::NodeIndexType;


pub trait AlignmentCosts: Copy {

    type AlignmentGraphType: AlignmentGraph;
    type QueueType<N, O>: AstarQueue<N, O>
        where N: NodeIndexType,
              O: OffsetType;

    fn new_alignment_graph(&self, aln_type: AlignmentType) -> Self::AlignmentGraphType;

    fn mismatch(&self) -> u8;
    
    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;
    
    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;

    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize;
}

pub trait GetAlignmentCosts {
    type Costs: AlignmentCosts;
    
    fn get_costs(&self) -> &Self::Costs;
}


/// Type of alignment to perform
///
/// Global alignment forces end-to-end alignment of the graph and the query
/// Ends-free alignment allows for free indels at the beginning of the query, end of the query,
/// beginning of the graph, or the end of the graph. The maximum
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum AlignmentType {
    /// Perform global alignment of the query sequence to the graph
    Global,

    /// Allow free indels at the beginning or end (optionally up to a given maximum)
    EndsFree {
        qry_free_begin: Bound<usize>,
        qry_free_end: Bound<usize>,
        graph_free_begin: Bound<usize>,
        graph_free_end: Bound<usize>
    }
}

#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub enum Score {
    Score(NonMaxU32),  // Use non-max, such that the maximum value can be used for Unvisited
    
    #[default]
    Unvisited
}

impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> Ordering {
        match self {
            Self::Score(score) => match other {
                Self::Score(other_score) => score.cmp(other_score),
                Self::Unvisited => Ordering::Less,
            },
            Self::Unvisited => match other {
                Self::Score(_) => Ordering::Greater,
                Self::Unvisited => Ordering::Equal
            }
        }
    }
}

impl Add<usize> for Score {
    type Output = Self;

    fn add(self, rhs: usize) -> Self::Output {
        match self {
            Self::Score(score) => Self::Score(NonMaxU32::new(u32::from(score) + rhs as u32).unwrap()),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl Add<u8> for Score {
    type Output = Self;

    fn add(self, rhs: u8) -> Self::Output {
        match self {
            Self::Score(score) => Self::Score(NonMaxU32::new(u32::from(score) + rhs as u32).unwrap()),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl AddAssign<u8> for Score {
    fn add_assign(&mut self, rhs: u8) {
        match self {
            Self::Score(score) => *score = NonMaxU32::new(u32::from(*score) + rhs as u32).unwrap(),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl Sub<usize> for Score {
    type Output = Self;

    fn sub(self, rhs: usize) -> Self::Output {
        match self {
            Self::Score(score) => Self::Score(NonMaxU32::new(u32::from(score) - rhs as u32).unwrap()),
            Self::Unvisited => panic!("Can't subtract from Score::Unvisited!")
        }
    }
}

impl SubAssign<usize> for Score {
    fn sub_assign(&mut self, rhs: usize) {
        match self {
            Self::Score(score) => *score = NonMaxU32::new(u32::from(*score) - rhs as u32).unwrap(),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl Sub<u8> for Score {
    type Output = Self;

    fn sub(self, rhs: u8) -> Self::Output {
        match self {
            Self::Score(score) => Self::Score(NonMaxU32::new(u32::from(score) - rhs as u32).unwrap()),
            Self::Unvisited => panic!("Can't subtract from Score::Unvisited!")
        }
    }
}

impl SubAssign<u8> for Score {
    fn sub_assign(&mut self, rhs: u8) {
        match self {
            Self::Score(score) => *score = NonMaxU32::new(u32::from(*score) - rhs as u32).unwrap(),
            Self::Unvisited => panic!("Can't add to Score::Unvisited!")
        }
    }
}

impl From<Score> for u32 {
    fn from(value: Score) -> Self {
        match value {
            Score::Score(score) => score.into(),
            Score::Unvisited => panic!("Trying to convert Score::Unvisited!")
        }
    }
}

impl Display for Score {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Score(score) => Display::fmt(score, f),
            Self::Unvisited => Display::fmt("unvisited", f)
        }
    }
}
