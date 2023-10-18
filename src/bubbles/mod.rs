use smallvec::SmallVec;
use crate::aligner::state::Score;

pub mod finder;
pub mod index;

pub type ReachedBubbleExits<O> = Vec<SmallVec<[(O, Score); 4]>>;
