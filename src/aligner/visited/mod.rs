pub mod interval_tree;
pub mod intervals;
pub mod transformations;
pub mod merging;

pub use interval_tree::IntervalTree;
pub use intervals::{AlignState, Backtrace, IntervalSmallVec, VisitedIntervalTree,
                    VisitedInterval, VisitedIntervalType, VisitedIntervalData, VisitedIntervalOnDiag};
pub use merging::{SortedContainerMerger, merge_intervals};