pub mod gap_linear;
pub mod gap_affine;

use std::fmt::Display;
use crate::aligner::offsets::{Diag, OffsetType};
use crate::graphs::AlignableGraph;
use crate::aligner::visited::VisitedIntervalType;

pub use gap_affine::GapAffine;
use crate::aligner::Alignment;
use crate::debug::DebugOutputWriter;


pub trait AlignmentCosts: Copy {
    type LayerComputeType<'a, O>: LayerCompute<'a, O>
    where
        O: OffsetType;

    fn new_layer_compute<'a, G, O>(&self, graph: &'a G) -> Self::LayerComputeType<'a, O>
    where
        G: AlignableGraph,
        O: OffsetType;
}

pub trait AlignmentCostsEdit: AlignmentCosts {
    fn mismatch(&self) -> i64;
}

pub trait AlignmentCostsLinear: AlignmentCostsEdit {
    fn gap_extend(&self) -> i64;
}

pub trait AlignmentCostsAffine: AlignmentCostsLinear {
    fn gap_open(&self) -> i64;
}

pub trait AlignmentCostsTwoPiece: AlignmentCostsAffine {
    fn gap_extend2(&self) -> i64;
    fn gap_open2(&self) -> i64;
}

pub trait AlignmentStateCompute<N, O, Ix>: Display {

}

/// A trait that defines operations on the alignment state tree, and thus provides an interface
/// to generate new alignment states based on current states.
pub trait LayerCompute<'a, O>: Display
where
    O: OffsetType,
{
    fn with_debug(self, debug_writer: &'a DebugOutputWriter) -> Self;
    fn extend<G: AlignableGraph>(&mut self, graph: &G, seq: &[u8]) -> Option<(Diag, VisitedIntervalType<O>)>;
    fn next_layer<G: AlignableGraph>(&mut self, graph: &G, seq: &[u8], new_score: i64);
    fn backtrace<G: AlignableGraph>(&self, graph: &G, seq: &[u8], end_diag: Diag, end_interval: &VisitedIntervalType<O>) -> Alignment<G::NodeIndex>;
}
