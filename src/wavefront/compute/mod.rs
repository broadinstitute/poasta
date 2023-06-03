pub mod gap_linear;
pub mod gap_affine;

use crate::graph::POAGraph;
use crate::alignment::Alignment;
use crate::wavefront::fr_points::{OffsetPrimitive, FRPoint, ExtendCandidate};

pub trait WFCompute<Offset: OffsetPrimitive> : Default {
    fn reached_point(&self, point: &FRPoint<Offset>) -> bool;

    fn extend_candidates(&self) -> Vec<ExtendCandidate<Offset>>;
    fn extend(&mut self, candidate: &ExtendCandidate<Offset>);
    fn next<Seq: Eq + Clone>(&mut self, graph: &POAGraph<Seq>, new_score: i64);
    fn backtrace<Seq: Eq + Clone>(&self, graph: &POAGraph<Seq>, end_point: &FRPoint<Offset>) -> Alignment;
}
