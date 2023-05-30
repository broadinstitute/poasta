pub mod gap_linear;
pub mod gap_affine;

use crate::graph::POAGraph;
use crate::wavefront::{OffsetPrimitive, FRPoint};

pub trait WFCompute<Offset: OffsetPrimitive> : Default {
    fn reached_point(&self, point: &FRPoint<Offset>) -> bool;

    fn extend_candidates(&self) -> Vec<FRPoint<Offset>>;
    fn extend(&mut self, point: &FRPoint<Offset>);
    fn next<Seq: Eq + Clone>(&mut self, graph: &POAGraph<Seq>, new_score: i64);
}
