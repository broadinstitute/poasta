pub mod gap_linear;
pub mod gap_affine;

use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use crate::graph::POAGraph;
use crate::alignment::Alignment;
use crate::wavefront::fr_points::{OffsetPrimitive, FRPoint, ExtendCandidate};

pub trait WFCompute: Default {
    type OffsetType: OffsetPrimitive;

    fn reached_end(&self, graph: &POAGraph, seq_length: usize) -> Option<FRPoint<Self::OffsetType>>;

    fn reset(&mut self);
    fn extend_candidates(&self) -> Vec<FRPoint<Self::OffsetType>>;
    fn extend(&mut self, candidate: &ExtendCandidate<Self::OffsetType>) -> bool;
    fn next(&mut self, graph: &POAGraph, seq_len: usize, new_score: i64);

    fn backtrace(&self, graph: &POAGraph, end_point: FRPoint<Self::OffsetType>) -> Alignment;

    fn write_csv(&self, writer: &mut BufWriter<File>) -> Result<(), Box<dyn Error>>;
}
