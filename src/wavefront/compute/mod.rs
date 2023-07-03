pub mod gap_linear;
pub mod gap_affine;

use crate::graph::POAGraph;
use crate::alignment::Alignment;
use crate::debug::DebugOutputWriter;
use crate::wavefront::offsets::{OffsetPrimitive, OffsetContainer};

pub trait WFCompute: Default {
    type OffsetType: OffsetPrimitive;

    fn furthest_offsets(&self) -> Option<&OffsetContainer<Self::OffsetType>>;
    fn reached_end(&self, graph: &POAGraph, seq_length: usize) -> Option<usize>;
    fn reset(&mut self);

    fn extend_candidates(&self) -> Vec<(usize, Self::OffsetType)>;
    fn is_further(&self, node: usize, offset: usize) -> bool;
    fn update_extended_path(&mut self, start_node: usize, path: Vec<(usize, usize)>);
    fn next(&mut self, graph: &POAGraph, seq_len: usize, new_score: i64);

    fn backtrace(&self, graph: &POAGraph, sequence: &[u8], end_node: usize) -> Alignment;

    fn log_debug_data(&self, debug: &DebugOutputWriter);
}
