pub mod offsets;
pub mod aln_graph;
pub mod scoring;
pub mod queue;
pub mod alignment;
mod dfa;
pub mod astar;
pub mod heuristic;
pub mod config;

use std::ops::Bound;
use std::sync::Arc;
use nonmax::NonMaxU32;
use crate::graphs::AlignableRefGraph;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentType, Score};

pub use alignment::{AlignedPair, Alignment};
use crate::aligner::astar::{astar_alignment, AstarResult};
use crate::aligner::config::AlignmentConfig;
use crate::bubbles::index::BubbleIndex;
use crate::debug::DebugOutputWriter;


pub enum AlignmentMode<O> {
    /// Perform global alignment of the query sequence to the graph
    Global,

    /// Allow free indels at the beginning or end (optionally up to a given maximum)
    EndsFree {
        qry_free_begin: Bound<O>,
        qry_free_end: Bound<O>,
        graph_free_begin: Bound<O>,
        graph_free_end: Bound<O>
    }
}


pub struct PoastaAligner<'a, C>
where
    C: AlignmentConfig,
{
    config: C,
    aln_type: AlignmentType,
    debug_writer: Option<&'a DebugOutputWriter>
}


impl<'a, C> PoastaAligner<'a, C>
    where C: AlignmentConfig,
{
    pub fn new(config: C, aln_type: AlignmentType) -> Self {
        Self {
            config,
            aln_type,
            debug_writer: None
        }
    }

    pub fn new_with_debug(config: C, aln_type: AlignmentType, debug_writer: &'a DebugOutputWriter) -> Self {
        Self {
            config,
            aln_type,
            debug_writer: Some(debug_writer)
        }
    }

    pub fn align_with_existing_bubbles<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        existing_bubbles: Arc<BubbleIndex<G::NodeIndex>>,
    ) -> AstarResult<G::NodeIndex>
        where O: OffsetType,
              G: AlignableRefGraph,
    {
        self.align_internal::<O, _>(ref_graph, seq, Some(existing_bubbles), true)
    }
    
    pub fn align_no_pruning<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
    ) -> AstarResult<G::NodeIndex>
        where O: OffsetType,
              G: AlignableRefGraph,
    {
        self.align_internal::<O, _>(ref_graph, seq, None, false)
    }

    fn align_internal<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        existing_bubbles: Option<Arc<BubbleIndex<G::NodeIndex>>>,
        enable_pruning: bool,
    ) -> AstarResult<G::NodeIndex>
    where
        O: OffsetType,
        G: AlignableRefGraph,
    {
        astar_alignment::<O, _, _, _, _, _>(
            &self.config,
            ref_graph,
            seq,
            self.aln_type,
            self.debug_writer,
            existing_bubbles,
            enable_pruning,
        )
    }
    
    pub fn align<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
    ) -> AstarResult<G::NodeIndex>
    where
        O: OffsetType,
        G: AlignableRefGraph,
    {
        // Handle edge case: empty graph (only start/end nodes)
        if ref_graph.node_count() == 0 {
            // For empty graph, return a result with appropriate score
            // In ends-free alignment, aligning against empty reference should be possible
            let score = if seq.is_empty() {
                Score::Score(NonMaxU32::new(0).unwrap()) // Perfect match: empty query against empty reference
            } else {
                // Query has content but reference is empty - treat as all insertions
                // For ends-free alignment, this should have some cost but be valid
                Score::Score(NonMaxU32::new(seq.len() as u32 * 4).unwrap()) // Rough cost estimate
            };
            
            return AstarResult {
                score,
                alignment: Alignment::new(), // Empty alignment
                num_queued: 0,
                num_visited: 0,
                num_pruned: 0,
            };
        }
        
        self.align_internal::<O, _>(ref_graph, seq, None, true)
    }
}
