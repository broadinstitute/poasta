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
use std::rc::Rc;
use crate::graphs::AlignableRefGraph;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentType;

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

    pub fn align<O, G, Seq>(
        &self,
        ref_graph: &G,
        seq: &Seq,
    ) -> AstarResult<G::NodeIndex>
        where O: OffsetType,
              G: AlignableRefGraph,
              Seq: AsRef<[u8]>,
    {
        self.align_u8::<O, _>(ref_graph, seq.as_ref(), None, true)
    }

    pub fn align_with_existing_bubbles<O, G, Seq>(
        &self,
        ref_graph: &G,
        seq: &Seq,
        existing_bubbles: Rc<BubbleIndex<G::NodeIndex>>,
    ) -> AstarResult<G::NodeIndex>
        where O: OffsetType,
              G: AlignableRefGraph,
              Seq: AsRef<[u8]>,
    {
        self.align_u8::<O, _>(ref_graph, seq.as_ref(), Some(existing_bubbles), true)
    }
    
    pub fn align_no_pruning<O, G, Seq>(
        &self,
        ref_graph: &G,
        seq: &Seq,
    ) -> AstarResult<G::NodeIndex>
        where O: OffsetType,
              G: AlignableRefGraph,
              Seq: AsRef<[u8]>,
    {
        self.align_u8::<O, _>(ref_graph, seq.as_ref(), None, false)
    }

    fn align_u8<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        existing_bubbles: Option<Rc<BubbleIndex<G::NodeIndex>>>,
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
}
