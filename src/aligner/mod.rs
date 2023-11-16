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
use crate::graphs::AlignableRefGraph;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentType;

pub use alignment::{AlignedPair, Alignment};
use scoring::Score;
use crate::aligner::astar::astar_alignment;
use crate::aligner::config::AlignmentConfig;


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


pub struct PoastaAligner<C>
where
    C: AlignmentConfig,
{
    config: C,
    aln_type: AlignmentType,
}


impl<C> PoastaAligner<C>
    where C: AlignmentConfig,
{
    pub fn new(config: C, aln_type: AlignmentType) -> Self {
        Self {
            config,
            aln_type,
        }
    }

    pub fn align<'a, O, G, Seq>(
        &self,
        ref_graph: &'a G,
        seq: &'a Seq,
    ) -> (Score, Alignment<G::NodeIndex>)
        where O: OffsetType,
              G: AlignableRefGraph,
              Seq: AsRef<[u8]> + 'a,
    {
        self.align_u8::<O, _>(ref_graph, seq.as_ref())
    }

    fn align_u8<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
    ) -> (Score, Alignment<G::NodeIndex>)
    where
        O: OffsetType,
        G: AlignableRefGraph,
    {
        astar_alignment::<O, _, _, _, _, _>(
            &self.config,
            ref_graph,
            seq,
            self.aln_type,
        )
    }
}
