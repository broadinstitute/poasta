use std::ops::Bound;

use astar::{AstarAlignableGraph, AstarState};
use config::AlignerConfig;
use cost_models::AlignmentCostModel;

use crate::graph::alignment::GraphAlignment;

pub mod config;
pub mod cost_models;
pub(crate) mod astar;
pub(crate) mod fr_points;

/// Enum representing the kind of alignment to perform
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    /// Perform global alignment of the query sequence to the graph
    Global,

    /// Allow free indels at the beginning or end (optionally up to a given maximum),
    /// e.g., for semi-global alignment.
    EndsFree {
        qry_free_begin: Bound<usize>,
        qry_free_end: Bound<usize>,
        graph_free_begin: Bound<usize>,
        graph_free_end: Bound<usize>,
    },
}

pub struct PoastaAligner<C> {
    config: C,
}

impl<C> PoastaAligner<C>
where C: AlignerConfig,
{
    pub fn new(config: C) -> Self {
        Self {
            config,
        }
    }
    
    pub fn align<D, G, CM, AS>(&self, graph: &G, seq: impl AsRef<[u8]>) -> GraphAlignment<G::NodeType>
        where G: AstarAlignableGraph,
            CM: AlignmentCostModel<AstarState<G::NodeType, D> = AS>,
            C: AlignerConfig<CostModel=CM>,
            AS: AstarState
    {
        self.align_u8::<D, G, CM, AS>(graph, seq.as_ref())
    }
    
    fn align_u8<D, G, CM, AS>(&self, graph: &G, seq: &[u8]) -> GraphAlignment<G::NodeType> 
        where G: AstarAlignableGraph,
              CM: AlignmentCostModel<AstarState<G::NodeType, D> = AS>,
              C: AlignerConfig<CostModel=CM>,
              AS: AstarState
    {
        let mut astar_state = AS::init();
        
        GraphAlignment::new()
    }
}
