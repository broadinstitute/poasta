use super::offsets::OffsetType;
use crate::graphs::{AlignableRefGraph, NodeIndexType};
use std::fmt::Debug;
use crate::aligner::astar::AstarVisited;
use crate::aligner::scoring::{AlignmentCosts, Score};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignState {
    Match,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
}


/// A struct that represents a node in the alignment graph
///
/// It holds cursors to a node in the reference graph to which we are aligning,
/// and a query offset.
///
/// This struct does not store the alignment state (e.g., whether it's in Match,
/// insertion or deletion state). It is up to the specific alignment graph
/// implementations to retain that context.
///
/// ### See also
///
/// - [`AlignmentGraph`]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct AlignmentGraphNode<N, O>
where
    N: Copy + Clone + Debug + Eq,
    O: Copy + Clone + Debug + Eq,
{
    node: N,
    offset: O,
}


impl<N, O> AlignmentGraphNode<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{

    pub fn new(node: N, offset: O) -> Self {
        Self {
            node,
            offset,
        }
    }

    #[inline(always)]
    pub fn node(&self) -> N {
        self.node
    }

    #[inline(always)]
    pub fn offset(&self) -> O {
        self.offset
    }
}


pub trait AlignmentGraph {
    type CostModel: AlignmentCosts;

    fn get_costs(&self) -> &Self::CostModel;

    fn initial_states<G, O>(&self, ref_graph: &G) -> Vec<AlignmentGraphNode<G::NodeIndex, O>>
        where G: AlignableRefGraph,
              O: OffsetType;

    fn is_end<G, O>(&self, ref_graph: &G, seq: &[u8], node: &AlignmentGraphNode<G::NodeIndex, O>, aln_state: AlignState) -> bool
        where G: AlignableRefGraph,
              O: OffsetType;

    fn expand_all<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        state: AlignState,
        f: F
    )
        where V: AstarVisited<G::NodeIndex, O>,
              G: AlignableRefGraph,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState);

    fn expand_ref_graph_end<V, N, O, F>(&self, visited_data: &mut V, parent: &AlignmentGraphNode<N, O>, score: Score, f: F)
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);

    fn expand_query_end<V, N, O, F>(&self, visited_data: &mut V, parent: &AlignmentGraphNode<N, O>, child: N, score: Score, f: F)
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);

    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: Score,
        f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);
}
