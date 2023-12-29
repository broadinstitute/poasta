use std::rc::Rc;
use crate::aligner::astar::AstarVisited;
use crate::aligner::heuristic::{Dijkstra, AstarHeuristic, MinimumGapCostAffine};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, GapAffine};
use crate::aligner::scoring::gap_affine::AffineAstarData;
use crate::bubbles::index::BubbleIndex;
use crate::graphs::{AlignableRefGraph, NodeIndexType};

pub trait AlignmentConfig {
    type Costs: AlignmentCosts;
    type AstarData<N, O>: AstarVisited<N, O>
        where N: NodeIndexType,
              O: OffsetType;
    type Heuristic<N, O>: AstarHeuristic<N, O>;

    fn init_alignment<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType
    ) -> (
        <Self::Costs as AlignmentCosts>::AlignmentGraphType,
        Self::AstarData<G::NodeIndex, O>,
        Self::Heuristic<G::NodeIndex, O>
    )
        where G: AlignableRefGraph,
              O: OffsetType;

    fn init_alignment_with_existing_bubbles<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType,
        bubble_index: Rc<BubbleIndex<G::NodeIndex>>,
    ) -> (
        <Self::Costs as AlignmentCosts>::AlignmentGraphType,
        Self::AstarData<G::NodeIndex, O>,
        Self::Heuristic<G::NodeIndex, O>
    )
        where G: AlignableRefGraph,
              O: OffsetType;
}


pub struct AffineDijkstra(pub GapAffine);

impl AlignmentConfig for AffineDijkstra {
    type Costs = GapAffine;
    type AstarData<N, O> = AffineAstarData<N, O>
        where N: NodeIndexType,
              O: OffsetType;
    type Heuristic<N, O> = Dijkstra<N, O>;

    fn init_alignment<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType
    ) -> (
        <Self::Costs as AlignmentCosts>::AlignmentGraphType,
        Self::AstarData<G::NodeIndex, O>,
        Self::Heuristic<G::NodeIndex, O>
    )
        where G: AlignableRefGraph,
              O: OffsetType
    {
        let aln_graph = self.0.new_alignment_graph(aln_type);

        let bubble_index = Rc::new(BubbleIndex::new(ref_graph));

        let astar_data = AffineAstarData::new(self.0, ref_graph, seq, bubble_index);
        let heuristic = Dijkstra::default();

        (aln_graph, astar_data, heuristic)
    }

    fn init_alignment_with_existing_bubbles<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType,
        bubble_index: Rc<BubbleIndex<G::NodeIndex>>,
    ) -> (
        <Self::Costs as AlignmentCosts>::AlignmentGraphType,
        Self::AstarData<G::NodeIndex, O>,
        Self::Heuristic<G::NodeIndex, O>
    )
        where G: AlignableRefGraph,
              O: OffsetType
    {
        let aln_graph = self.0.new_alignment_graph(aln_type);

        let astar_data = AffineAstarData::new(self.0, ref_graph, seq, bubble_index);
        let heuristic = Dijkstra::default();

        (aln_graph, astar_data, heuristic)
    }
}

pub struct AffineMinGapCost(pub GapAffine);

impl AlignmentConfig for AffineMinGapCost {
    type Costs = GapAffine;
    type AstarData<N, O> = AffineAstarData<N, O>
        where N: NodeIndexType,
              O: OffsetType;
    type Heuristic<N, O> = MinimumGapCostAffine<N, O>;

    fn init_alignment<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType
    ) -> (
        <Self::Costs as AlignmentCosts>::AlignmentGraphType,
        Self::AstarData<G::NodeIndex, O>,
        Self::Heuristic<G::NodeIndex, O>
    )
        where G: AlignableRefGraph,
              O: OffsetType
    {
        let aln_graph = self.0.new_alignment_graph(aln_type);

        let bubble_index = Rc::new(BubbleIndex::new(ref_graph));

        let astar_data = AffineAstarData::new(self.0, ref_graph, seq, bubble_index.clone());
        let heuristic = MinimumGapCostAffine::new(self.0, bubble_index, seq.len());

        (aln_graph, astar_data, heuristic)
    }

    fn init_alignment_with_existing_bubbles<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType,
        bubble_index: Rc<BubbleIndex<G::NodeIndex>>,
    ) -> (
        <Self::Costs as AlignmentCosts>::AlignmentGraphType,
        Self::AstarData<G::NodeIndex, O>,
        Self::Heuristic<G::NodeIndex, O>
    )
        where G: AlignableRefGraph,
              O: OffsetType
    {
        let aln_graph = self.0.new_alignment_graph(aln_type);

        let astar_data = AffineAstarData::new(self.0, ref_graph, seq, bubble_index.clone());
        let heuristic = MinimumGapCostAffine::new(self.0, bubble_index.clone(), seq.len());

        (aln_graph, astar_data, heuristic)
    }
}