use crate::aligner::astar::AstarVisited;
use crate::aligner::heuristic::{Dijkstra, AstarHeuristic, MinimumGapCostAffine};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, GapAffine};
use crate::aligner::scoring::gap_affine::AffineAstarData;
use crate::bubbles::index::BubbleIndexBuilder;
use crate::graphs::{AlignableRefGraph, NodeIndexType};

pub trait AlignmentConfig {
    type Costs: AlignmentCosts;
    type AstarData<N, O>: AstarVisited<N, O>
        where N: NodeIndexType,
              O: OffsetType;
    type Heuristic: AstarHeuristic;


    fn init_alignment<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType
    ) -> (<Self::Costs as AlignmentCosts>::AlignmentGraphType, Self::AstarData<G::NodeIndex, O>, Self::Heuristic)
        where G: AlignableRefGraph,
              O: OffsetType;
}


pub struct AffineDijkstra(pub GapAffine);

impl AlignmentConfig for AffineDijkstra {
    type Costs = GapAffine;
    type AstarData<N, O> = AffineAstarData<N, O>
        where N: NodeIndexType,
              O: OffsetType;
    type Heuristic = Dijkstra;

    fn init_alignment<O, G>(
        &self,
        ref_graph: &G,
        _: &[u8],
        aln_type: AlignmentType
    ) -> (<Self::Costs as AlignmentCosts>::AlignmentGraphType, Self::AstarData<G::NodeIndex, O>, Self::Heuristic)
        where G: AlignableRefGraph,
              O: OffsetType
    {
        let aln_graph = self.0.new_alignment_graph(aln_type);

        let (bubble_index, _) = BubbleIndexBuilder::new(ref_graph)
            .build();

        let astar_data = AffineAstarData::new(self.0, ref_graph, bubble_index);
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
    type Heuristic = MinimumGapCostAffine;

    fn init_alignment<O, G>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_type: AlignmentType
    ) -> (<Self::Costs as AlignmentCosts>::AlignmentGraphType, Self::AstarData<G::NodeIndex, O>, Self::Heuristic)
        where G: AlignableRefGraph,
              O: OffsetType
    {
        let aln_graph = self.0.new_alignment_graph(aln_type);

        let (bubble_index, dist_to_end) = BubbleIndexBuilder::new(ref_graph)
            .build();

        let astar_data = AffineAstarData::new(self.0, ref_graph, bubble_index);
        let heuristic = MinimumGapCostAffine::new(self.0, dist_to_end, seq.len());

        (aln_graph, astar_data, heuristic)
    }
}