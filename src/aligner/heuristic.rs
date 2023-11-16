use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, GapAffine};
use crate::bubbles::index::BubbleIndex;
use crate::graphs::NodeIndexType;

pub trait AstarHeuristic {
    fn h<N, O>(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType;
}

pub trait WithBubbleIndex<N>
    where N: NodeIndexType,
{
    fn get_bubble_index(&self) -> &BubbleIndex<N>;
}

/// A* heuristic that always returns 0, such that
/// A* reduces to standard Dijkstra's algorithm.
#[derive(Default)]
pub struct Dijkstra;


impl AstarHeuristic for Dijkstra {
    fn h<N, O>(&self, _: &AlignmentGraphNode<N, O>, _: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType
    {
        0
    }
}


pub struct MinimumGapCostAffine {
    costs: GapAffine,
    dist_to_end: Vec<usize>,
    seq_length: usize,
}

impl MinimumGapCostAffine {
    pub fn new(costs: GapAffine, dist_to_end: Vec<usize>, seq_length: usize) -> Self {
        Self {
            costs,
            dist_to_end,
            seq_length,
        }
    }
}

impl AstarHeuristic for MinimumGapCostAffine {

    fn h<N, O>(&self, aln_node: &AlignmentGraphNode<N, O>, mut aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType,
    {
        let dist_to_exit =  self.dist_to_end[aln_node.node().index()].saturating_sub(1);
        let target_offset = aln_node.offset().as_usize() + dist_to_exit;

        let mut gap_length;
        if target_offset < self.seq_length {
            gap_length = self.seq_length - target_offset;

            // requires insertions, so if in deletion or match state, we need to open a new gap
            if aln_state != AlignState::Insertion {
                aln_state = AlignState::Match;
            }
        } else {
            gap_length = target_offset - self.seq_length;

            // requires deletions, so if in insertion or match state, we need to open a new gap
            if aln_state != AlignState::Deletion {
                aln_state = AlignState::Match;
            }
        }

        self.costs.gap_cost(aln_state, gap_length)
    }
}
