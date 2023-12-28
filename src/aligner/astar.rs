use std::fmt::Write;
use crate::aligner::Alignment;
use crate::aligner::aln_graph::{AlignmentGraph, AlignmentGraphNode, AlignState};
use crate::aligner::config::AlignmentConfig;
use crate::aligner::dfa::{DepthFirstGreedyAlignment, ExtendResult};
use crate::aligner::heuristic::AstarHeuristic;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, Score};
use crate::bubbles::index::BubbleIndex;
use crate::debug::DebugOutputWriter;
use crate::debug::messages::DebugOutputMessage;
use crate::errors::PoastaError;
use crate::graphs::{AlignableRefGraph, NodeIndexType};


#[derive(Clone, Debug)]
pub struct AstarQueuedItem<N, O>(pub Score, pub AlignmentGraphNode<N, O>, pub AlignState)
    where N: NodeIndexType,
          O: OffsetType;

impl<N, O> AstarQueuedItem<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    #[inline]
    pub fn score(&self) -> Score {
        self.0
    }

    #[inline]
    pub fn aln_node(&self) -> AlignmentGraphNode<N, O> {
        self.1
    }

    #[inline]
    pub fn aln_state(&self) -> AlignState {
        self.2
    }
}

pub trait AstarQueue<N, O>: Default
    where N: NodeIndexType,
          O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>>;

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score, h: usize);
}

pub trait AstarVisited<N, O>
    where
        N: NodeIndexType,
        O: OffsetType,
{
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score;
    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score);

    fn visit(&mut self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState);

    fn prune(&self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool;

    fn update_score_if_lower(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
        parent: &AlignmentGraphNode<N, O>,
        parent_state: AlignState,
        score: Score,
    ) -> bool;

    fn backtrace<G>(&self, ref_graph: &G, aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
        where G: AlignableRefGraph<NodeIndex=N>;

    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), PoastaError>;
}

pub struct AstarResult<N>
    where N: NodeIndexType,
{
    pub score: Score,
    pub alignment: Alignment<N>,

    pub num_queued: usize,
    pub num_visited: usize,
    pub num_pruned: usize,
}

impl<N> Default for AstarResult<N>
    where N: NodeIndexType
{
    fn default() -> Self {
        Self {
            score: Score::Score(0),
            alignment: Vec::default(),

            num_queued: 0,
            num_visited: 0,
            num_pruned: 0
        }
    }
}


pub fn astar_alignment<O, C, Costs, AG, Q, G>(
    config: &C,
    ref_graph: &G,
    seq: &[u8],
    aln_type: AlignmentType,
    debug_writer: Option<&DebugOutputWriter>,
    existing_bubbles: Option<(BubbleIndex<G::NodeIndex>, Vec<(usize, usize)>)>
) -> AstarResult<G::NodeIndex>
    where O: OffsetType,
          C: AlignmentConfig<Costs=Costs>,
          Costs: AlignmentCosts<AlignmentGraphType=AG, QueueType<G::NodeIndex, O>=Q>,
          AG: AlignmentGraph,
          Q: AstarQueue<G::NodeIndex, O>,
          G: AlignableRefGraph,
{
    let (aln_graph, mut visited_data, heuristic) =
        if let Some((bubble_index, dist_to_end)) = existing_bubbles {
            config.init_alignment_with_existing_bubbles(ref_graph, seq, aln_type, bubble_index, dist_to_end)
        } else {
            config.init_alignment(ref_graph, seq, aln_type)
        };

    let mut result = AstarResult::default();
    let mut queue = Q::default();
    for initial_state in aln_graph.initial_states(ref_graph).into_iter() {
        let h = heuristic.h(&initial_state, AlignState::Match);
        queue.queue_aln_state(initial_state, AlignState::Match, Score::Score(0), h);
        visited_data.set_score(&initial_state, AlignState::Match, Score::Score(0));

        result.num_queued += 1;
    }

    let (end_score, end_node) = 'main: loop {
        let Some(AstarQueuedItem(score, aln_node, aln_state)) = queue.pop_aln_state() else {
            panic!("Could not align sequence! Empty queue before reaching end!")
        };

        if score > visited_data.get_score(&aln_node, aln_state) {
            continue;
        }

        eprintln!("---- FRONT ---- [Score: {score}] {aln_node:?} {aln_state:?}");
        eprintln!("- is end node? {:?} == {:?}", ref_graph.end_node(), aln_node.node());
        if aln_graph.is_end(ref_graph, seq, &aln_node, aln_state) {
            result.num_visited += 1;
            break (score, aln_node);
        }

        if visited_data.prune(score, &aln_node, aln_state) {
            // eprintln!("PRUNE {aln_node:?} ({aln_state:?}), score: {score:?}");
            result.num_pruned += 1;
            continue;
        }

        visited_data.visit(score, &aln_node, aln_state);
        result.num_visited += 1;

        // Perform depth-first greedy aligning of matches between graph and query
        if aln_state == AlignState::Match {
            let mut dfa = DepthFirstGreedyAlignment::new(ref_graph, seq, score, &aln_node);

            while let Some(end_point) = dfa.extend(&mut visited_data) {
                // eprintln!("DFA END: {end_point:?}");
                match end_point {
                    ExtendResult::RefGraphEnd(parent, child) => {
                        if aln_graph.is_end(ref_graph, seq, &child, AlignState::Match) {
                            break 'main (score, child);
                        }

                        // We can potentially still open an insertion from the parent, i.e., the
                        // last node in the graph that was not the end node.
                        aln_graph.expand_ref_graph_end(&mut visited_data, &parent, score,
                                                       |score_delta, succ, succ_state| {
                            let h = heuristic.h(&succ, succ_state);
                            result.num_queued += 1;
                            queue.queue_aln_state(succ, succ_state, score+score_delta, h);
                        });
                    },
                    ExtendResult::QueryEnd(parent, child) => {
                        // At query end, we can still open a deletion
                        aln_graph.expand_query_end(&mut visited_data, &parent, child, score,
                                                   |score_delta, succ, succ_state| {
                            let h = heuristic.h(&succ, succ_state);
                            result.num_queued += 1;
                            queue.queue_aln_state(succ, succ_state, score+score_delta, h);
                        })
                    },
                    ExtendResult::Mismatch(parent, child) => {
                        // In addition to mismatch state, also expand indel states
                        aln_graph.expand_mismatch(&mut visited_data, &parent, &child, score,
                                                  |score_delta, succ, succ_state| {
                            let h = heuristic.h(&succ, succ_state);
                            result.num_queued += 1;
                            queue.queue_aln_state(succ, succ_state, score+score_delta, h);
                        })
                    }
                }
            }

            result.num_visited += dfa.get_num_visited();
        } else {
            // For states other than the match state, don't use depth-first matching and simply
            // queue the neighbors.
            aln_graph.expand_all(&mut visited_data, ref_graph, seq, score, &aln_node, aln_state,
                                 |score_delta, succ, succ_state| {
                let h = heuristic.h(&succ, succ_state);
                result.num_queued += 1;
                queue.queue_aln_state(succ, succ_state, score+score_delta, h);
            })
        }

        if let Some(debug) = debug_writer {
            debug.log(DebugOutputMessage::new_from_astar_data(&visited_data));
        }
    };

    if let Some(debug) = debug_writer {
        debug.log(DebugOutputMessage::new_from_astar_data(&visited_data));
    }

    result.score = end_score;
    result.alignment = visited_data.backtrace(ref_graph, &end_node);

    result
}


