pub mod offsets;
pub mod state;
pub mod scoring;
pub mod queue;
pub mod alignment;

use smallvec::SmallVec;
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::bubbles::index::BubbleIndex;
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, Score, StateGraph, StateGraphNode, StateGraphSuccessors};
use crate::aligner::scoring::AlignmentCosts;

use crate::debug::DebugOutputWriter;
pub use alignment::{AlignedPair, Alignment};
use crate::aligner::queue::AlignStateQueue;
use crate::bubbles::ReachedBubbleExits;
use crate::debug::messages::DebugOutputMessage;

pub struct PoastaAligner<'a, C>
where
    C: AlignmentCosts
{
    costs: C,
    debug_output: Option<&'a DebugOutputWriter>,
}

impl<'a, C> PoastaAligner<'a, C>
where
    C: AlignmentCosts
{
    pub fn new(costs: C) -> Self {
        Self {
            costs,
            debug_output: None,
        }
    }

    pub fn new_with_debug_output(costs: C, debug_writer: &'a DebugOutputWriter) -> Self {
        PoastaAligner {
            costs,
            debug_output: Some(debug_writer),
        }
    }

    pub fn align<O, G, Seq>(
        &mut self,
        graph: &G,
        bubble_index: &BubbleIndex<G::NodeIndex, O>,
        seq: &Seq
    ) -> (Score, Alignment<G::NodeIndex>)
    where
        O: OffsetType,
        G: AlignableGraph,
        Seq: AsRef<[u8]>
    {
        self.align_u8(graph, bubble_index, seq.as_ref())
    }

    fn align_u8<O, G, SG, S>(&mut self, graph: &G, bubble_index: &BubbleIndex<G::NodeIndex, O>, seq: &[u8]) -> (Score, Alignment<G::NodeIndex>)
    where
        O: OffsetType,
        G: AlignableGraph,
        C: AlignmentCosts<StateGraphType<G::NodeIndex, O>=SG>,
        SG: StateGraph<G::NodeIndex, O, StateNode=S>,
        S: StateGraphNode<G::NodeIndex, O>
    {
        let max_offset = O::max_value().as_usize();

        assert!(seq.len() - 1 < max_offset, "Sequence is too long for Offset integer type!");

        let mut state_graph = self.costs.new_state_graph::<G, O>(graph);

        let mut queue = AlignStateQueue::new();
        let start_state = S::new(graph.start_node(), O::zero(), AlignState::Start);
        queue.queue_state(start_state, 0);

        let reached_end;
        let mut bubble_exits_reached: ReachedBubbleExits<O> = vec![
            SmallVec::default(); graph.node_count_with_start()
        ];
        let mut score = Score::Score(0usize);
        'main: loop {
            let Some(front) = queue.pop_front() else {
                panic!("Empty queue?");
            };

            if front.is_empty() {
                score += 1usize;
                queue.next_score();
                continue;
            }

            // eprintln!("----- SCORE: {} ---------", score);

            for state in &front {
                // eprintln!("Current state: {:?}", state);

                if self.reached_end(graph, seq, state) {
                    reached_end = (score, state.clone());
                    // eprintln!("- reached end.");
                    break 'main;
                }

                if score > state_graph.get_score(state) {
                    // eprintln!("- Already found another lower scoring path to this state.");
                    continue;
                }

                if !self.can_improve_alignment(bubble_index, &bubble_exits_reached, state, score) {
                    // The exit of a bubble the current graph node is part of has already been reached,
                    // and no path from this state can improve on that score.
                    // eprintln!("- Ignoring state {:?} because it can't improve over the current bubble exit score.", state);
                    continue;
                }

                // Update reached bubble exits
                match state.state() {
                    AlignState::Match | AlignState::Mismatch => {
                        if bubble_index.is_exit(state.node()) {
                            bubble_exits_reached[state.node().index()]
                                .push((state.offset(), score))
                        }
                    },
                    _ => ()
                }

                // Process successor states
                let mut sucesssors = StateGraphSuccessors::new(graph, bubble_index, seq, score, state);
                while let Some(leaf) = sucesssors.queue_next(&mut state_graph, &mut queue, &mut bubble_exits_reached) {
                    // eprintln!("- EXTEND match leaf {:?}", leaf);

                    if self.reached_end(graph, seq, &leaf) {
                        reached_end = (score, leaf);
                        break 'main;
                    }
                }

                // Ensure extension of indels
                if let Some((ins_ext, score_delta)) = state_graph.extend_insertion(state, score) {
                    // eprintln!("- EXPAND queue INS_EXT {:?} (score: {}+{}={})", ins_ext, score, score_delta, score + score_delta);
                    queue.queue_state(ins_ext, score_delta)
                }

                for child in graph.successors(state.node()) {
                    if let Some((del_ext, score_delta)) = state_graph.extend_deletion(state, child, score) {
                        // eprintln!("- EXPAND queue DEL_EXT {:?} (score: {}+{}={})", del_ext, score, score_delta, score + score_delta);
                        queue.queue_state(del_ext, score_delta)
                    }
                }
            }

            if let Some(debug_out) = self.debug_output {
                match DebugOutputMessage::new_from_sg(&state_graph, score) {
                    Ok(msg) => debug_out.log(msg),
                    Err(e) => eprintln!("WARNING! could not output state graph to TSV: {}", e)
                }
            }

            queue.next_score();
            score += 1usize;
            // eprintln!();
        }

        if let Some(debug_out) = self.debug_output {
            match DebugOutputMessage::new_from_sg(&state_graph, score) {
                Ok(msg) => debug_out.log(msg),
                Err(e) => eprintln!("WARNING! could not output state graph to TSV: {}", e)
            }
        }

        let (end_score, end_state) = reached_end;
        let alignment = state_graph.backtrace(&end_state);

        // let dp_cells = ((seq.len() + 1) * graph.node_count()) * 3;
        // eprintln!("States explored: {:?}, {:.1}% of DP cells ({})", state_tree.num_nodes(), (state_tree.num_nodes() as f64 * 100.0) / dp_cells as f64, dp_cells);

        (end_score, alignment)
    }

    fn reached_end<G, S, N, O>(&self, graph: &G, seq: &[u8], state: &S) -> bool
    where
        G: AlignableGraph<NodeIndex=N>,
        N: NodeIndexType,
        O: OffsetType,
        S: StateGraphNode<N, O>
    {
        state.offset().as_usize() == seq.len() && graph.is_end(state.node())
    }

    fn can_improve_alignment<S, N, O>(
        &self,
        bubble_index: &BubbleIndex<N, O>,
        bubble_exits_reached: &ReachedBubbleExits<O>,
        state: &S,
        current_score: Score,
    ) -> bool
    where
        S: StateGraphNode<N, O>,
        N: NodeIndexType,
        O: OffsetType,
    {
        if !bubble_index.node_is_part_of_bubble(state.node()) {
            return true;
        }

        for bubble in bubble_index.get_node_bubbles(state.node()) {
            let reached = &bubble_exits_reached[bubble.bubble_exit.index()];
            if reached.is_empty() {
                continue;
            }

            for (reached_exit_offset, score) in reached.iter() {
                // eprintln!("State {:?} is part of a bubble with exit {:?}, reached at offset {:?} and score: {:?}",
                //     state, bubble.bubble_exit, reached_exit_offset, score);

                let target_offset = state.offset() + bubble.dist_to_exit;

                // Case 1: Could we reach the (bubble exit, reached_offset) with a lower score
                //         by opening or extending a deletion?
                if target_offset <= *reached_exit_offset {
                    let offset_diff = *reached_exit_offset - state.offset();
                    let vertical_gap = offset_diff - bubble.dist_to_exit;

                    let del_open_score_from_exit = self.costs.gap_cost(AlignState::Match, vertical_gap.as_usize());

                    if *score + del_open_score_from_exit <= current_score {
                        // eprintln!("- Current score: {}, score + gap: {} >= {}",
                        //           current_score, current_score + del_open_score_from_exit, score);
                        // eprintln!("- can't improve (vertical)");
                        return false;
                    }
                }

                // Case 2: Could we reach the (bubble exit, target offset) with a lower score
                //         compared to the case where we would open an insertion from the bubble exit?
                if target_offset > *reached_exit_offset && state.node() != bubble.bubble_exit {
                    let horizontal_gap = target_offset - *reached_exit_offset;

                    let ins_open_score_from_exit = self.costs.gap_cost(AlignState::Match, horizontal_gap.as_usize());
                    if *score + ins_open_score_from_exit <= current_score {
                        // eprintln!("- Bubble exit score: {}, bubble exit score + insertion: {} <= {}",
                        //           score, *score + ins_open_score_from_exit, current_score);
                        // eprintln!("- can't improve (horizontal)");
                        return false;
                    }
                }
            }
        }

        true
    }
}
