pub mod offsets;
pub mod state;
pub mod extend;
pub mod scoring;
pub mod queue;
pub mod alignment;
pub mod visited;

use std::cell::RefCell;
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::state::{AlignState, StateTreeNode, Backtrace, TreeIndexType};
use crate::aligner::scoring::AlignmentCosts;
use crate::aligner::queue::{AlignStateQueue, PathShift, VisitedPath};
use crate::aligner::extend::PathExtender;

use crate::debug::DebugOutputWriter;
use crate::debug::messages::DebugOutputMessage;

pub use alignment::{AlignedPair, Alignment};

enum ExtendResult<Ix: TreeIndexType> {
    NewExtendedNodes(Vec<Ix>),
    ReachedEnd(Ix)
}

use ExtendResult::{NewExtendedNodes, ReachedEnd};
use scoring::AlignmentStateTree;


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

    pub fn align<O, Ix, G, S, N>(&mut self, graph: &G, sequence: &S) -> (usize, Alignment<N>)
    where
        O: OffsetType,
        Ix: TreeIndexType,
        G: AlignableGraph<NodeIndex=N>,
        S: AsRef<[u8]>,
        N: NodeIndexType,
    {
        let seq = sequence.as_ref();
        let max_offset = O::max_value().as_usize();

        assert!(seq.len() - 1 < max_offset, "Sequence is too long for Offset integer type!");

        let mut queue = AlignStateQueue::new();
        let mut state_tree: <C as AlignmentCosts>::StateTreeType<N, O, Ix> = self.costs.to_new_state_tree(graph);

        // Add graph start nodes to queue
        for start_node in graph.start_nodes().iter() {
            let start_state = StateTreeNode::new_start(*start_node);
            let new_ix = state_tree.add_node(start_state);
            queue.queue_endpoint(0, new_ix);
        }

        let mut score = 0;
        let reached_end_state;
        'main: loop {
            let Some(mut current) = queue.pop_current() else {
                panic!("Empty queue?")
            };

            if current.is_empty() {
                score += 1;
                continue;
            }

            // eprintln!("MARK VISITED score: {score}");
            self.mark_paths_as_visited(&mut queue, graph, seq, &mut state_tree, current.visited_paths());

            // Close indels for current score, and add to current queue
            // eprintln!("CLOSE INDELS score: {score}");
            let new_states = state_tree.close_indels_for(current.endpoints());
            current.queue_additional(new_states);

            // Try to extend the alignment along matching sequence in the graph
            // eprintln!("EXTEND score: {score}");
            match self.extend(graph, seq, &mut state_tree, current.endpoints_mut(), &mut queue) {
                ReachedEnd(end) => {
                    reached_end_state = (end, score);
                    break;
                },
                NewExtendedNodes(additional_nodes) => current.queue_additional(additional_nodes)
            }

            // If the end not reached yet, expand into next alignment states, including mismatches
            // and indels. New states to explore are queued per score, such that lower scores are
            // explored first.
            // eprintln!("EXPAND score: {score}");
            for state_ix in current.endpoints() {
                if let Some(end) = state_tree.generate_next(&mut queue, graph, seq.len(), *state_ix) {
                    reached_end_state = (end, score + self.costs.mismatch() as usize);
                    break 'main;
                }
            }

            score += 1;
        }

        let (end_state, end_score) = reached_end_state;
        let alignment = self.backtrace(&state_tree, end_state);

        if let Some(debug) = self.debug_output {
            debug.log(DebugOutputMessage::new_from_state_tree(&state_tree));
        }

        let dp_cells = ((seq.len() + 1) * graph.node_count()) * 3;
        eprintln!("States explored: {:?}, {:.1}% of DP cells ({})", state_tree.num_nodes(), (state_tree.num_nodes() as f64 * 100.0) / dp_cells as f64, dp_cells);

        (end_score, alignment)
    }

    fn mark_paths_as_visited<O, Ix, G, N, T>(&mut self, queue: &mut AlignStateQueue<N, O, Ix>,
                                             graph: &G, seq: &[u8], tree: &mut T, paths: &[VisitedPath<N, O>])
    where
        O: OffsetType,
        Ix: TreeIndexType,
        G: AlignableGraph<NodeIndex=N>,
        N: NodeIndexType,
        T: AlignmentStateTree<N, O, Ix>
    {
        for path in paths {
            // eprintln!("Mark path as visited: {path:?}");

            let new_path = match path.shift {
                PathShift::Right(_) => self.mark_path_as_visited_ins(graph, seq, tree, path),
                PathShift::Down(_) => self.mark_path_as_visited_del(graph, seq, tree, path)
            };

            if let Some(new) = new_path {
                if self.costs.gap_extend2() > 0 && self.costs.gap_extend2() != self.costs.gap_extend() {
                    queue.queue_visited_path(self.costs.gap_extend2() - 1, new.clone())
                }

                if self.costs.gap_extend() > 0 {
                    queue.queue_visited_path(self.costs.gap_extend() - 1, new);
                }

            }
        }
    }

    fn mark_path_as_visited_del<O, Ix, G, N, T>(&self, graph: &G, seq: &[u8], tree: &mut T, path: &VisitedPath<N, O>) -> Option<VisitedPath<N, O>>
    where
        O: OffsetType,
        Ix: TreeIndexType,
        G: AlignableGraph<NodeIndex=N>,
        N: NodeIndexType,
        T: AlignmentStateTree<N, O, Ix>,
    {
        let PathShift::Down(shift_down) = path.shift else {
            panic!("Called mark_path_as_visited_del() on a path with shift != PathShift::Down(O)");
        };

        let start: (N, O, O, RefCell<G::SuccessorIterator<'_>>) = (
            path.start_node, path.start_offset, O::one(),
            RefCell::new(graph.successors(path.start_node))
        );

        let mut stack = vec![start];
        let next_valid_child = |tree: &T, offset: O, dist_from_start: O, succ: &RefCell<G::SuccessorIterator<'_>>| -> Option<N> {
            // eprintln!("Max dist {dist_from_start:?} > {:?}", path.length.saturating_sub(&O::one()) + shift_down);
            // Exclude the path end-point, as it will be expanded on in another step
            if (path.start_offset + dist_from_start).as_usize() >= seq.len()
                    || dist_from_start >= path.length.saturating_sub(&O::one()) + shift_down {
                // eprintln!("{offset:?} >= {} || {dist_from_start:?} >= {:?}", seq.len(),
                //           path.length.saturating_sub(&O::one()) + shift_down);
                return None;
            }

            while let Some(child) = succ.borrow_mut().next() {
                let child_seq_offset = path.start_offset + dist_from_start;

                // eprint!("Checking {child:?}, offset: {:?}... ", child_seq_offset);
                // eprint!(" dist {dist_from_start:?} <= path_len: {:?}, symbol equal: {:?} ...", dist_from_start <= path.length,
                //         graph.is_symbol_equal(child, seq[child_seq_offset.as_usize() - 1]));
                if dist_from_start <= path.length
                    && graph.is_symbol_equal(child, seq[child_seq_offset.as_usize() - 1])
                {
                    // eprintln!("({child:?}, {child_seq_offset:?})");
                    return Some(child);
                }

                // If extending beyond the original path end point, make sure we don't go into already visited
                // states
                if dist_from_start > path.length
                    && !tree.visited(child, offset, AlignState::Match)
                {
                    // eprintln!("Going beyond. ({child:?}, {offset:?})");
                    return Some(child);
                }
                // eprintln!(" Nothing.")
            }

            None
        };

        let mut has_non_skipped_nodes = false;
        while !stack.is_empty() {
            let (parent_node, parent_offset, dist_from_start, succ) = stack.last().unwrap();
            // eprintln!("Parent: ({:?}, {:?}), dist: {:?}", parent_node, parent_offset, dist_from_start);

            if let Some(child) = next_valid_child(tree, *parent_offset, *dist_from_start, succ) {
                // eprintln!("Forward.");
                let child_offset;

                if *dist_from_start >= shift_down {
                    has_non_skipped_nodes = true;

                    tree.mark_visited(child, *parent_offset, AlignState::Deletion);
                    tree.mark_visited(child, *parent_offset, AlignState::Mismatch);

                    // eprintln!("Child: ({:?}, {:?}), dist: {:?}", child, *parent_offset, dist_from_start);

                    child_offset = parent_offset.increase_one();
                } else {
                    // eprintln!("Skip: ({:?}, {:?}), dist: {:?}", child, *parent_offset, dist_from_start);
                    child_offset = *parent_offset;
                }

                let child_succ = RefCell::new(graph.successors(child));
                let child_on_stack = (child, child_offset, dist_from_start.increase_one(), child_succ);
                stack.push(child_on_stack);
            } else {
                stack.pop();
            }
        }

        if has_non_skipped_nodes {
            Some(path.to_next())
        } else {
            None
        }
    }

    fn mark_path_as_visited_ins<O, Ix, G, N, T>(&self, graph: &G, seq: &[u8], tree: &mut T, path: &VisitedPath<N, O>) -> Option<VisitedPath<N, O>>
        where
            O: OffsetType,
            Ix: TreeIndexType,
            G: AlignableGraph<NodeIndex=N>,
            N: NodeIndexType,
            T: AlignmentStateTree<N, O, Ix>,
    {
        let PathShift::Right(shift_right) = path.shift else {
            panic!("Called mark_path_as_visited_ins() on a path with shift != PathShift::Right(O)");
        };

        let start = (
            path.start_node, path.start_offset, O::one(),
            RefCell::new(graph.successors(path.start_node))
        );
        tree.mark_visited(path.start_node, path.start_offset + shift_right, AlignState::Insertion);
        tree.mark_visited(path.start_node, path.start_offset + shift_right, AlignState::Mismatch);

        let mut stack = vec![start];

        let next_valid_child = |offset: O, dist_from_start: O, succ: &RefCell<G::SuccessorIterator<'_>>| -> Option<N> {
            if dist_from_start >= path.length.saturating_sub(&O::one())
                    || offset.as_usize() >= seq.len() {
                return None;
            }

            while let Some(child) = succ.borrow_mut().next() {
                let child_offset = offset.increase_one();
                if graph.is_symbol_equal(child, seq[child_offset.as_usize() - 1]) {
                    return Some(child);
                }
            }

            None
        };

        while !stack.is_empty() {
            let (_, parent_offset, dist_from_start, succ) = stack.last().unwrap();

            if let Some(child) = next_valid_child(*parent_offset, *dist_from_start, succ) {
                let child_offset = parent_offset.increase_one();

                tree.mark_visited(child, child_offset + shift_right, AlignState::Deletion);
                tree.mark_visited(child, child_offset + shift_right, AlignState::Mismatch);

                let child_succ = RefCell::new(graph.successors(child));
                let child_on_stack = (child, child_offset, dist_from_start.increase_one(), child_succ);
                stack.push(child_on_stack);
            } else {
                stack.pop();
            }
        }

        if (path.start_offset + shift_right).as_usize() < seq.len() {
            Some(path.to_next())
        } else {
            None
        }
    }

    fn extend<O, Ix, G, N, T>(
        &mut self,
        graph: &G,
        seq: &[u8],
        tree: &mut T,
        end_points: &mut [Ix],
        queue: &mut AlignStateQueue<N, O, Ix>,
    ) -> ExtendResult<Ix>
    where
        O: OffsetType,
        Ix: TreeIndexType,
        G: AlignableGraph<NodeIndex=N>,
        N: NodeIndexType,
        C: AlignmentCosts<StateTreeType<N, O, Ix> = T>,
        T: AlignmentStateTree<N, O, Ix>,
    {
        // Before extending, see if we already reached an end state
        if let Some(end) = end_points.iter()
            .find(|ix| {
                let node = tree.get_node(**ix);
                match node.state() {
                    AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                        node.offset().as_usize() == seq.len() && graph.is_end(node.node())
                    },
                    _ => false
                }
            })
        {
            return ReachedEnd(*end);
        }

        // Let's try to extend along matching sequence in the graph
        let mut additional_states = Vec::with_capacity(end_points.len());
        for node_ix in end_points.iter_mut() {
            let node = tree.get_node(*node_ix);

            match node.state() {
                AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                    let path_extender = PathExtender::new(graph, seq, tree, *node_ix);

                    let mut first = true;
                    for new_path in path_extender {
                        if first {
                            *node_ix = new_path.end();
                            first = false;
                        } else {
                            additional_states.push(new_path.end())

                        }

                        if self.costs.gap_extend() > 0 {
                            let visited_path_del = VisitedPath::from_extended_path(&new_path, AlignState::Deletion);
                            let visited_path_ins = VisitedPath::from_extended_path(&new_path, AlignState::Insertion);

                            let score_delta = self.costs.gap_open() + self.costs.gap_extend();

                            // We do score_delta -1 because we popped the current score from the queue in the main loop
                            queue.queue_visited_path(score_delta - 1, visited_path_del);
                            queue.queue_visited_path(score_delta - 1, visited_path_ins);
                        }

                        if self.costs.gap_extend2() > 0 {
                            let visited_path_del = VisitedPath::from_extended_path(&new_path, AlignState::Deletion2);
                            let visited_path_ins = VisitedPath::from_extended_path(&new_path, AlignState::Insertion2);

                            let score_delta = self.costs.gap_open2() + self.costs.gap_extend2();

                            // We do score_delta -1 because we popped the current score from the queue in the main loop
                            queue.queue_visited_path(score_delta - 1, visited_path_del);
                            queue.queue_visited_path(score_delta - 1, visited_path_ins);
                        }
                    }
                },
                AlignState::Deletion | AlignState::Deletion2 | AlignState::Insertion | AlignState::Insertion2 => ()
            }
        }

        // Check if one of our extended paths have reached the end
        if let Some(end) = end_points.iter().chain(additional_states.iter())
            .find(|ix| {
                let node = tree.get_node(**ix);
                match node.state() {
                    AlignState::Start | AlignState::Match | AlignState::Mismatch => {
                        node.offset().as_usize() == seq.len() && graph.is_end(node.node())
                    },
                    _ => false
                }
            })
        {
            return ReachedEnd(*end);
        }

        NewExtendedNodes(additional_states)
    }

    fn backtrace<O, Ix, N, T>(
        &self,
        tree: &T,
        end_state: Ix
    ) -> Alignment<N>
    where
        O: OffsetType,
        Ix: TreeIndexType,
        N: NodeIndexType,
        T: AlignmentStateTree<N, O, Ix>,
        C: AlignmentCosts<StateTreeType<N, O, Ix> = T>,
    {
        let mut curr = Some(end_state);
        let mut alignment = Alignment::new();

        while let Some(n) = curr {
            let state = tree.get_node(n);

            let Some(bt) = state.backtrace() else {
                break;
            };

            match state.state() {
                AlignState::Match | AlignState::Mismatch => {
                    match bt {
                        Backtrace::Step(_) => {
                            alignment.push(AlignedPair { rpos: Some(state.node()), qpos: Some(state.offset().as_usize() - 1) });
                        },
                        // On indel close we don't have to do anything, the next iteration will take care of the indel
                        Backtrace::ClosedIndel(_) => (),
                    }
                },
                AlignState::Insertion => {
                    alignment.push(AlignedPair { rpos: None, qpos: Some(state.offset().as_usize() - 1) });
                },
                AlignState::Deletion => {
                    alignment.push(AlignedPair { rpos: Some(state.node()), qpos: None });
                },
                AlignState::Start | AlignState::Insertion2 | AlignState::Deletion2 =>
                    panic!("Unexpected align state in backtrace!")
            }

            let prev_node = bt.prev();
            curr = match tree.get_node(prev_node).state() {
                AlignState::Start => None,
                _ => Some(prev_node)
            }
        }

        alignment.reverse();
        alignment
    }
}
