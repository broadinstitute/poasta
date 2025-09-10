use std::marker::PhantomData;

use crate::aligner::cost_models::AlignmentCostModel;
use crate::aligner::fr_points::Score;
use crate::aligner::traits::AlignableGraph;
use crate::aligner::utils::AlignedPair;
use super::AstarState;


pub fn create<C, G, F>(mut astar_state: C::AstarStateType<G>, initial_aln_states: &[C::Item], heuristic_func: F)
    -> AstarRunnable<C, G, F>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
    F: Fn(&C::AstarStateType<G>, &C::Item) -> usize,
{
    for aln_state in initial_aln_states {
        astar_state.update_if_further(aln_state, 0);
        astar_state.queue_item(aln_state.clone(), heuristic_func(&astar_state, aln_state));
    }

    AstarRunnable {
        state: astar_state,
        heuristic_func,
        dummy: PhantomData,
    }
}


pub struct AstarRunnable<C, G, F>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    state: C::AstarStateType<G>,
    heuristic_func: F,
    dummy: PhantomData<(C, G)>
}

impl<C, G, F> AstarRunnable<C, G, F>
where
    C: AlignmentCostModel,
    G: AlignableGraph,
    F: Fn(&C::AstarStateType<G>, &C::Item) -> usize,
    
{
    pub fn pop_front(&mut self) -> Option<C::Item> {
        self.state.pop_front()
    }
    
    pub fn is_end(&self, graph: &G, item: &C::Item) -> bool {
        self.state.is_end(graph, item)
    }

    pub fn get_score(&self, item: &C::Item) -> Score {
        self.state.get_score(item)
    }
    
    pub fn relax(&mut self, graph: &G, seq: &[u8], item: &C::Item) {
        self.state.relax(graph, seq, item, |state, e| {
            (self.heuristic_func)(state, e)
        });
    }
    
    pub fn backtrace(&self, graph: &G, end: &C::Item) -> Vec<AlignedPair<G::NodePosType>> {
        self.state.backtrace(graph, end)
    }
    
}

