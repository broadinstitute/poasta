//! A* heuristics: admissible lower bounds on the cost still to pay.
//!
//! [`AstarHeuristic`] is what the search calls. [`HeuristicConfig`] is the
//! small value an engine carries to build one per alignment. This module
//! holds both traits and Dijkstra (`h = 0`). `min_gap.rs` holds the
//! minimum-gap-cost heuristic, the default, and `chain/` the gap-chaining
//! seed heuristic.

pub mod chain;
mod min_gap;

#[cfg(feature = "viz")]
pub use chain::{ChainExplain, MatchId, MatchInfo, Point, Potential, Staircase, SubPoint};
pub use chain::{SeedParams, SmoothHeuristic, SmoothHeuristicConfig};
pub use min_gap::{MinGapCostHeuristic, MinimumGapCost, MinimumGapCostAffine};

use std::marker::PhantomData;
use std::sync::Arc;

use crate::align::AlignState;
use crate::align::AlignmentMode;
use crate::align::astar::aln_graph::AlignmentGraphNode;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::offsets::OffsetType;
use crate::align::traits::AlignableGraph;
use crate::graph::bubbles::index::BubbleIndex;
use crate::graph::traits::GraphNodeId;

pub trait AstarHeuristic<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    /// Get ready for one alignment. Called once per `align`, before the
    /// first `h`. A heuristic keeps its buffers between calls and resets only
    /// what the previous alignment used.
    fn prepare<G>(
        &mut self,
        graph: &G,
        query: &[u8],
        bubbles: &Arc<BubbleIndex<N>>,
        mode: AlignmentMode,
    ) where
        G: AlignableGraph<Node = N>;

    /// Takes `&mut self`: a heuristic may update its own state on a query.
    /// `graph` is the one `prepare` saw.
    fn h<G>(
        &mut self,
        graph: &G,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> usize
    where
        G: AlignableGraph<Node = N>;

    /// Called for every expanded `Match` state: popped ones with
    /// `h_limit = usize::MAX`, and the ones greedy extension passes with the
    /// popped state's h. A heuristic that deletes seed matches acts only if
    /// `h` here is at most `h_limit`, so that the state is at its optimal
    /// score. Returns whether it changed anything.
    fn on_expand<G>(
        &mut self,
        _graph: &G,
        _aln_node: &AlignmentGraphNode<N, O>,
        _h_limit: usize,
    ) -> bool
    where
        G: AlignableGraph<Node = N>,
    {
        false
    }

    /// Whether `h` can grow during a search. The engine then re-scores
    /// popped states and pushes back those whose priority went up.
    fn is_dynamic(&self) -> bool {
        false
    }
}

/// Heuristic that always returns 0, reducing A* to plain Dijkstra: exact and
/// admissible, but explores every state within the optimal-score radius.
pub struct Dijkstra<N, O>(PhantomData<(N, O)>);

// Manual impl: `#[derive(Default)]` would incorrectly require `N: Default, O:
// Default` even though `PhantomData` doesn't need it.
impl<N, O> Default for Dijkstra<N, O> {
    fn default() -> Self {
        Self(PhantomData)
    }
}

impl<N, O> AstarHeuristic<N, O> for Dijkstra<N, O>
where
    N: GraphNodeId,
    O: OffsetType,
{
    fn prepare<G>(
        &mut self,
        _graph: &G,
        _query: &[u8],
        _bubbles: &Arc<BubbleIndex<N>>,
        _mode: AlignmentMode,
    ) where
        G: AlignableGraph<Node = N>,
    {
    }

    fn h<G>(
        &mut self,
        _graph: &G,
        _aln_node: &AlignmentGraphNode<N, O>,
        _aln_state: AlignState,
    ) -> usize
    where
        G: AlignableGraph<Node = N>,
    {
        0
    }
}

/// Picks a heuristic. The config is the small, `Copy` value that callers
/// and the CLI pass around, and the engine's `H` type parameter.
/// `instantiate` makes the stateful heuristic that the engine keeps across
/// alignments.
pub trait HeuristicConfig {
    type Heuristic<N: GraphNodeId, O: OffsetType, C: AlignmentCostModel + Copy>: AstarHeuristic<N, O>;

    fn instantiate<N, O, C>(&self, costs: C) -> Self::Heuristic<N, O, C>
    where
        N: GraphNodeId,
        O: OffsetType,
        C: AlignmentCostModel + Copy;
}

/// `h ≡ 0`: A* reduces to Dijkstra. Exact and admissible, but explores every
/// state within the optimal-score radius.
#[derive(Default, Clone, Copy)]
pub struct DijkstraHeuristic;

impl HeuristicConfig for DijkstraHeuristic {
    type Heuristic<N: GraphNodeId, O: OffsetType, C: AlignmentCostModel + Copy> = Dijkstra<N, O>;

    fn instantiate<N, O, C>(&self, _costs: C) -> Dijkstra<N, O>
    where
        N: GraphNodeId,
        O: OffsetType,
        C: AlignmentCostModel + Copy,
    {
        Dijkstra::default()
    }
}
