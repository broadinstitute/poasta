use std::collections::{HashSet};
use rustc_hash::FxHashSet;
use crate::aligner::offsets::OffsetType;
use crate::graphs::{AlignableGraph, NodeIndexType};

pub trait VisitedSet<N, O> {
    fn visited(&self, node: N, offset: O) -> bool;
    fn mark_visited(&mut self, node: N, offset: O);
}

impl<N, O> VisitedSet<N, O> for HashSet<(N, O)>
    where
        N: NodeIndexType,
        O: OffsetType
{
    fn visited(&self, node: N, offset: O) -> bool {
        let key = (node, offset);
        self.contains(&key)
    }

    fn mark_visited(&mut self, node: N, offset: O) {
        self.insert((node, offset));
    }
}

pub struct VisitedSetPerNode<O> {
    visited: Vec<FxHashSet<O>>
}

impl<O> VisitedSetPerNode<O>
where
    O: OffsetType,
{
    pub fn new<G>(graph: &G) -> Self
    where
        G: AlignableGraph,
    {
        Self {
            visited: vec![FxHashSet::default(); graph.node_count()]
        }
    }
}

impl<N, O> VisitedSet<N, O> for VisitedSetPerNode<O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    fn visited(&self, node: N, offset: O) -> bool {
        self.visited[node.index()].contains(&offset)
    }

    fn mark_visited(&mut self, node: N, offset: O) {
        self.visited[node.index()].insert(offset);
    }
}

mod ahash_set {
    use ahash::AHashSet;
    use crate::aligner::offsets::OffsetType;
    use crate::graphs::NodeIndexType;

    impl<N, O> super::VisitedSet<N, O> for AHashSet<(N, O)>
        where
            N: NodeIndexType,
            O: OffsetType
    {
        fn visited(&self, node: N, offset: O) -> bool {
            let key = (node, offset);
            self.contains(&key)
        }

        fn mark_visited(&mut self, node: N, offset: O) {
            self.insert((node, offset));
        }
    }
}





