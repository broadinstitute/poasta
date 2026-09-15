//! Superbubble detection.

use std::cmp::min;
use std::iter::Rev;
use std::ops::Range;

use crate::graph::traits::GraphWithNodeOrdering;

/// Identify superbubbles in a directed acyclic graph.
///
/// This algorithm is described in Gärtner et al. \[1\].
///
/// 1. Gärtner, Fabian, Lydia Müller, and Peter F. Stadler. "Superbubbles Revisited."
///    Algorithms for Molecular Biology 13, no. 1 (December 1, 2018): 16.
///    <https://doi.org/10.1186/s13015-018-0134-3>.
///
/// The algorithm works on a DFS-based topological order, which is what
/// `node_rank` already is, so every per-node table is a `Vec` indexed by rank.
pub struct SuperbubbleFinder<'a, G>
where
    G: GraphWithNodeOrdering,
{
    graph: &'a G,
    /// By rank: the lowest rank among the node's predecessors, `-1` if none.
    out_parent: Vec<i64>,
    /// By rank: the highest rank among the node's successors, `i64::MAX` if none.
    out_child: Vec<i64>,
}

impl<'a, G> SuperbubbleFinder<'a, G>
where
    G: GraphWithNodeOrdering,
{
    pub fn new(graph: &'a G) -> Self {
        let n = graph.node_count();
        let mut out_parent = vec![-1i64; n];
        let mut out_child = vec![i64::MAX; n];
        for rank in 0..n {
            let node = graph.rank_to_node(rank);
            if let Some(min_pred) = graph.predecessors(node).map(|p| graph.node_rank(p) as i64).min() {
                out_parent[rank] = min_pred;
            }
            if let Some(max_succ) = graph.successors(node).map(|s| graph.node_rank(s) as i64).max() {
                out_child[rank] = max_succ;
            }
        }

        Self {
            graph,
            out_parent,
            out_child,
        }
    }

    pub fn iter(&self) -> SuperbubbleIterator<'_, G> {
        SuperbubbleIterator::new(self)
    }
}

pub struct SuperbubbleIterator<'a, G>
where
    G: GraphWithNodeOrdering,
{
    finder: &'a SuperbubbleFinder<'a, G>,
    // Indexed by rank. Every entry is written when its node is processed and
    // only read for nodes processed earlier, so the initial value is never
    // seen. `i64::MAX` is still the harmless choice since the only update is
    // a `min`.
    out_parent_map: Vec<i64>,
    stack: Vec<G::NodeType>,
    curr: Rev<Range<usize>>,
    candidate_exit: Option<G::NodeType>,
}

impl<'a, G> SuperbubbleIterator<'a, G>
where
    G: GraphWithNodeOrdering,
{
    fn new(finder: &'a SuperbubbleFinder<'a, G>) -> Self {
        Self {
            finder,
            out_parent_map: vec![i64::MAX; finder.graph.node_count()],
            stack: Vec::default(),
            curr: (0..finder.graph.node_count()).rev(),
            candidate_exit: None,
        }
    }
}

impl<'a, G> Iterator for SuperbubbleIterator<'a, G>
where
    G: GraphWithNodeOrdering,
{
    type Item = (G::NodeType, G::NodeType);

    fn next(&mut self) -> Option<Self::Item> {
        let finder = self.finder;
        let graph = finder.graph;

        while let Some(curr) = self.curr.next() {
            let mut to_return = None;

            let n = graph.rank_to_node(curr);
            let furthest_child = finder.out_child[curr];

            if furthest_child == curr as i64 + 1 {
                if let Some(c) = self.candidate_exit {
                    self.stack.push(c);
                }

                self.candidate_exit = Some(graph.rank_to_node(curr + 1))
            } else {
                while let Some(candidate) = self.candidate_exit {
                    if furthest_child <= graph.node_rank(candidate) as i64 {
                        break;
                    }

                    self.candidate_exit = self.stack.pop();

                    if let Some(c) = self.candidate_exit {
                        let new_out_parent = min(
                            self.out_parent_map[graph.node_rank(candidate)],
                            self.out_parent_map[graph.node_rank(c)],
                        );
                        self.out_parent_map[graph.node_rank(c)] = new_out_parent;
                    }
                }
            }

            if let Some(candidate) = self.candidate_exit {
                if self.out_parent_map[graph.node_rank(candidate)] as usize == curr {
                    to_return = Some((n, candidate));

                    self.candidate_exit = self.stack.pop();
                    if let Some(c) = self.candidate_exit {
                        let new_out_parent = min(
                            self.out_parent_map[graph.node_rank(candidate)],
                            self.out_parent_map[graph.node_rank(c)],
                        );
                        self.out_parent_map[graph.node_rank(c)] = new_out_parent;
                    }
                }
            }

            self.out_parent_map[curr] = finder.out_parent[curr];

            if let Some(c) = self.candidate_exit {
                let new_out_parent = min(self.out_parent_map[curr], self.out_parent_map[graph.node_rank(c)]);
                self.out_parent_map[graph.node_rank(c)] = new_out_parent;
            }

            if to_return.is_some() {
                return to_return;
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::SuperbubbleFinder;
    use crate::graph::mock::{create_test_graph1, create_test_graph2};
    use crate::graph::traits::GraphWithNodeOrdering;
    use rustc_hash::FxHashSet;

    #[test]
    fn test_superbubble_finder() {
        let g1 = create_test_graph1();
        let finder = SuperbubbleFinder::new(&g1);
        let bubbles1: FxHashSet<_> = finder
            .iter()
            .filter(|(s, t)| *s != g1.end_node() && *t != g1.end_node())
            .map(|(n1, n2)| (g1[n1].0, g1[n2].0))
            .collect();

        assert_eq!(
            bubbles1,
            FxHashSet::from_iter([(1, 2), (2, 3), (4, 5), (5, 6), (7, 8), (8, 9),])
        );

        let g2 = create_test_graph2();
        let finder = SuperbubbleFinder::new(&g2);
        let bubbles2: FxHashSet<_> = finder
            .iter()
            .map(|(n1, n2)| (g2[n1].0, g2[n2].0))
            .collect();

        assert_eq!(
            bubbles2,
            FxHashSet::from_iter([(8, 15), (11, 12), (5, 7), (3, 8), (1, 3)])
        );
    }
}
