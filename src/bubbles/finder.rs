use std::cmp::min;
use std::iter::Rev;
use std::ops::Range;
use rustc_hash::FxHashMap;
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::graphs::tools::rev_postorder_nodes;

/// Identify superbubbles in a directed acyclic graph
///
/// This algorithm is described in Gaertner et al. \[1\].
///
/// 1. Gärtner, Fabian, Lydia Müller, and Peter F. Stadler. “Superbubbles Revisited.”
///    Algorithms for Molecular Biology 13, no. 1 (December 1, 2018): 16.
///    https://doi.org/10.1186/s13015-018-0134-3.
pub struct SuperbubbleFinder<'a, G>
where
    G: AlignableGraph
{
    graph: &'a G,
    rev_postorder: Vec<usize>,
    inv_rev_postorder: Vec<G::NodeIndex>,
    out_parent: FxHashMap<G::NodeIndex, i64>,
    out_child: FxHashMap<G::NodeIndex, i64>,
}

impl<'a, G> SuperbubbleFinder<'a, G>
where
    G: AlignableGraph,
{
    pub fn new(graph: &'a G) -> Self {
        let inv_rev_postorder = rev_postorder_nodes(graph);
        let mut rev_postorder = vec![0; inv_rev_postorder.len()];
        for (postorder, node_ix) in inv_rev_postorder.iter().enumerate() {
            rev_postorder[node_ix.index()] = postorder
        }

        let out_parent = graph.all_nodes()
            .map(|n| {
                let min_pred = graph.predecessors(n)
                    .map(|p| rev_postorder[p.index()] as i64)
                    .min()
                    .unwrap_or(-1);

                (n, min_pred)
            })
            .collect();

        let out_child = graph.all_nodes()
            .map(|n| {
                let max_succ = graph.successors(n)
                    .map(|s| rev_postorder[s.index()] as i64)
                    .max()
                    .unwrap_or(i64::MAX);

                (n, max_succ)
            })
            .collect();

        Self {
            graph,
            rev_postorder,
            inv_rev_postorder,
            out_parent,
            out_child,
        }
    }

    #[inline(always)]
    pub(crate) fn get_inv_rev_postorder(&self) -> &[G::NodeIndex] {
        &self.inv_rev_postorder
    }

    pub fn iter(&self) -> SuperbubbleIterator<'_, G> {
        SuperbubbleIterator::new(self)
    }
}


pub struct SuperbubbleIterator<'a, G>
where
    G: AlignableGraph,
{
    finder: &'a SuperbubbleFinder<'a, G>,
    out_parent_map: FxHashMap<G::NodeIndex, i64>,
    stack: Vec<G::NodeIndex>,
    curr: Rev<Range<usize>>,
    candidate_exit: Option<G::NodeIndex>,
}

impl<'a, G> SuperbubbleIterator<'a, G>
where
    G: AlignableGraph,
{
    fn new(finder: &'a SuperbubbleFinder<'a, G>) -> Self {
        Self {
            finder,
            out_parent_map: FxHashMap::default(),
            stack: Vec::default(),
            curr: (0..finder.graph.node_count_with_start_and_end()).rev(),
            candidate_exit: None,
        }
    }
}


impl<'a, G> Iterator for SuperbubbleIterator<'a, G>
where
    G: AlignableGraph,
{
    type Item = (G::NodeIndex, G::NodeIndex);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(curr) = self.curr.next() {
            let mut to_return = None;

            let n = self.finder.inv_rev_postorder[curr];
            let furthest_child = self.finder.out_child[&n];

            if furthest_child == curr as i64 + 1 {
                if let Some(c) = self.candidate_exit {
                    self.stack.push(c);
                }

                self.candidate_exit = Some(self.finder.inv_rev_postorder[curr + 1])
            } else {
                while let Some(candidate) = self.candidate_exit {
                    if furthest_child <= self.finder.rev_postorder[candidate.index()] as i64 {
                        break;
                    }

                    self.candidate_exit = self.stack.pop();

                    if let Some(c) = self.candidate_exit {
                        let new_out_parent = min(self.out_parent_map[&candidate], self.out_parent_map[&c]);
                        self.out_parent_map.entry(c)
                            .and_modify(|v| *v = new_out_parent)
                            .or_insert(new_out_parent);
                    }
                }
            }

            if let Some(candidate) = self.candidate_exit {
                if self.out_parent_map[&candidate] as usize == curr {
                    to_return = Some((n, candidate));

                    self.candidate_exit = self.stack.pop();
                    if let Some(c) = self.candidate_exit {
                        let new_out_parent = min(self.out_parent_map[&candidate], self.out_parent_map[&c]);
                        self.out_parent_map.entry(c)
                            .and_modify(|v| *v = new_out_parent)
                            .or_insert(new_out_parent);
                    }
                }
            }

            self.out_parent_map.entry(n)
                .and_modify(|v| *v = self.finder.out_parent[&n])
                .or_insert(self.finder.out_parent[&n]);

            if let Some(c) = self.candidate_exit {
                let new_out_parent = min(self.out_parent_map[&n], self.out_parent_map[&c]);
                self.out_parent_map.entry(c)
                    .and_modify(|v| *v = new_out_parent)
                    .or_insert(new_out_parent);

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
    use rustc_hash::FxHashSet;
    use crate::bubbles::finder::SuperbubbleFinder;
    use crate::graphs::mock::{create_test_graph1, create_test_graph2};


    #[test]
    fn test_superbubble_finder() {
        let g1 = create_test_graph1();
        let finder = SuperbubbleFinder::new(&g1);
        let bubbles1: FxHashSet<_> = finder.iter()
            .map(|(n1, n2)| (g1[n1], g1[n2]))
            .collect();

        assert_eq!(bubbles1, FxHashSet::from_iter([
            (1, 2),
            (2, 3),
            (4, 5),
            (5, 6),
            (7, 8),
            (8, 9)
        ]));

        let g2 = create_test_graph2();
        let finder = SuperbubbleFinder::new(&g2);
        let bubbles2: FxHashSet<_> = finder.iter()
            .map(|(n1, n2)| (g2[n1], g2[n2]))
            .collect();

        assert_eq!(bubbles2, FxHashSet::from_iter([
            (8, 14),
            (11, 12),
            (5, 7),
            (3, 8),
            (1, 3)
        ]));
    }
}