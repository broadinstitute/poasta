use std::cmp::min;
use std::iter::Rev;
use std::ops::Range;
use rustc_hash::FxHashMap;
use crate::graphs::AlignableGraph;
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
    rev_postorder: FxHashMap<G::NodeIndex, usize>,
    inv_rev_postorder: Vec<G::NodeIndex>,
    out_parent: FxHashMap<G::NodeIndex, i64>,
    out_child: FxHashMap<G::NodeIndex, i64>,
    out_parent_map: FxHashMap<G::NodeIndex, i64>,
    stack: Vec<G::NodeIndex>,
    curr: Rev<Range<usize>>,
    candidate_exit: Option<G::NodeIndex>,
}

impl<'a, G> SuperbubbleFinder<'a, G>
where
    G: AlignableGraph,
{
    pub fn new(graph: &'a G) -> Self {
        let inv_rev_postorder = rev_postorder_nodes(graph);
        let rev_postorder: FxHashMap<G::NodeIndex, usize> = inv_rev_postorder.iter()
            .enumerate()
            .map(|(order, node)| (*node, order))
            .collect();

        let out_parent = graph.all_nodes()
            .map(|n| {
                let min_pred = graph.predecessors(n)
                    .map(|p| rev_postorder[&p] as i64)
                    .min()
                    .unwrap_or(-1);

                (n, min_pred)
            })
            .collect();

        let out_child = graph.all_nodes()
            .map(|n| {
                let max_succ = graph.successors(n)
                    .map(|p| rev_postorder[&p] as i64)
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
            out_parent_map: FxHashMap::default(),
            stack: Vec::default(),
            curr: (0..graph.node_count_with_start()).rev(),
            candidate_exit: None,
        }
    }
}

impl<'a, G> Iterator for SuperbubbleFinder<'a, G>
where
    G: AlignableGraph,
{
    type Item = (G::NodeIndex, G::NodeIndex);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(curr) = self.curr.next() {
            let mut to_return = None;

            let n = self.inv_rev_postorder[curr];
            let furthest_child = self.out_child[&n];

            if furthest_child == curr as i64 + 1 {
                if let Some(c) = self.candidate_exit {
                    self.stack.push(c);
                }

                self.candidate_exit = Some(self.inv_rev_postorder[curr + 1])
            } else {
                while let Some(candidate) = self.candidate_exit {
                    if furthest_child <= self.rev_postorder[&candidate] as i64 {
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
                .and_modify(|v| *v = self.out_parent[&n])
                .or_insert(self.out_parent[&n]);

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
    use ahash::HashSet;
    use petgraph::graph::{DiGraph, NodeIndex, NodeIndices, Neighbors};
    use petgraph::Incoming;
    use rustc_hash::FxHashMap;
    use crate::bubbles::finder::SuperbubbleFinder;
    use crate::graphs::AlignableGraph;

    type NIx = usize;

    type MockGraph = DiGraph<i64, (), NIx>;

    impl AlignableGraph for MockGraph {
        type NodeIndex = NodeIndex<NIx>;

        type NodeIterator<'a> = NodeIndices<NIx>
            where Self: 'a;

        type PredecessorIterator<'a> = Neighbors<'a, (), NIx>
            where Self: 'a;
        type SuccessorIterator<'a> = Neighbors<'a, (), NIx>
            where Self: 'a;

        fn all_nodes(&self) -> Self::NodeIterator<'_> {
            self.node_indices()
        }

        fn node_count(&self) -> usize {
            self.node_count()
        }

        fn node_count_with_start(&self) -> usize {
            self.node_count()
        }

        fn edge_count(&self) -> usize {
            todo!()
        }

        fn start_node(&self) -> Self::NodeIndex {
            self.node_indices().next().unwrap()
        }

        fn predecessors(&self, node: Self::NodeIndex) -> Self::PredecessorIterator<'_> {
            self.neighbors_directed(node, Incoming)
        }

        fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_> {
            self.neighbors(node)
        }

        fn is_end(&self, _: Self::NodeIndex) -> bool {
            false
        }

        fn get_symbol(&self, _: Self::NodeIndex) -> char {
            '-'
        }

        fn is_symbol_equal(&self, _: Self::NodeIndex, _: u8) -> bool {
            false
        }
    }

    fn create_test_graph1() -> MockGraph {
        let mut nmap = FxHashMap::default();
        let mut g = MockGraph::default();

        for i in 1..=9 {
            let nix = g.add_node(i);
            nmap.insert(i, nix);
        }

        let edges = [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (3, 7),
            (7, 8),
            (8, 9)
        ];

        for (s, t) in edges.iter() {
            g.add_edge(nmap[s], nmap[t], ());
        }

        g
    }

    fn create_test_graph2() -> MockGraph {
        let mut nmap = FxHashMap::default();
        let mut g = MockGraph::default();

        for i in 1..=15 {
            let nix = g.add_node(i);
            nmap.insert(i, nix);
        }

        let edges = [
            (1, 2),
            (1, 3),
            (2, 3),
            (3, 4),
            (3, 5),
            (3, 11),
            (4, 8),
            (5, 6),
            (5, 9),
            (6, 7),
            (6, 10),
            (7, 8),
            (8, 13),
            (8, 14),
            (9, 10),
            (10, 7),
            (11, 12),
            (12, 8),
            (13, 14),
            (13, 15),
            (15, 14)
        ];

        for (s, t) in edges.iter() {
            g.add_edge(nmap[s], nmap[t], ());
        }

        g
    }

    #[test]
    fn test_superbubble_finder() {
        let g1 = create_test_graph1();
        let finder = SuperbubbleFinder::new(&g1);
        let bubbles1: HashSet<_> = finder
            .map(|(n1, n2)| (g1[n1], g1[n2]))
            .collect();

        assert_eq!(bubbles1, HashSet::from_iter([
            (1, 2),
            (2, 3),
            (4, 5),
            (5, 6),
            (7, 8),
            (8, 9)
        ]));

        let g2 = create_test_graph2();
        let finder = SuperbubbleFinder::new(&g2);
        let bubbles2: HashSet<_> = finder
            .map(|(n1, n2)| (g2[n1], g2[n2]))
            .collect();

        assert_eq!(bubbles2, HashSet::from_iter([
            (8, 14),
            (11, 12),
            (5, 7),
            (3, 8),
            (1, 3)
        ]));
    }
}