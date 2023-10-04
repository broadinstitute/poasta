//! A module containing a mock graph struct useful for creating
//! test graphs in unit tests

use rustc_hash::FxHashMap;
use petgraph::graph::{DiGraph, NodeIndex, NodeIndices, Neighbors};
use petgraph::Incoming;

use crate::graphs::AlignableGraph;

pub(crate) type NIx = u32;


pub(crate) type MockGraph = DiGraph<i64, (), NIx>;

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

pub(crate) fn create_test_graph1() -> MockGraph {
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

pub(crate) fn create_test_graph2() -> MockGraph {
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
