//! A module containing a mock graph struct useful for creating
//! test graphs in unit tests

use std::ops::{Deref, DerefMut};

use petgraph::graph::{DiGraph, Neighbors, NodeIndex, NodeIndices};
use petgraph::{Incoming, Outgoing};
use rustc_hash::FxHashMap;

use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

use super::traits::GraphWithNodeLengths;
use super::utils::rev_postorder_nodes;

pub(crate) type NIx = u32;

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub(crate) struct NodeData(pub i64, pub usize);

impl NodeData {
    fn new(label: i64) -> Self {
        NodeData(label, 0)
    }
}

#[derive(Default)]
pub(crate) struct MockGraph(DiGraph<NodeData, (), NIx>, Vec<NodeIndex<NIx>>);

impl MockGraph {
    fn post_process_graph(&mut self) {
        let rev_postorder = rev_postorder_nodes(self);

        for (rank, node) in rev_postorder.iter().enumerate() {
            self.0[*node].1 = rank;
        }
    }
}

impl Deref for MockGraph {
    type Target = DiGraph<NodeData, (), NIx>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for MockGraph {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl GraphBase for MockGraph {
    type NodeType = NodeIndex<NIx>;

    type NodeIter<'a>
        = NodeIndices<NIx>
    where
        Self: 'a;

    type Successors<'a>
        = Neighbors<'a, (), NIx>
    where
        Self: 'a;
    type Predecessors<'a>
        = Neighbors<'a, (), NIx>
    where
        Self: 'a;

    fn all_nodes_iter(&self) -> Self::NodeIter<'_> {
        self.node_indices()
    }

    fn node_count(&self) -> usize {
        self.0.node_count()
    }

    fn is_empty(&self) -> bool {
        self.0.node_count() == 0
    }

    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.neighbors_directed(node, Incoming)
    }

    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.neighbors(node)
    }

    fn in_degree(&self, node: Self::NodeType) -> usize {
        self.neighbors_directed(node, Incoming).count()
    }

    fn out_degree(&self, node: Self::NodeType) -> usize {
        self.neighbors(node).count()
    }
}

impl GraphWithNodeOrdering for MockGraph {
    fn start_node(&self) -> Self::NodeType {
        self.node_indices()
            .find(|n| self.neighbors_directed(*n, Incoming).count() == 0)
            .unwrap()
    }

    fn end_node(&self) -> Self::NodeType {
        self.node_indices()
            .find(|n| self.neighbors_directed(*n, Outgoing).count() == 0)
            .unwrap()
    }

    fn node_rank(&self, node: Self::NodeType) -> usize {
        self.0[node].1
    }

    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType {
        self.1[node_rank]
    }
}

impl GraphWithNodeLengths for MockGraph {
    fn node_length(&self, _: Self::NodeType) -> usize {
        1
    }
}

pub(crate) fn create_test_graph1() -> MockGraph {
    let mut nmap = FxHashMap::default();
    let mut g = MockGraph::default();

    for i in 1..=9 {
        let nix = g.add_node(NodeData::new(i));
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
        (8, 9),
    ];

    for (s, t) in edges.iter() {
        g.add_edge(nmap[s], nmap[t], ());
    }

    // Create a mock "end node"
    let end_node = g.add_node(NodeData::new(10));
    for n in g.node_indices() {
        if n != end_node && g.successors(n).count() == 0 {
            g.add_edge(n, end_node, ());
        }
    }

    g.post_process_graph();

    g
}

pub(crate) fn create_test_graph2() -> MockGraph {
    let mut nmap = FxHashMap::default();
    let mut g = MockGraph::default();

    for i in 1..=15 {
        let nix = g.add_node(NodeData::new(i));
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
        (8, 15),
        (9, 10),
        (10, 7),
        (11, 12),
        (12, 8),
        (13, 14),
        (13, 15),
        (14, 15),
    ];

    for (s, t) in edges.iter() {
        g.add_edge(nmap[s], nmap[t], ());
    }

    g.post_process_graph();

    g
}
