pub use petgraph::graph::IndexType;
use petgraph::stable_graph::{Neighbors, NodeIndices};
use petgraph::visit::EdgeRef;
use petgraph::{Incoming, Outgoing};

use crate::align::engine::{AlignResult, AlignedPair};
use crate::errors::PoastaError;
use crate::graph::alignment::AddAlignment;

use super::traits::{GraphBase, GraphWithNodeOrdering};
use super::utils::rev_postorder_nodes;

pub(crate) mod graph_impl {
    use petgraph::graph::IndexType;
    use petgraph::{stable_graph::StableDiGraph, visit::GraphBase};

    pub type POAGraphImpl<Ix> = StableDiGraph<POANodeData<Ix>, POAEdgeData, Ix>;
    pub type POANodeIndex<Ix> = <POAGraphImpl<Ix> as GraphBase>::NodeId;

    #[derive(Debug)]
    pub struct POANodeData<Ix>
    where
        Ix: IndexType,
    {
        pub symbol: u8,
        pub aligned_nodes: Vec<POANodeIndex<Ix>>,
        pub rank: Ix,
        /// Shortest path length (in edges) from the start sentinel.
        pub depth_min: Ix,
        /// Longest path length (in edges) from the start sentinel.
        pub depth_max: Ix,
    }

    impl<Ix> POANodeData<Ix>
    where
        Ix: IndexType,
    {
        pub fn new(symbol: u8) -> Self {
            POANodeData {
                symbol,
                aligned_nodes: Vec::new(),
                rank: Ix::new(0),
                depth_min: Ix::new(0),
                depth_max: Ix::new(0),
            }
        }
    }

    #[derive(Debug, Default)]
    pub struct POAEdgeData {
        pub weight: usize,
        pub sequence_ids: Vec<usize>,
    }

    impl POAEdgeData {
        pub fn new_with_seq_id(seq_id: usize, weight: usize) -> Self {
            POAEdgeData {
                weight,
                sequence_ids: vec![seq_id],
            }
        }

        pub fn new_for_start_or_end() -> Self {
            POAEdgeData {
                weight: 0,
                sequence_ids: Vec::default(),
            }
        }
    }
}

use crate::align::traits::AlignableGraph;
use crate::graph::traits::GraphNodeId;
pub use graph_impl::POANodeIndex;
use graph_impl::{POAEdgeData, POANodeData};

/// A sequence aligned to the POA graph.
///
/// Stores the sequence name and the start node in the graph.
#[derive(Debug, Clone)]
pub struct Sequence<Ix>(pub(crate) String, pub(crate) Ix)
where
    Ix: IndexType;

impl<Ix> Sequence<Ix>
where
    Ix: IndexType,
{
    pub fn name(&self) -> &String {
        &self.0
    }

    pub fn start_node(&self) -> Ix {
        self.1
    }
}

/// A partial order alignment graph
///
/// A POA graph represents a multiple sequence alignment as a directed acyclic graph. Character-labeled
/// nodes represent the symbols in the input sequences and edges indicate which symbols are adjacent in
/// at least one sequence.
///
/// Each sequence used as input for graph construction can be reconstructed by tracing the path of nodes,
/// and the sequence start node is stored in `sequences`.
#[derive(Debug, Default)]
pub struct POAGraph<Ix = u32>
where
    Ix: IndexType,
{
    pub(crate) graph: graph_impl::POAGraphImpl<Ix>,
    pub sequences: Vec<Sequence<POANodeIndex<Ix>>>,
    topological_sorted: Vec<POANodeIndex<Ix>>,
    start_node: POANodeIndex<Ix>,
    end_node: POANodeIndex<Ix>,
}

impl<Ix> POAGraph<Ix>
where
    Ix: IndexType,
{
    pub fn new() -> Self {
        let mut graph = graph_impl::POAGraphImpl::<Ix>::default();
        let start_node = graph.add_node(POANodeData::new(b'#'));
        let end_node = graph.add_node(POANodeData::new(b'$'));
        let _initial_edge =
            graph.add_edge(start_node, end_node, POAEdgeData::new_for_start_or_end());

        Self {
            graph,
            sequences: Vec::default(),
            topological_sorted: vec![start_node, end_node],
            start_node,
            end_node,
        }
    }

    pub fn node_symbol(&self, node: POANodeIndex<Ix>) -> u8 {
        self.graph[node].symbol
    }

    pub(crate) fn add_edge(
        &mut self,
        s: POANodeIndex<Ix>,
        t: POANodeIndex<Ix>,
        sequence_id: usize,
        weight: usize,
    ) {
        // If edge exists, update sequence ID and weight of the existing one
        if let Some(e) = self.graph.find_edge(s, t) {
            let edge_data = self.graph.edge_weight_mut(e).unwrap();
            edge_data.sequence_ids.push(sequence_id);
            edge_data.weight += weight;
        } else {
            self.graph
                .add_edge(s, t, POAEdgeData::new_with_seq_id(sequence_id, weight));
        }
    }

    pub fn add_nodes_for_sequence(
        &mut self,
        sequence: &[u8],
        weights: &[usize],
        start: usize,
        end: usize,
    ) -> Option<(POANodeIndex<Ix>, POANodeIndex<Ix>)> {
        if start == end {
            return None;
        }

        let mut first_node = None;
        let mut prev = None;
        for pos in start..end {
            let curr_node = self.graph.add_node(POANodeData::new(sequence[pos]));

            if first_node.is_none() {
                first_node = Some(curr_node);
            }

            if let Some(prev_node) = prev {
                self.add_edge(
                    prev_node,
                    curr_node,
                    self.sequences.len(),
                    weights[pos - 1] + weights[pos],
                )
            }

            prev = Some(curr_node)
        }

        Some((first_node.unwrap(), prev.unwrap()))
    }

    pub(crate) fn post_process(&mut self) -> Result<(), PoastaError<Ix>> {
        self.topological_sorted.clear();

        // By repeatedly immediately calling next() on the Edges iterator returned by edges(), we
        // ensure that the returned EdgeIndex is always valid. If using normal iteration, the removal
        // of an edge might invalidate following edge indices.
        while let Some(e) = self.graph.edges(self.start_node).next() {
            self.graph.remove_edge(e.id());
        }

        while let Some(e) = self.graph.edges_directed(self.end_node, Incoming).next() {
            self.graph.remove_edge(e.id());
        }

        // Connect nodes with no incoming edges to the start node
        let all_nodes: Vec<POANodeIndex<Ix>> = self.graph.node_indices().collect();
        for node in all_nodes.iter() {
            if *node != self.start_node
                && *node != self.end_node
                && self.graph.neighbors_directed(*node, Incoming).count() == 0
            {
                self.graph
                    .add_edge(self.start_node, *node, POAEdgeData::new_for_start_or_end());
            }
        }

        // Connect nodes with no outgoing edges to the end node
        for node in all_nodes.iter() {
            if *node != self.end_node
                && *node != self.start_node
                && self.graph.neighbors_directed(*node, Outgoing).count() == 0
            {
                self.graph
                    .add_edge(*node, self.end_node, POAEdgeData::new_for_start_or_end());
            }
        }

        self.topological_sorted = rev_postorder_nodes(self);

        // Update the rank field on each node so that node_rank() returns the correct value.
        for (rank, &node) in self.topological_sorted.iter().enumerate() {
            self.graph[node].rank = Ix::new(rank);
        }

        // Relax depth_min / depth_max in topological order. Start sentinel
        // stays at (0, 0); each other node takes min/max over its
        // predecessors' depths plus one.
        self.graph[self.start_node].depth_min = Ix::new(0);
        self.graph[self.start_node].depth_max = Ix::new(0);
        for &v in &self.topological_sorted {
            if v == self.start_node {
                continue;
            }
            let mut dmin = usize::MAX;
            let mut dmax = 0usize;
            let mut has_pred = false;
            let mut preds = self.graph.neighbors_directed(v, Incoming).detach();
            while let Some(u) = preds.next_node(&self.graph) {
                has_pred = true;
                let du_min = self.graph[u].depth_min.index();
                let du_max = self.graph[u].depth_max.index();
                if du_min + 1 < dmin {
                    dmin = du_min + 1;
                }
                if du_max + 1 > dmax {
                    dmax = du_max + 1;
                }
            }
            if !has_pred {
                dmin = 0;
            }
            self.graph[v].depth_min = Ix::new(dmin);
            self.graph[v].depth_max = Ix::new(dmax);
        }

        Ok(())
    }

    /// Number of real nodes on the shortest start→end path.
    pub fn l_min_real(&self) -> usize {
        self.graph[self.end_node]
            .depth_min
            .index()
            .saturating_sub(1)
    }

    /// Number of real nodes on the longest start→end path.
    pub fn l_max_real(&self) -> usize {
        self.graph[self.end_node]
            .depth_max
            .index()
            .saturating_sub(1)
    }
}

impl<T> GraphNodeId for T
where
    T: IndexType,
{
    #[inline(always)]
    fn index(&self) -> usize {
        self.index()
    }
}

impl<Ix> GraphBase for POAGraph<Ix>
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    type NodeIter<'a> = NodeIndices<'a, POANodeData<Ix>, Ix>;
    type Successors<'a> = Neighbors<'a, POAEdgeData, Ix>;
    type Predecessors<'a> = Neighbors<'a, POAEdgeData, Ix>;

    #[inline]
    fn all_nodes_iter(&self) -> Self::NodeIter<'_> {
        self.graph.node_indices()
    }

    #[inline]
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.graph.node_count() <= 2 // Graph with only start and end nodes is considered empty
    }

    #[inline]
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.graph.neighbors(node)
    }

    #[inline]
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }

    #[inline]
    fn out_degree(&self, node: Self::NodeType) -> usize {
        self.graph.neighbors(node).count()
    }

    #[inline]
    fn in_degree(&self, node: Self::NodeType) -> usize {
        self.graph.neighbors_directed(node, Incoming).count()
    }
}

impl<'a, Ix> GraphBase for &'a POAGraph<Ix>
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    type NodeIter<'b>
        = NodeIndices<'b, POANodeData<Ix>, Ix>
    where
        'a: 'b;
    type Successors<'b>
        = Neighbors<'b, POAEdgeData, Ix>
    where
        'a: 'b;
    type Predecessors<'b>
        = Neighbors<'b, POAEdgeData, Ix>
    where
        'a: 'b;

    #[inline]
    fn all_nodes_iter(&self) -> Self::NodeIter<'_> {
        self.graph.node_indices()
    }

    #[inline]
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.graph.node_count() <= 2 // Graph with only start and end nodes is considered empty
    }

    #[inline]
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.graph.neighbors(node)
    }

    #[inline]
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }

    #[inline]
    fn out_degree(&self, node: Self::NodeType) -> usize {
        self.graph.neighbors(node).count()
    }

    #[inline]
    fn in_degree(&self, node: Self::NodeType) -> usize {
        self.graph.neighbors_directed(node, Incoming).count()
    }
}

impl<Ix> GraphWithNodeOrdering for POAGraph<Ix>
where
    Ix: IndexType,
{
    #[inline(always)]
    fn start_node(&self) -> Self::NodeType {
        self.start_node
    }

    #[inline(always)]
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }

    #[inline(always)]
    fn node_rank(&self, node: Self::NodeType) -> usize {
        self.graph[node].rank.index()
    }

    #[inline(always)]
    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType {
        self.topological_sorted[node_rank]
    }
}

impl<Ix> GraphWithNodeOrdering for &POAGraph<Ix>
where
    Ix: IndexType,
{
    #[inline(always)]
    fn start_node(&self) -> Self::NodeType {
        self.start_node
    }

    #[inline(always)]
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }

    #[inline(always)]
    fn node_rank(&self, node: Self::NodeType) -> usize {
        self.graph[node].rank.index()
    }

    #[inline(always)]
    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType {
        self.topological_sorted[node_rank]
    }
}

impl<Ix> AlignableGraph for POAGraph<Ix>
where
    Ix: IndexType,
{
    type Node = POANodeIndex<Ix>;

    fn get_node_symbol(&self, p: Self::Node) -> u8 {
        self.graph[p].symbol
    }
}

impl<Ix> AlignableGraph for &POAGraph<Ix>
where
    Ix: IndexType,
{
    type Node = POANodeIndex<Ix>;

    fn get_node_symbol(&self, p: Self::Node) -> u8 {
        self.graph[p].symbol
    }
}

/// Update the graph with an aligned sequence to the graph
impl<Ix> AddAlignment<AlignResult<POAGraph<Ix>>> for POAGraph<Ix>
where
    Ix: IndexType,
{
    type Error = PoastaError<Ix>;

    fn add_alignment(
        &mut self,
        sequence_name: &str,
        sequence: &[u8],
        alignment_opt: Option<&AlignResult<POAGraph<Ix>>>,
        weights: &[usize],
    ) -> Result<(), PoastaError<Ix>> {
        if sequence.len() != weights.len() {
            return Err(PoastaError::WeightsUnequalSize(
                sequence.len(),
                weights.len(),
            ));
        }

        if alignment_opt.is_none() {
            // No aligned bases, just add unaligned nodes
            let (nfirst, _) = self
                .add_nodes_for_sequence(sequence, weights, 0, sequence.len())
                .unwrap();
            self.sequences
                .push(Sequence(sequence_name.to_owned(), nfirst));
            self.post_process()?;

            return Ok(());
        }

        let aln_result = alignment_opt.unwrap();

        // Check start and end of alignment
        let valid_ix: Vec<usize> = aln_result
            .alignment
            .iter()
            .filter_map(|e| e.query_pos())
            .filter(|qpos| *qpos < sequence.len())
            .collect();

        if valid_ix.is_empty() {
            return Err(PoastaError::InvalidAlignment);
        }

        // Add unaligned bases
        let first = valid_ix.first().unwrap();
        let last = valid_ix.last().unwrap();

        let mut nodes_unaligned_begin = self.add_nodes_for_sequence(sequence, weights, 0, *first);

        let mut prev = if let Some((_, begin_n2)) = nodes_unaligned_begin {
            Some(begin_n2)
        } else {
            None
        };

        let nodes_unaligned_end =
            self.add_nodes_for_sequence(sequence, weights, last + 1, sequence.len());

        // Add aligned bases
        for AlignedPair(rpos, qpos) in &aln_result.alignment {
            if qpos.is_none() {
                continue;
            }

            let q = qpos.unwrap();
            let mut curr: Option<POANodeIndex<Ix>> = None;
            let qsymbol = sequence[q];

            if let Some(r) = rpos {
                // We got an aligned pair
                let rsymbol = self.graph[*r].symbol;
                if rsymbol == qsymbol {
                    curr = Some(*r);
                } else {
                    // Aligned to a node with a different symbol
                    // Check if that node is already aligned to other nodes in the graph with that symbol
                    for other_ix in &self.graph[*r].aligned_nodes {
                        if self.graph[*other_ix].symbol == qsymbol {
                            curr = Some(*other_ix);
                            break;
                        }
                    }

                    if curr.is_none() {
                        // Even the selected node does not have any matching aligning nodes, create a new node with this symbol
                        let new_node = self.graph.add_node(POANodeData::new(qsymbol));
                        curr = Some(new_node);

                        // Add this new node to the `aligned_nodes` in the other existing nodes
                        let other_nodes = self.graph[*r].aligned_nodes.clone();
                        for other_ix in &other_nodes {
                            self.graph[*other_ix].aligned_nodes.push(new_node);
                            self.graph[new_node].aligned_nodes.push(*other_ix);
                        }

                        self.graph[*r].aligned_nodes.push(new_node);
                        self.graph[new_node].aligned_nodes.push(*r);
                    }
                }
            } else {
                // It's an insertion
                let new_node = self.graph.add_node(POANodeData::new(qsymbol));
                curr = Some(new_node);
            }

            if nodes_unaligned_begin.is_none() {
                nodes_unaligned_begin = Some((curr.unwrap(), curr.unwrap()));
            }

            // `curr` should be set by now. Add edge from previous node if exists.
            if let Some(ref p) = prev {
                self.add_edge(
                    *p,
                    curr.unwrap(),
                    self.sequences.len(),
                    weights[q - 1] + weights[q],
                );
            }

            prev = curr;
        }

        if let Some((unaligned_end, _)) = nodes_unaligned_end {
            self.add_edge(
                prev.unwrap(),
                unaligned_end,
                self.sequences.len(),
                weights[*last] + weights[*last + 1],
            );
        }

        self.sequences.push(Sequence(
            sequence_name.to_owned(),
            nodes_unaligned_begin.unwrap().0,
        ));

        self.post_process()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn l_min_l_max_empty_graph() {
        let g: POAGraph<u32> = POAGraph::new();
        assert_eq!(g.l_min_real(), 0);
        assert_eq!(g.l_max_real(), 0);
    }

    #[test]
    fn l_min_l_max_linear_graph() {
        let mut g: POAGraph<u32> = POAGraph::new();
        g.add_alignment("s0", b"ACGTA", None, &[1; 5]).unwrap();
        assert_eq!(g.l_min_real(), 5);
        assert_eq!(g.l_max_real(), 5);
    }

    #[test]
    fn l_min_l_max_bubble_graph() {
        // Two paths through a bubble: one of length 4, one of length 5.
        let mut g: POAGraph<u32> = POAGraph::new();
        g.add_alignment("s0", b"ACGTA", None, &[1; 5]).unwrap();
        g.add_alignment("s1", b"ACGA", None, &[1; 4]).unwrap();
        assert_eq!(g.l_min_real(), 4);
        assert_eq!(g.l_max_real(), 5);
    }
}
