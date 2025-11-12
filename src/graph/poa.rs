pub use petgraph::graph::IndexType;
use petgraph::stable_graph::{Neighbors, NodeIndices};
use petgraph::visit::EdgeRef;
use petgraph::{Incoming, Outgoing};

use crate::errors::{GraphError, PoastaError};
use crate::graph::alignment::{AddAlignment, SeqAlignedPair, SeqAlignment};

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
        pub aligned_nodes: Vec<Ix>,
        pub rank: Ix,
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

use crate::align::traits::{AlignableGraph, AlignableGraphNodePos};
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
        let mut graph = graph_impl::POAGraphImpl::new();
        let start_node = graph.add_node(POANodeData::new(b'#'));
        let end_node = graph.add_node(POANodeData::new(b'$'));
        let initial_edge =
            graph.add_edge(start_node, end_node, POAEdgeData::new_for_start_or_end());

        Self {
            graph,
            sequences: Vec::default(),
            topological_sorted: vec![start_node, end_node],
            start_node,
            end_node,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.graph.node_count() == 2 // Graph with only start and end nodes is considered empty
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
                .add_edge(s, t, POAEdgeData::new(sequence_id, weight));
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

        self.topological_sorted = rev_postorder_nodes(&self.graph);

        Ok(())
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
        self.graph[node].rank
    }

    #[inline(always)]
    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType {
        self.topological_sorted[node_rank]
    }
}


/// Update the graph with an aligned sequence to the graph
impl<Ix> AddAlignment<SeqAlignment<Ix>> for POAGraph<Ix>
where 
    Ix: IndexType,
{
    type Error = PoastaError<Ix>;
    
    fn add_alignment(
        &mut self, 
        sequence_name: &str, 
        sequence: &[u8], 
        alignment_opt: Option<&SeqAlignment<Ix>>, 
        weights: &[usize]
    ) -> Result<(), PoastaError<Ix>> {
        if sequence.len() != weights.len() {
            return Err(GraphError::WeightsUnequalSize(
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

        let alignment = alignment_opt.unwrap();

        // Check start and end of alignment
        let valid_ix: Vec<usize> = alignment
            .iter()
            .filter_map(|e| e.qpos)
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
        for SeqAlignedPair { rpos, qpos } in alignment {
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