use std::fmt::{Display, Formatter};
use petgraph::graph::{IndexType as PetgraphIndexType};
use petgraph::{Incoming, Outgoing};
use petgraph::algo::toposort;
use petgraph::prelude::{NodeIndex, StableDiGraph};
use petgraph::stable_graph::{Neighbors, NodeIndices};
use petgraph::visit::{GraphBase, EdgeRef};

use serde::{Deserialize, Serialize};
use serde::de::DeserializeOwned;

use crate::errors::PoastaError;
use crate::aligner::alignment::{AlignedPair, Alignment};
use crate::graphs::{AlignableRefGraph, NodeIndexType};
use crate::io::graph::format_as_dot;

/// A sequence aligned to the POA graph.
///
/// Stores the sequence name and the start node in the graph.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sequence<N>(pub(crate) String, pub(crate) N)
where
    N: NodeIndexType
;

impl<N> Sequence<N>
where
    N: NodeIndexType
{
    pub fn name(&self) -> &String {
        &self.0
    }

    pub fn start_node(&self) -> N {
        self.1
    }
}


#[derive(Debug, Serialize, Deserialize)]
pub struct POANodeData<N>
where
    N: NodeIndexType
{
    pub symbol: u8,
    pub aligned_nodes: Vec<N>,
}

impl<N> POANodeData<N>
where
    N: NodeIndexType
{
    pub(crate) fn new(symbol: u8) -> Self {
        POANodeData {
            symbol,
            aligned_nodes: Vec::new(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct POAEdgeData {
    pub weight: usize,
    pub sequence_ids: Vec<usize>,
}

impl POAEdgeData {
    fn new(sequence_id: usize, weight: usize) -> Self {
        POAEdgeData {
            weight,
            sequence_ids: vec![sequence_id],
        }
    }

    fn new_for_start_or_end() -> Self {
        POAEdgeData {
            weight: 0, sequence_ids: vec![]
        }
    }
}

pub type POAGraphType<Ix> = StableDiGraph<POANodeData<NodeIndex<Ix>>, POAEdgeData, Ix>;
pub(crate) type POANodeIndex<Ix> = <POAGraphType<Ix> as GraphBase>::NodeId;

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct POAGraph<Ix = u32>
where
    Ix: PetgraphIndexType
{
    pub(crate) graph: POAGraphType<Ix>,
    pub sequences: Vec<Sequence<POANodeIndex<Ix>>>,
    topological_sorted: Vec<POANodeIndex<Ix>>,
    start_node: POANodeIndex<Ix>,
    end_node: POANodeIndex<Ix>,
}

impl<Ix> POAGraph<Ix>
where
    Ix: PetgraphIndexType + DeserializeOwned
{
    pub fn new() -> Self {
        let mut graph = POAGraphType::<Ix>::default();
        let start_node = graph.add_node(POANodeData::new(b'#'));
        let end_node = graph.add_node(POANodeData::new(b'$'));

        POAGraph {
            graph,
            sequences: Vec::new(),
            topological_sorted: Vec::new(),
            start_node,
            end_node,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.node_count() == 0
    }

    pub(crate) fn add_edge(&mut self, s: POANodeIndex<Ix>, t: POANodeIndex<Ix>, sequence_id: usize, weight: usize) {
        // If edge exists, update sequence ID and weight of the existing one
        if let Some(e) = self.graph.find_edge(s, t) {
            let mut edge_data = self.graph.edge_weight_mut(e).unwrap();
            edge_data.sequence_ids.push(sequence_id);
            edge_data.weight += weight;
        } else {
            self.graph.add_edge(s, t, POAEdgeData::new(sequence_id, weight));
        }
    }

    pub fn add_nodes_for_sequence<T: AsRef<[u8]>>(
        &mut self,
        sequence: T,
        weights: &[usize],
        start: usize,
        end: usize,
    ) -> Option<(POANodeIndex<Ix>, POANodeIndex<Ix>)> {
        let seq = sequence.as_ref();

        if start == end {
            return None;
        }

        let mut first_node = None;
        let mut prev = None;
        for pos in start..end {
            let curr_node = self.graph.add_node(POANodeData::new(seq[pos]));

            if first_node.is_none() {
                first_node = Some(curr_node);
            }

            if let Some(prev_node) = prev {
                self.add_edge(prev_node, curr_node, self.sequences.len(), weights[pos - 1] + weights[pos])
            }

            prev = Some(curr_node)
        }

        Some((first_node.unwrap(), prev.unwrap()))
    }

    pub fn add_alignment_with_weights<T: AsRef<[u8]>>(
        &mut self,
        sequence_name: &str,
        sequence: T,
        alignment_opt: Option<&Alignment<POANodeIndex<Ix>>>,
        weights: &[usize]
    ) -> Result<(), PoastaError> {
        let seq = sequence.as_ref();

        if seq.len() != weights.len() {
            return Err(PoastaError::WeightsUnequalSize(seq.len(), weights.len()))
        }

        if alignment_opt.is_none() {
            // No aligned bases, just add unaligned nodes
            let (nfirst, _) = self.add_nodes_for_sequence(
                seq, weights, 0, seq.len()).unwrap();
            self.sequences.push(Sequence(sequence_name.to_owned(), nfirst));
            self.post_process()?;

            return Ok(())
        }

        let alignment = alignment_opt.unwrap();

        // Check start and end of alignment
        let valid_ix: Vec<usize> = alignment.iter()
            .filter_map(|e| e.qpos)
            .filter(|qpos| *qpos < seq.len()).collect();

        if valid_ix.is_empty() {
            return Err(PoastaError::InvalidAlignment);
        }

        // Add unaligned bases
        let first = valid_ix.first().unwrap();
        let last = valid_ix.last().unwrap();

        let mut nodes_unaligned_begin = self.add_nodes_for_sequence(
            seq, weights, 0, *first);

        let mut prev = if let Some((_, begin_n2)) = nodes_unaligned_begin {
            Some(begin_n2)
        } else {
            None
        };

        let nodes_unaligned_end = self.add_nodes_for_sequence(
            seq, weights, last+1, seq.len());

        // Add aligned bases
        for AlignedPair {rpos, qpos} in alignment {
            if qpos.is_none() {
                continue;
            }

            let q = qpos.unwrap();
            let mut curr: Option<POANodeIndex<Ix>> = None;
            let qsymbol = seq[q];

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
                        let new_node = self.graph.add_node(POANodeData::new(qsymbol.clone()));
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
                self.add_edge(*p, curr.unwrap(), self.sequences.len(), weights[q - 1] + weights[q]);
            }

            prev = curr;
        }

        if let Some((unaligned_end, _)) = nodes_unaligned_end {
            self.add_edge(prev.unwrap(), unaligned_end, self.sequences.len(), weights[*last] + weights[*last + 1]);
        }

        self.sequences.push(Sequence(sequence_name.to_owned(), nodes_unaligned_begin.unwrap().0));

        self.post_process()?;

        Ok(())
    }

    pub(crate) fn post_process(&mut self) -> Result<(), PoastaError> {
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
        let all_nodes: Vec<NodeIndex<Ix>> = self.graph.node_indices().collect();
        for node in all_nodes.iter() {
            if *node != self.start_node
                && *node != self.end_node
                && self.graph.neighbors_directed(*node, Incoming).count() == 0
            {
                self.graph.add_edge(self.start_node, *node, POAEdgeData::new_for_start_or_end());
            }
        }

        // Connect nodes with no outgoing edges to the end node
        for node in all_nodes.iter() {
            if *node != self.end_node
                && *node != self.start_node
                && self.graph.neighbors_directed(*node, Outgoing).count() == 0
            {
                self.graph.add_edge(*node, self.end_node, POAEdgeData::new_for_start_or_end());
            }
        }

        self.topological_sorted = toposort(&self.graph, None)?;

        Ok(())
    }

    pub fn get_node_ranks(&self) -> Vec<usize> {
        let mut ranks = vec![0; self.topological_sorted.len()];
        for (rank, node) in self.topological_sorted.iter().enumerate() {
            ranks[NodeIndexType::index(node)] = rank;
        }

        ranks
    }

}

impl<Ix> AlignableRefGraph for POAGraph<Ix>
where
    Ix: PetgraphIndexType + DeserializeOwned,
{
    type NodeIndex = POANodeIndex<Ix>;
    type NodeIterator<'a> = NodeIndices<'a, POANodeData<POANodeIndex<Ix>>, Ix>;
    type PredecessorIterator<'a> = Neighbors<'a, POAEdgeData, Ix>;
    type SuccessorIterator<'a> = Neighbors<'a, POAEdgeData, Ix>;

    #[inline]
    fn all_nodes(&self) -> Self::NodeIterator<'_> {
        self.graph.node_indices()
    }

    #[inline]
    fn node_count(&self) -> usize {
        self.graph.node_count() - 2
    }

    #[inline]
    fn node_count_with_start_and_end(&self) -> usize {
        self.graph.node_count()
    }

    #[inline]
    fn edge_count(&self) -> usize {
        // Exclude edges from start and end node
        self.graph.edge_count()
            - self.graph.neighbors_directed(self.start_node, Outgoing).count()
            - self.graph.neighbors_directed(self.end_node, Incoming).count()
    }

    #[inline]
    fn start_node(&self) -> Self::NodeIndex {
        self.start_node
    }

    #[inline]
    fn end_node(&self) -> Self::NodeIndex {
        self.end_node
    }

    #[inline]
    fn predecessors(&self, node: Self::NodeIndex) -> Self::PredecessorIterator<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }

    #[inline]
    fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_> {
        self.graph.neighbors(node)
    }

    #[inline]
    fn in_degree(&self, node: Self::NodeIndex) -> usize {
        self.graph.neighbors_directed(node, Incoming).count()
    }

    #[inline]
    fn out_degree(&self, node: Self::NodeIndex) -> usize {
        self.graph.neighbors_directed(node, Outgoing).count()
    }

    #[inline]
    fn is_end(&self, node: Self::NodeIndex) -> bool {
        self.end_node == node
    }

    #[inline]
    fn get_symbol(&self, node: Self::NodeIndex) -> char {
        char::from(self.graph[node].symbol)
    }

    #[inline]
    fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool {
        self.end_node == node || self.graph[node].symbol == symbol
    }

    #[inline]
    fn get_node_ordering(&self) -> Vec<usize> {
        self.get_node_ranks()
    }
}

impl<Ix> Display for POAGraph<Ix>
where
    Ix: PetgraphIndexType
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        format_as_dot(f, self)
    }
}

/// Wrapper with fixed node index types to make serialization to disk easier
#[derive(Serialize, Deserialize)]
pub enum POAGraphWithIx {
    U8(POAGraph<u8>),
    U16(POAGraph<u16>),
    U32(POAGraph<u32>),
    USIZE(POAGraph<usize>)
}

impl Display for POAGraphWithIx {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::U8(g) => format_as_dot(f, g),
            Self::U16(g) => format_as_dot(f, g),
            Self::U32(g) => format_as_dot(f, g),
            Self::USIZE(g) => format_as_dot(f, g)
        }
    }
}