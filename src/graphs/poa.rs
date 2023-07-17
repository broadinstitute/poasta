use std::fmt::{Display, Formatter};
use petgraph::graph::{IndexType as PetgraphIndexType};
use petgraph::Incoming;
use petgraph::algo::toposort;
use petgraph::prelude::{NodeIndex, StableDiGraph};
use petgraph::stable_graph::Neighbors;
use petgraph::visit::GraphBase;

use serde::{Deserialize, Serialize};
use serde::de::DeserializeOwned;

use crate::errors::PoastaError;
use crate::aligner::alignment::{AlignedPair, Alignment};
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::io::graph::format_as_dot;

/// A sequence aligned to the POA graph.
///
/// Stores the sequence name and the start node in the graph.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sequence<N>(String, N)
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
    pub rank: usize
}

impl<N> POANodeData<N>
where
    N: NodeIndexType
{
    fn new(symbol: u8) -> Self {
        POANodeData {
            symbol,
            aligned_nodes: Vec::new(),
            rank: 0
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

    fn new_for_start() -> Self {
        POAEdgeData {
            weight: 0, sequence_ids: vec![]
        }
    }
}

pub type POAGraphType<Ix> = StableDiGraph<POANodeData<NodeIndex<Ix>>, POAEdgeData, Ix>;
type POANodeIndex<Ix> = <POAGraphType<Ix> as GraphBase>::NodeId;

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct POAGraph<Ix = u32>
where
    Ix: PetgraphIndexType
{
    pub graph: POAGraphType<Ix>,
    pub sequences: Vec<Sequence<POANodeIndex<Ix>>>,
    topological_sorted: Vec<POANodeIndex<Ix>>,
    start_node: Vec<POANodeIndex<Ix>>,
    end_nodes: Vec<POANodeIndex<Ix>>,
}

impl<Ix> POAGraph<Ix>
where
    Ix: PetgraphIndexType + DeserializeOwned
{
    pub fn new() -> Self {
        POAGraph {
            graph: POAGraphType::<Ix>::default(),
            sequences: Vec::new(),
            topological_sorted: Vec::new(),
            start_node: Vec::new(),
            end_nodes: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.graph.node_count() == 0
    }

    fn add_edge(&mut self, s: POANodeIndex<Ix>, t: POANodeIndex<Ix>, sequence_id: usize, weight: usize) {
        // If edge exists, update sequence ID and weight of the existing one
        if let Some(e) = self.graph.find_edge(s, t) {
            let mut edge_data = self.graph.edge_weight_mut(e).unwrap();
            edge_data.sequence_ids.push(self.sequences.len());
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

    fn post_process(&mut self) -> Result<(), PoastaError> {
        self.topological_sorted.clear();

        for start_ix in self.start_node.iter() {
            self.graph.remove_node(*start_ix);
        }

        // Create a special "start" node that has outgoing edges to all other nodes without other
        // incoming edges
        let start_node = self.graph.add_node(POANodeData::new(b'#'));
        let all_nodes: Vec<NodeIndex<Ix>> = self.graph.node_indices().collect();
        for node in all_nodes.into_iter() {
            if node != start_node && self.graph.neighbors_directed(node, Incoming).count() == 0 {
                self.graph.add_edge(start_node, node, POAEdgeData::new_for_start());
            }
        }
        self.start_node.push(start_node);

        self.topological_sorted = toposort(&self.graph, None)?;

        for (rank, node) in self.topological_sorted.iter().enumerate() {
            self.graph[*node].rank = rank;

            if self.graph.neighbors(*node).count() == 0 {
                self.end_nodes.push(*node);
            }
        }

        Ok(())
    }

    pub fn max_rank(&self) -> usize {
        self.graph.node_count()
    }

    pub fn get_node_by_rank(&self, rank: usize) -> POANodeIndex<Ix> {
        self.topological_sorted[rank]
    }

    pub fn get_node_rank(&self, node: POANodeIndex<Ix>) -> usize {
        self.graph[node].rank
    }

    pub fn predecessors(&self, rank: usize) -> impl Iterator<Item=usize> + '_ {
        self.graph.neighbors_directed(self.topological_sorted[rank], Incoming)
            .map(|v| self.graph[v].rank)
    }

    pub fn successors(&self, rank: usize) -> impl Iterator<Item=usize> + '_ {
        self.graph.neighbors(self.topological_sorted[rank])
            .map(|v| self.graph[v].rank)
    }

    pub fn is_neighbor_rank(&self, from_rank: usize, to_canditate_rank: usize) -> bool {
        self.successors(from_rank).any(|rank| rank == to_canditate_rank)
    }

    pub fn end_nodes(&self) -> &Vec<POANodeIndex<Ix>> {
        &self.end_nodes
    }
}

impl<Ix> AlignableGraph for POAGraph<Ix>
where
    Ix: PetgraphIndexType + DeserializeOwned,
{
    type NodeIndex = POANodeIndex<Ix>;
    type SuccessorIterator<'a> = Neighbors<'a, POAEdgeData, Ix>;

    fn start_nodes(&self) -> &Vec<Self::NodeIndex> {
        &self.start_node
    }

    fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_> {
        self.graph.neighbors(node)
    }

    fn is_end(&self, node: Self::NodeIndex) -> bool {
        self.end_nodes.iter().any(|v| *v == node)
    }

    fn get_symbol(&self, node: Self::NodeIndex) -> char {
        char::from(self.graph[node].symbol)
    }

    fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool {
        self.graph[node].symbol == symbol
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