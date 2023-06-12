use petgraph::prelude::*;
use petgraph::algo::toposort;

use crate::alignment::{Alignment, AlignedPair};

pub enum NodeByRank {
    Start,
    Node(NodeIndex)
}


#[derive(Debug)]
pub struct POANodeData {
    pub symbol: u8,
    pub aligned_nodes: Vec<NodeIndex>,
    pub rank: usize
}

impl POANodeData {
    fn new(symbol: u8) -> Self {
        POANodeData {
            symbol,
            aligned_nodes: Vec::new(),
            rank: 0
        }
    }
}

#[derive(Debug)]
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
}

type POAGraphType = DiGraph<POANodeData, POAEdgeData>;

pub struct POAGraph {
    pub graph: POAGraphType,
    pub sequences: Vec<NodeIndex>,
    topological_sorted: Vec<NodeIndex>,
    start_nodes: Vec<NodeIndex>,
    end_nodes: Vec<NodeIndex>
}

impl POAGraph {
    pub fn new() -> Self {
        POAGraph {
            graph: POAGraphType::new(),
            sequences: Vec::new(),
            topological_sorted: Vec::new(),
            start_nodes: Vec::new(),
            end_nodes: Vec::new(),
        }
    }

    fn add_edge(&mut self, s: NodeIndex, t: NodeIndex, sequence_id: usize, weight: usize) {
        // If edge exists, update sequence ID and weight of the existing one
        if let Some(e) = self.graph.find_edge(s, t) {
            let mut edge_data = self.graph.edge_weight_mut(e).unwrap();
            edge_data.sequence_ids.push(self.sequences.len());
            edge_data.weight += weight;
        } else {
            self.graph.add_edge(s, t, POAEdgeData::new(sequence_id, weight));
        }
    }

    pub fn add_nodes_for_sequence(
        &mut self,
        sequence: &[u8],
        weights: &[usize],
        start: usize,
        end: usize,
    ) -> Option<(NodeIndex, NodeIndex)> {
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
                self.add_edge(prev_node, curr_node, self.sequences.len(), weights[pos - 1] + weights[pos])
            }

            prev = Some(curr_node)
        }

        Some((first_node.unwrap(), prev.unwrap()))
    }

    pub fn add_alignment_with_weights(
        &mut self,
        sequence: &[u8],
        alignment_opt: Option<&Alignment>,
        weights: &[usize]
    ) -> Result<(), String> {
        if sequence.len() != weights.len() {
            return Err(format!(
                "Sequence length ({}) is not the same as the number of weights ({})!",
                sequence.len(), weights.len()
            ));
        }

        if alignment_opt.is_none() {
            // No aligned bases, just add unaligned nodes
            let (nfirst, _) = self.add_nodes_for_sequence(
                sequence, weights, 0, sequence.len()).unwrap();
            self.sequences.push(nfirst);
            self.post_process()?;

            return Ok(())
        }

        let alignment = alignment_opt.unwrap();

        // Check start and end of alignment
        let valid_ix: Vec<usize> = alignment.iter()
            .filter_map(|e| e.qpos)
            .filter(|qpos| *qpos < sequence.len()).collect();

        if valid_ix.is_empty() {
            return Err(String::from("No valid aligned positions!"));
        }

        // Add unaligned bases
        let first = valid_ix.first().unwrap();
        let last = valid_ix.last().unwrap();

        let mut nodes_unaligned_begin = self.add_nodes_for_sequence(
            sequence, weights, 0, *first);

        let mut prev = if let Some((_, begin_n2)) = nodes_unaligned_begin {
            Some(begin_n2)
        } else {
            None
        };

        let nodes_unaligned_end = self.add_nodes_for_sequence(
            sequence, weights, last+1, sequence.len());

        // Add aligned bases
        for AlignedPair {rpos, qpos} in alignment {
            if qpos.is_none() {
                continue;
            }

            let q = qpos.unwrap();
            let mut curr: Option<NodeIndex> = None;
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
                let new_node = self.graph.add_node(POANodeData::new(qsymbol.clone()));
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

        self.sequences.push(nodes_unaligned_begin.unwrap().0);

        self.post_process()?;

        Ok(())
    }

    fn post_process(&mut self) -> Result<(), String> {
        self.topological_sorted.clear();
        self.start_nodes.clear();

        self.topological_sorted = toposort(&self.graph, None)
            .map_err(|_| String::from("Graph contains a cycle!"))?;

        for (rank, node) in self.topological_sorted.iter().enumerate() {
            self.graph[*node].rank = rank + 1; // Rank zero is reserved for the special "start" node

            if self.graph.neighbors_directed(*node, Incoming).count() == 0 {
                self.start_nodes.push(*node);
            }

            if self.graph.neighbors(*node).count() == 0 {
                self.end_nodes.push(*node);
            }
        }

        Ok(())
    }

    pub fn max_rank(&self) -> usize {
        self.graph.node_count() + 1
    }

    pub fn get_node_by_rank(&self, rank: usize) -> NodeByRank {
        if rank == 0 {
            NodeByRank::Start
        } else {
            NodeByRank::Node(self.topological_sorted[rank - 1])
        }
    }

    pub fn get_node_rank(&self, node: NodeIndex) -> usize {
        self.graph[node].rank
    }

    pub fn neighbors_for_rank(&self, rank: usize) -> Box<dyn Iterator<Item=NodeIndex> + '_> {
        match self.get_node_by_rank(rank) {
            NodeByRank::Start => Box::new(self.start_nodes.clone().into_iter()),
            NodeByRank::Node(node) => Box::new(self.graph.neighbors(node))
        }
    }

    pub fn is_neighbor_rank(&self, from_rank: usize, to_canditate_rank: usize) -> bool {
        match self.get_node_by_rank(from_rank) {
            NodeByRank::Start => self.start_nodes.iter().any(|succ| self.graph[*succ].rank == to_canditate_rank),
            NodeByRank::Node(node) => self.graph.neighbors(node).any(|succ| self.graph[succ].rank == to_canditate_rank)
        }
    }

    pub fn is_symbol_equal(&self, rank: usize, symbol: u8) -> bool {
        match self.get_node_by_rank(rank) {
            NodeByRank::Start => false,
            NodeByRank::Node(node) => self.graph[node].symbol == symbol
        }
    }

    pub fn start_nodes(&self) -> &Vec<NodeIndex> {
        &self.start_nodes
    }

    pub fn end_nodes(&self) -> &Vec<NodeIndex> {
        &self.end_nodes
    }
}
