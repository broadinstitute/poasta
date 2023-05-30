use std::borrow::Borrow;
use std::collections::HashMap;
use petgraph::prelude::*;
use petgraph::algo::toposort;

pub trait Alphabet<T, C> {
    fn encode(&self, symbol: T) -> Option<C>;
    fn decode(&self, code: C) -> Option<T>;
}

pub struct ASCIIAlphabet {
    coder: Vec<Option<u8>>,
    decoder: Vec<Option<u8>>,
    num_codes: u8,
}

impl ASCIIAlphabet {
    /// Construct a new alphabet from the given ASCII symbols.
    pub fn new<C, T>(symbols: T) -> Self
        where
            C: Borrow<u8>,
            T: IntoIterator<Item = C>,
    {
        let mut alphabet = ASCIIAlphabet {
            coder: vec![None; u8::MAX as usize],
            decoder: vec![None; u8::MAX as usize],
            num_codes: 0,
        };

        for c in symbols.into_iter().map(|e| *e.borrow()) {
            // `c` is an ASCII character, use its ordinal value as index for `coder`, and assign
            // our own code to this character.
            if alphabet.coder[c as usize].is_none() {
                alphabet.coder[c as usize] = Some(alphabet.num_codes);
                alphabet.decoder[alphabet.num_codes as usize] = Some(c);
                alphabet.num_codes += 1;
            }
        }

        alphabet
    }
}

impl Alphabet<u8, u8> for ASCIIAlphabet {
    fn encode(&self, symbol: u8) -> Option<u8> {
        self.coder[symbol as usize]
    }

    fn decode(&self, code: u8) -> Option<u8> {
        self.decoder[code as usize]
    }
}

pub struct AlignedPair {
    pub rpos: Option<NodeIndex>,
    pub qpos: Option<usize>
}

impl AlignedPair {
    fn is_aligned(&self) -> bool {
        matches!((self.rpos, self.qpos), (Some(_), Some(_)))
    }

    fn is_indel(&self) -> bool {
        !self.is_aligned()
    }
}

pub type Alignment = Vec<AlignedPair>;

#[derive(Debug)]
pub struct POANode<A> {
    pub code: A,
    pub aligned_nodes: Vec<NodeIndex>,
    pub layer: usize
}

impl<A> POANode<A> {
    fn new(code: A) -> POANode<A> {
        POANode {
            code,
            aligned_nodes: Vec::new(),
            layer: 0
        }
    }
}

#[derive(Debug)]
pub struct POAEdge {
    pub weight: usize,
    pub sequence_ids: Vec<usize>,
}

impl POAEdge {
    fn new(sequence_id: usize, weight: usize) -> Self {
        POAEdge {
            weight,
            sequence_ids: vec![sequence_id],
        }
    }
}

type POAGraphType<A> = DiGraph<POANode<A>, POAEdge>;

pub struct POAGraph<A>
{
    pub graph: POAGraphType<A>,
    pub sequences: Vec<NodeIndex>,
    pub rank_to_node: Vec<NodeIndex>,
    node_to_rank: HashMap<NodeIndex, usize>,
}

impl<A: Eq + Clone> POAGraph<A>
{
    pub fn new() -> Self {
        POAGraph {
            graph: POAGraphType::new(),
            sequences: Vec::new(),
            rank_to_node: Vec::new(),
            node_to_rank: HashMap::new()
        }
    }

    fn add_edge(&mut self, s: NodeIndex, t: NodeIndex, sequence_id: usize, weight: usize) {
        // If edge exists, update sequence ID and weight of the existing one
        if let Some(e) = self.graph.find_edge(s, t) {
            let mut edge_data = self.graph.edge_weight_mut(e).unwrap();
            edge_data.sequence_ids.push(self.sequences.len());
            edge_data.weight += weight;
        } else {
            self.graph.add_edge(s, t, POAEdge::new(sequence_id, weight));
        }
    }

    pub fn add_nodes_for_sequence(
        &mut self,
        sequence: &[A],
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
            let c = &sequence[pos];
            let curr_node = self.graph.add_node(POANode::new(c.clone()));

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
        sequence: &[A],
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
            let (nfirst, nlast) = self.add_nodes_for_sequence(
                sequence, weights, 0, sequence.len()).unwrap();
            self.sequences.push(nfirst);
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
            let qcode = &sequence[q];

            if let Some(r) = rpos {
                // We got an aligned pair
                let ref_code = &self.graph[*r].code;
                if *ref_code == *qcode {
                    curr = Some(*r);
                } else {
                    // Aligned to a node with a different symbol
                    // Check if that node is already aligned to other nodes in the graph with that symbol
                    for other_ix in &self.graph[*r].aligned_nodes {
                        if self.graph[*other_ix].code == *qcode {
                            curr = Some(*other_ix);
                            break;
                        }
                    }

                    if curr.is_none() {
                        // Even the selected node does not have any matching aligning nodes, create a new node with this symbol
                        let new_node = self.graph.add_node(POANode::new(qcode.clone()));
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
                let new_node = self.graph.add_node(POANode::new(qcode.clone()));
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

        self.topological_sort()?;

        Ok(())
    }

    fn topological_sort(&mut self) -> Result<(), String> {
        self.rank_to_node.clear();
        self.node_to_rank.clear();

        self.rank_to_node = toposort(&self.graph, None)
            .map_err(|_| String::from("Graph contains a cycle!"))?;

        self.node_to_rank = self.rank_to_node.iter().enumerate().map(
            |(rank, node)| (node.clone(), rank)).collect();

        Ok(())
    }

    pub fn get_node_by_rank(&self, rank: usize) -> NodeIndex {
        self.rank_to_node[rank]
    }

    pub fn get_node_rank(&self, node: NodeIndex) -> usize {
        self.node_to_rank[&node]
    }

    pub fn is_neighbor(&self, n: NodeIndex, succ_candidate: NodeIndex) -> bool {
        self.graph.neighbors(n).any(|succ| succ == succ_candidate)
    }
}
