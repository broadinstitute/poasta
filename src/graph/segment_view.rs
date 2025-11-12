//! This module offers tools to transform a standard, character-labeled POA
//! graph into a view where non-branching paths ("segments") in the graph are
//! transformed to sequence-labeled nodes, more akin to a typical genomic
//! sequence graph.

use std::collections::VecDeque;

use petgraph::stable_graph::{Neighbors, NodeIndices};
use petgraph::visit::IntoEdgeReferences;
use petgraph::Direction::Incoming;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::poa::{IndexType, POAGraph, POANodeIndex};
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};
use crate::graph::utils::rev_postorder_nodes;

mod graph_impl {
    use petgraph::prelude::StableDiGraph;
    use petgraph::visit::GraphBase;

    use crate::graph::poa::{IndexType, POANodeIndex};

    pub type SegmentGraphImpl<Ix> = StableDiGraph<SegmentNode<Ix>, (), Ix>;
    pub type SegmentGraphNodeIndex<Ix> = <SegmentGraphImpl<Ix> as GraphBase>::NodeId;

    #[derive(Debug)]
    pub struct SegmentNode<Ix>
    where
        Ix: IndexType,
    {
        pub sequence: Vec<u8>,
        pub orig_node_start: POANodeIndex<Ix>,
        pub orig_node_end: POANodeIndex<Ix>,
        pub rank: Ix,
    }

    impl<Ix> SegmentNode<Ix>
    where
        Ix: IndexType,
    {
        pub fn new(
            sequence: &[u8],
            orig_node_start: POANodeIndex<Ix>,
            orig_node_end: POANodeIndex<Ix>,
        ) -> Self {
            Self {
                sequence: Vec::from(sequence),
                orig_node_start,
                orig_node_end,
                rank: Ix::default(),
            }
        }
    }
}

pub use graph_impl::SegmentGraphNodeIndex;
use graph_impl::{SegmentGraphImpl, SegmentNode};

/// Creates a "segment" view for a character-labeled POA graph
/// 
/// The view is constructed by traversing the POA graph and detecting non-branching
/// paths, combining those characters into a single sequence-labeled node.
pub struct SegmentGraphView<'a, Ix>
where
    Ix: IndexType,
{
    poa_graph: &'a POAGraph<Ix>,
    graph: SegmentGraphImpl<Ix>,
    node_to_segment: FxHashMap<POANodeIndex<Ix>, (SegmentGraphNodeIndex<Ix>, usize)>,
    start_node: SegmentGraphNodeIndex<Ix>,
    end_node: SegmentGraphNodeIndex<Ix>,
    toposorted: Vec<SegmentGraphNodeIndex<Ix>>,
}

impl<'a, Ix> SegmentGraphView<'a, Ix>
where
    Ix: IndexType,
{
    pub fn new(poa_graph: &'a POAGraph<Ix>) -> Self {
        // Find non-branching paths, and store them as nodes ("segments") in the
        // new graph.
        let mut segment_graph = SegmentGraphImpl::new();
        let mut start_node = segment_graph.add_node(SegmentNode::new(
            b"$",
            poa_graph.start_node(),
            poa_graph.start_node(),
        ));
        
        let mut end_node = segment_graph.add_node(SegmentNode::new(
            b"#",
            poa_graph.end_node(),
            poa_graph.end_node()
        ));
        
        // Keep track of segment start and end nodes, to be used later to infer
        // which edges to add between segments
        let mut segment_starts = FxHashMap::default();
        let mut segment_ends = FxHashMap::default();
        let mut node_to_segment = FxHashMap::default();
        
        segment_starts.insert(poa_graph.start_node(), start_node);
        segment_ends.insert(poa_graph.start_node(), start_node);
        segment_starts.insert(poa_graph.end_node(), end_node);
        segment_ends.insert(poa_graph.end_node(), end_node);
        
        // Prepare traversing the POA graph, constructing segment nodes
        let mut visited = FxHashSet::default();
        let mut queue = VecDeque::new();
        visited.insert(poa_graph.start_node());
        queue.push_back(poa_graph.start_node());
        let mut curr_segment_id = segment_graph.node_count();
        while let Some(front) = queue.pop_front() {
            let mut curr_node = front;
            
            // Don't include the start node in a segment
            if front != poa_graph.start_node() {
                // Start building segment sequence
                let mut segment_seq = vec![poa_graph.node_symbol(front)];
                let mut curr_out_degree = poa_graph.out_degree(curr_node);
                let mut curr_within_seg_pos = 0usize;
                
                let new_segment_node = segment_graph.add_node(SegmentNode::new(
                    "",
                    curr_node,
                    curr_node,
                ));
                node_to_segment.insert(front, (new_segment_node, curr_within_seg_pos));
                segment_starts.insert(front, new_segment_node);
                
                // Follow successors when there's a single outgoing edge, and the
                // successor node only has a single incoming edge
                while curr_out_degree == 1 {
                    let next_node = poa_graph.successors(curr_node).next().unwrap();
                    let in_degree_next = poa_graph.in_degree(next_node);
                    curr_within_seg_pos += 1;
    
                    if in_degree_next == 1 && next_node != poa_graph.end_node() {
                        segment_seq.push(poa_graph.get_symbol(next_node));
                        node_to_segment.insert(next_node, (new_segment_node, curr_within_seg_pos));
                    } else {
                        break;
                    }
    
                    curr_node = next_node;
                    curr_out_degree = poa_graph.out_degree(curr_node);
                }
                
                segment_graph[new_segment_node].sequence = segment_seq;
                segment_graph[new_segment_node].orig_node_end = curr_node;
                segment_ends.insert(curr_node, new_segment_node);
                
                visited.insert(curr_node);
            }
            
            // Queue successor nodes
            for succ in poa_graph.successors(curr_node) {
                if !visited.contains(&succ) && succ != poa_graph.end_node() {
                    visited.insert(succ);
                    queue.push_back(succ);
                }
            }
        }
        
        // Add edges between segments
        for edge in poa_graph.graph.edge_references() {
            if segment_ends.contains_key(&edge.source()) && segment_starts.contains_key(&edge.target()) {
                let src = segment_ends[&edge.source()];
                let target = segment_starts[&edge.target()];
                
                segment_graph.add_edge(src, target, ());
            }
        }
        
        let mut segment_view = Self {
            poa_graph,
            graph: segment_graph,
            node_to_segment,
            start_node,
            end_node,
            toposorted: Vec::default(),
        };
        
        segment_view.toposorted = rev_postorder_nodes(segment_graph);
        segment_view
    }
}


impl<'a, Ix> GraphBase for SegmentGraphView<'a, Ix> 
where 
    Ix: IndexType,
{
    type NodeType = SegmentGraphNodeIndex<Ix>;
    type NodeIter<'b> = NodeIndices<'b, SegmentGraphNodeIndex<Ix>, Ix>
        where Self: 'b;
    type Successors<'b> = Neighbors<'b, (), Ix>
        where Self: 'b;
    type Predecessors<'b> = Neighbors<'b, (), Ix>
            where Self: 'b;
    
    fn all_nodes_iter(&self) -> Self::NodeIter<'_> {
        self.graph.node_indices()
    }
    
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }
    
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.graph.neighbors(node)
    }
    
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.graph.neighbors_directed(node, Incoming)
    }
    
    fn out_degree(&self, node: Self::NodeType) -> usize {
        self.graph.neighbors(node).count()
    }
    
    fn in_degree(&self, node: Self::NodeType) -> usize {
        self.graph.neighbors_directed(node, Incoming).count()
    }
}

impl<'a, Ix> GraphWithNodeOrdering for SegmentGraphView<'a, Ix>
where 
    Ix: IndexType,
{
    fn start_node(&self) -> Self::NodeType {
        self.start_node
    }
    
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }
    
    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType {
        self.toposorted[node_rank]
    }
    
    fn node_rank(&self, node: Self::NodeType) -> usize {
        self.graph[node].rank
    }
}