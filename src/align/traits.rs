use std::error::Error;
use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};


pub trait AlignableGraph: 
    GraphWithNodeOrdering<NodeType = Self::Node> 
{
    type Node: GraphNodeId; // Mostly here to constrain subtrait associated types
    
    fn get_node_symbol(&self, p: Self::Node) -> u8;
}

pub trait AlignmentEngine<ToAlign> {
    type Graph: AlignableGraph;
    type Success;
    type Error: Error;
    
    fn align(&self, graph: &Self::Graph, to_align: ToAlign) -> Result<Self::Success, Self::Error>;
}
