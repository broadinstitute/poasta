use std::fmt;
use crate::graph::traits::{GraphNodeId, GraphWithNodeLengths, GraphWithStartEnd};

pub trait AlignableGraphNodePos: fmt::Debug + Clone + Copy + PartialEq + Eq {
    type NodeType: GraphNodeId;
    
    fn new(node: Self::NodeType, pos: usize) -> Self;
    fn node(&self) -> Self::NodeType;
    fn pos(&self) -> usize;
}

pub trait AlignableGraph: 
    GraphWithStartEnd<NodeType = Self::Node> 
    + GraphWithNodeLengths<NodeType = Self::Node>
{
    type Node: GraphNodeId; // Mostly here to constrain subtrait associated types
    type NodePosType: AlignableGraphNodePos<NodeType = Self::Node>;
    
    fn node_seq(&self, node: Self::Node) -> &[u8];

    fn get_node_symbol(&self, p: Self::NodePosType) -> u8 {
        self.node_seq(p.node())[p.pos()]
    }
}