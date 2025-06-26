use std::fmt;
use crate::graph::traits::{GraphNodeId, GraphBase, GraphWithStartEnd};

pub trait AlignableGraphNodePos: fmt::Debug + Clone + Copy + PartialEq + Eq {
    type NodeType: GraphNodeId;
    
    fn new(node: Self::NodeType, pos: usize) -> Self;
    fn node(&self) -> Self::NodeType;
    fn pos(&self) -> usize;
}

pub trait AlignableGraph: GraphWithStartEnd {
    type NodePosType: AlignableGraphNodePos<NodeType = Self::NodeType>;
    
    fn node_seq(&self, node: <Self as GraphBase>::NodeType) -> &[u8];

    fn node_len(&self, node: <Self as GraphBase>::NodeType) -> usize {
        self.node_seq(node).len()
    }
    
    fn get_node_symbol(&self, p: Self::NodePosType) -> u8 {
        self.node_seq(p.node())[p.pos()]
    }
}