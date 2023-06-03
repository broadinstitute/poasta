use petgraph::graph::NodeIndex;

#[derive(Clone, Debug)]
pub struct AlignedPair {
    /// Represents the node rank in the graph
    pub rpos: Option<NodeIndex>,

    /// Query sequence position
    pub qpos: Option<usize>
}

impl AlignedPair {
    pub fn new(rpos: Option<NodeIndex>, qpos: Option<usize>) -> Self {
        Self { rpos, qpos }
    }

    pub fn is_aligned(&self) -> bool {
        matches!((self.rpos, self.qpos), (Some(_), Some(_)))
    }

    pub fn is_indel(&self) -> bool {
        !self.is_aligned()
    }
}

pub type Alignment = Vec<AlignedPair>;
