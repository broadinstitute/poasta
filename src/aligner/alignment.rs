use crate::graphs::NodeIndexType;

#[derive(Clone, Debug)]
pub struct AlignedPair<N>
where
    N: NodeIndexType
{
    /// Represents the node rank in the graph
    pub rpos: Option<N>,

    /// Query sequence position
    pub qpos: Option<usize>
}

impl<N> AlignedPair<N>
where
    N: NodeIndexType
{
    pub fn new(rpos: Option<N>, qpos: Option<usize>) -> Self {
        Self { rpos, qpos }
    }

    pub fn is_aligned(&self) -> bool {
        matches!((self.rpos, self.qpos), (Some(_), Some(_)))
    }

    pub fn is_indel(&self) -> bool {
        !self.is_aligned()
    }
}

pub type Alignment<N> = Vec<AlignedPair<N>>;
