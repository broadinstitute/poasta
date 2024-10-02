
use super::astar::{AlignableGraphRef, AlignableGraphNodePos};

/// An aligned pair of residues. The first element represent
/// the position with a node of the graph, and the second element
/// represents the query sequence position.
///
/// In case of on insertion or deletion, set one of the elements to `None`.
#[derive(Debug, Clone, Copy)]
pub struct AlignedPair<P>(pub(crate) Option<P>, pub(crate) Option<usize>);

impl<P> AlignedPair<P>
where
    P: AlignableGraphNodePos,
{
    pub fn new(node_pos: Option<P>, query_pos: Option<usize>) -> Self {
        AlignedPair(node_pos, query_pos)
    }
    
    #[inline(always)]
    pub fn node_pos(&self) -> Option<P> {
        self.0
    }
    
    #[inline(always)]
    pub fn node(&self) -> Option<P::NodeType> {
        self.0.map(|pos| pos.node())
    }
    
    #[inline(always)]
    pub fn within_node_pos(&self) -> Option<usize> {
        self.0.map(|pos| pos.pos())
    }

    #[inline(always)]
    pub fn query_pos(&self) -> Option<usize> {
        self.1
    }
    
    pub fn is_aligned(&self) -> bool {
        self.0.is_some() && self.1.is_some()
    }
}


pub fn print_alignment<G>(graph: G, seq: &[u8], aln: &[AlignedPair<G::NodePosType>]) -> String
where
    G: AlignableGraphRef,
{
    let mut graph_chars = Vec::with_capacity(aln.len());
    let mut aln_chars = Vec::with_capacity(aln.len());
    let mut query_chars = Vec::with_capacity(aln.len());

    for pair in aln {
        if pair.is_aligned() {
            let node_seq = graph.node_seq(pair.node().unwrap());
            let qry = seq[pair.query_pos().unwrap()];

            let node_pos = pair.within_node_pos().unwrap();
            graph_chars.push(node_seq[node_pos]);
            aln_chars.push(if node_seq[node_pos] == qry { b'|' } else { b'*' });
            query_chars.push(qry);
        } else if let Some(gpos) = pair.node_pos() {
            let node_seq = graph.node_seq(gpos.node());
            graph_chars.push(node_seq[gpos.pos()]);
            aln_chars.push(b' ');
            query_chars.push(b'-');
        } else if let Some(qpos) = pair.query_pos() {
            let qry = seq[qpos];
            graph_chars.push(b'-');
            aln_chars.push(b' ');
            query_chars.push(qry);
        }
    }

    format!(
        "{}\n{}\n{}",
        String::from_utf8(graph_chars).unwrap(),
        String::from_utf8(aln_chars).unwrap(),
        String::from_utf8(query_chars).unwrap(),
    )
}