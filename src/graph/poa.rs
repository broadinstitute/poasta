//! The POA graph: a store that sequences are added to, and a read-only view
//! in CSR form, in topological order, that alignment and output read.
//!
//! The two have their own node types. A [`NodeId`] names a node in the store
//! and never changes. A [`POANodeIndex`] is a node's rank in the view, which
//! `post_process` renumbers after every added sequence.

use rustc_hash::FxHashMap;

use crate::align::traits::{AlignOutput, AlignableGraph, AlignedPair};
use crate::errors::PoastaError;
use crate::graph::alignment::AddAlignment;

pub use super::traits::IndexType;
use super::traits::{GraphBase, GraphWithNodeOrdering};

macro_rules! node_type {
    ($(#[$doc:meta])* $name:ident) => {
        $(#[$doc])*
        #[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
        #[repr(transparent)]
        pub struct $name<Ix = u32>(Ix);

        impl<Ix: IndexType> $name<Ix> {
            #[inline(always)]
            pub fn new(x: usize) -> Self {
                Self(Ix::new(x))
            }

            #[inline(always)]
            pub fn index(&self) -> usize {
                self.0.index()
            }
        }

        impl<Ix: IndexType> IndexType for $name<Ix> {
            #[inline(always)]
            fn new(x: usize) -> Self {
                Self(Ix::new(x))
            }

            #[inline(always)]
            fn index(&self) -> usize {
                self.0.index()
            }

            #[inline(always)]
            fn max() -> Self {
                Self(<Ix as IndexType>::max())
            }
        }
    };
}

node_type!(
    /// A node of a [`POAGraph`] as alignment sees it: its rank in the
    /// topological order. Valid until the graph changes.
    POANodeIndex
);

node_type!(
    /// A node of a [`POAGraph`] by its place in the store. Nodes are never
    /// deleted, so a `NodeId` stays valid as the graph grows.
    NodeId
);

/// What an edge of a [`POAGraph`] carries.
#[derive(Debug, Default)]
pub struct POAEdgeData {
    /// The sum of the base weights of every sequence along this edge. Each
    /// sequence adds the weights of both bases, so with weights of 1 a
    /// sequence adds 2. To count sequences, use `sequence_ids.len()`.
    pub weight: usize,
    /// The index in `POAGraph::sequences` of every sequence along this edge,
    /// in ascending order.
    pub sequence_ids: Vec<usize>,
}

/// What the edges between a sentinel and a real node carry.
static SENTINEL_EDGE: POAEdgeData = POAEdgeData {
    weight: 0,
    sequence_ids: Vec::new(),
};

#[derive(Debug)]
struct POAEdge<Ix> {
    source: NodeId<Ix>,
    target: NodeId<Ix>,
    data: POAEdgeData,
}

/// A sequence aligned to the POA graph.
///
/// Stores the sequence name and the start node in the graph.
#[derive(Debug, Clone)]
pub struct Sequence<Ix>(pub(crate) String, pub(crate) Ix)
where
    Ix: IndexType;

impl<Ix> Sequence<Ix>
where
    Ix: IndexType,
{
    pub fn name(&self) -> &String {
        &self.0
    }

    pub fn start_node(&self) -> Ix {
        self.1
    }
}

/// Lists per node in CSR form: node v has `items[start[v]..start[v + 1]]`,
/// with the edge each one comes from in `edges`. `Ix::max()` marks an edge
/// to or from a sentinel.
#[derive(Debug, Default)]
struct Adjacency<N, Ix> {
    start: Vec<usize>,
    items: Vec<N>,
    edges: Vec<Ix>,
}

impl<N: IndexType, Ix: IndexType> Adjacency<N, Ix> {
    #[inline(always)]
    fn range(&self, v: usize) -> std::ops::Range<usize> {
        self.start[v]..self.start[v + 1]
    }

    /// Size the lists for the given degrees. `cursor` gets each list's start.
    fn reset(&mut self, degrees: &[usize], cursor: &mut Vec<usize>) {
        self.start.clear();
        self.start.push(0);
        let mut total = 0;
        for &d in degrees {
            total += d;
            self.start.push(total);
        }
        self.items.clear();
        self.items.resize(total, N::default());
        self.edges.clear();
        self.edges.resize(total, <Ix as IndexType>::max());
        cursor.clear();
        cursor.extend_from_slice(&self.start[..degrees.len()]);
    }
}

/// A partial order alignment graph
///
/// A POA graph represents a multiple sequence alignment as a directed acyclic graph. Character-labeled
/// nodes represent the symbols in the input sequences and edges indicate which symbols are adjacent in
/// at least one sequence.
///
/// Each sequence used as input for graph construction can be reconstructed by tracing the path of nodes,
/// and the sequence start node is stored in `sequences`.
///
/// Nodes and edges live in a store that only grows, indexed by [`NodeId`].
/// `post_process` rebuilds the view from it, indexed by rank: successors and
/// predecessors in CSR form, the edges between the sentinels and the real
/// nodes, the symbols and the depths. The graph traits read the view, so
/// changes to the store show only after `post_process`.
#[derive(Debug)]
pub struct POAGraph<Ix = u32>
where
    Ix: IndexType,
{
    symbols: Vec<u8>,
    /// Per node, the nodes in the same MSA column.
    aligned_nodes: Vec<Vec<NodeId<Ix>>>,
    /// Edges between real nodes, in the order they were first added.
    edges: Vec<POAEdge<Ix>>,
    edge_ix: FxHashMap<(NodeId<Ix>, NodeId<Ix>), usize>,
    pub sequences: Vec<Sequence<NodeId<Ix>>>,
    start_id: NodeId<Ix>,
    end_id: NodeId<Ix>,

    node_of_rank: Vec<NodeId<Ix>>,
    rank_of: Vec<POANodeIndex<Ix>>,
    end_node: POANodeIndex<Ix>,
    rank_symbols: Vec<u8>,
    /// A node lists its newest edge first.
    succ: Adjacency<POANodeIndex<Ix>, Ix>,
    pred: Adjacency<POANodeIndex<Ix>, Ix>,
    /// Shortest and longest path length (in edges) from the start sentinel.
    depth_min: Vec<Ix>,
    depth_max: Vec<Ix>,

    /// Reused by `post_process`: the view by `NodeId`, and DFS state.
    id_succ: Adjacency<NodeId<Ix>, Ix>,
    id_pred: Adjacency<NodeId<Ix>, Ix>,
    visited_scratch: Vec<bool>,
    dfs_stack: Vec<(usize, usize)>,
    scratch: ViewScratch,
}

/// Buffers `post_process` reuses between calls.
#[derive(Debug, Default)]
struct ViewScratch {
    out_deg: Vec<usize>,
    in_deg: Vec<usize>,
    out_cur: Vec<usize>,
    in_cur: Vec<usize>,
    degrees: Vec<usize>,
    sources: Vec<usize>,
    sinks: Vec<usize>,
}

impl<Ix> Default for POAGraph<Ix>
where
    Ix: IndexType,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<Ix> POAGraph<Ix>
where
    Ix: IndexType,
{
    pub fn new() -> Self {
        let mut graph = Self {
            symbols: vec![b'#', b'$'],
            aligned_nodes: vec![Vec::new(), Vec::new()],
            edges: Vec::new(),
            edge_ix: FxHashMap::default(),
            sequences: Vec::new(),
            start_id: NodeId::new(0),
            end_id: NodeId::new(1),
            node_of_rank: Vec::new(),
            rank_of: Vec::new(),
            end_node: POANodeIndex::new(1),
            rank_symbols: Vec::new(),
            succ: Adjacency::default(),
            pred: Adjacency::default(),
            depth_min: Vec::new(),
            depth_max: Vec::new(),
            id_succ: Adjacency::default(),
            id_pred: Adjacency::default(),
            visited_scratch: Vec::new(),
            dfs_stack: Vec::new(),
            scratch: ViewScratch::default(),
        };
        graph.post_process();
        graph
    }

    #[inline]
    pub fn node_symbol(&self, node: POANodeIndex<Ix>) -> u8 {
        self.rank_symbols[node.index()]
    }

    /// The store's node at a rank.
    #[inline]
    pub fn node_id(&self, node: POANodeIndex<Ix>) -> NodeId<Ix> {
        self.node_of_rank[node.index()]
    }

    /// The rank of a node in the store.
    #[inline]
    pub fn node(&self, id: NodeId<Ix>) -> POANodeIndex<Ix> {
        self.rank_of[id.index()]
    }

    /// Every node in the store, by id, sentinels included.
    pub fn node_ids(&self) -> impl Iterator<Item = NodeId<Ix>> + use<Ix> {
        (0..self.symbols.len()).map(NodeId::new)
    }

    pub fn is_sentinel(&self, id: NodeId<Ix>) -> bool {
        id == self.start_id || id == self.end_id
    }

    /// The symbol of a node in the store, also before `post_process`.
    pub(crate) fn stored_symbol(&self, id: NodeId<Ix>) -> u8 {
        self.symbols[id.index()]
    }

    /// The nodes in the same MSA column as `id`.
    pub(crate) fn aligned_nodes(&self, id: NodeId<Ix>) -> &[NodeId<Ix>] {
        &self.aligned_nodes[id.index()]
    }

    /// Shortest and longest path length (in edges) from the start sentinel.
    pub(crate) fn depth(&self, node: POANodeIndex<Ix>) -> (usize, usize) {
        (
            self.depth_min[node.index()].index(),
            self.depth_max[node.index()].index(),
        )
    }

    fn edge_data(&self, e: Ix) -> &POAEdgeData {
        if e == <Ix as IndexType>::max() {
            &SENTINEL_EDGE
        } else {
            &self.edges[e.index()].data
        }
    }

    /// Outgoing edges of `node` in successor order, as (target, data).
    /// `node` is a rank. The edges to the end sentinel have weight 0 and no
    /// sequence ids.
    pub fn out_edges(
        &self,
        node: POANodeIndex<Ix>,
    ) -> impl Iterator<Item = (POANodeIndex<Ix>, &POAEdgeData)> + '_ {
        let r = self.succ.range(node.index());
        self.succ.items[r.clone()]
            .iter()
            .zip(&self.succ.edges[r])
            .map(|(&w, &e)| (w, self.edge_data(e)))
    }

    /// Incoming edges of `node` in predecessor order, as (source, data).
    /// `node` is a rank. The edges from the start sentinel have weight 0 and
    /// no sequence ids.
    pub fn in_edges(
        &self,
        node: POANodeIndex<Ix>,
    ) -> impl Iterator<Item = (POANodeIndex<Ix>, &POAEdgeData)> + '_ {
        let r = self.pred.range(node.index());
        self.pred.items[r.clone()]
            .iter()
            .zip(&self.pred.edges[r])
            .map(|(&u, &e)| (u, self.edge_data(e)))
    }

    /// Every edge between real nodes, as (source, target, data), in the order
    /// they were first added.
    pub(crate) fn real_edges(
        &self,
    ) -> impl Iterator<Item = (NodeId<Ix>, NodeId<Ix>, &POAEdgeData)> + '_ {
        self.edges.iter().map(|e| (e.source, e.target, &e.data))
    }

    /// Every edge, sentinel ones included, by rank of the source.
    pub(crate) fn all_edges(
        &self,
    ) -> impl Iterator<Item = (POANodeIndex<Ix>, POANodeIndex<Ix>, &POAEdgeData)> + '_ {
        self.all_nodes_iter()
            .flat_map(move |v| self.out_edges(v).map(move |(w, d)| (v, w, d)))
    }

    pub(crate) fn add_node(&mut self, symbol: u8) -> NodeId<Ix> {
        assert!(
            self.symbols.len() < <Ix as IndexType>::max().index(),
            "POAGraph::add_node: node index overflow"
        );
        self.symbols.push(symbol);
        self.aligned_nodes.push(Vec::new());
        NodeId::new(self.symbols.len() - 1)
    }

    /// Put `a` and `b` in the same MSA column.
    pub(crate) fn link_aligned(&mut self, a: NodeId<Ix>, b: NodeId<Ix>) {
        self.aligned_nodes[a.index()].push(b);
        self.aligned_nodes[b.index()].push(a);
    }

    pub(crate) fn add_edge(
        &mut self,
        s: NodeId<Ix>,
        t: NodeId<Ix>,
        sequence_id: usize,
        weight: u32,
    ) {
        // If edge exists, update sequence ID and weight of the existing one
        if let Some(&e) = self.edge_ix.get(&(s, t)) {
            let edge_data = &mut self.edges[e].data;
            edge_data.sequence_ids.push(sequence_id);
            edge_data.weight += weight as usize;
        } else {
            // `Ix::max()` marks sentinel edges, so no real edge may have it.
            assert!(
                self.edges.len() < <Ix as IndexType>::max().index(),
                "POAGraph::add_edge: edge index overflow"
            );
            self.edge_ix.insert((s, t), self.edges.len());
            self.edges.push(POAEdge {
                source: s,
                target: t,
                data: POAEdgeData {
                    weight: weight as usize,
                    sequence_ids: vec![sequence_id],
                },
            });
        }
    }

    pub fn add_nodes_for_sequence(
        &mut self,
        sequence: &[u8],
        weights: &[u32],
        start: usize,
        end: usize,
    ) -> Option<(NodeId<Ix>, NodeId<Ix>)> {
        if start == end {
            return None;
        }

        let mut first_node = None;
        let mut prev = None;
        for pos in start..end {
            let curr_node = self.add_node(sequence[pos]);

            if first_node.is_none() {
                first_node = Some(curr_node);
            }

            if let Some(prev_node) = prev {
                self.add_edge(
                    prev_node,
                    curr_node,
                    self.sequences.len(),
                    weights[pos - 1] + weights[pos],
                )
            }

            prev = Some(curr_node)
        }

        Some((first_node.unwrap(), prev.unwrap()))
    }

    /// Rebuild the view from the store. Nothing is checked: a cyclic store,
    /// which only a cyclic GFA gives, yields an order that is not
    /// topological.
    pub(crate) fn post_process(&mut self) {
        let n = self.symbols.len();
        let (start, end) = (self.start_id.index(), self.end_id.index());
        let s = &mut self.scratch;

        // The lists by NodeId first. Real nodes without incoming edges hang
        // off the start sentinel, and those without outgoing edges lead to
        // the end sentinel. The sentinels list them by descending id, as
        // petgraph did. A graph without real nodes gets one edge from start
        // to end.
        for deg in [&mut s.out_deg, &mut s.in_deg] {
            deg.clear();
            deg.resize(n, 0);
        }
        for e in &self.edges {
            s.out_deg[e.source.index()] += 1;
            s.in_deg[e.target.index()] += 1;
        }
        let is_real = |v: usize| v != start && v != end;
        s.sources.clear();
        s.sources
            .extend((0..n).rev().filter(|&v| is_real(v) && s.in_deg[v] == 0));
        s.sinks.clear();
        s.sinks
            .extend((0..n).rev().filter(|&v| is_real(v) && s.out_deg[v] == 0));
        if s.sources.is_empty() {
            s.sources.push(end);
        }
        s.out_deg[start] += s.sources.len();
        s.in_deg[end] += s.sinks.len();
        for &v in &s.sources {
            s.in_deg[v] += 1;
        }
        for &v in &s.sinks {
            s.out_deg[v] += 1;
        }
        self.id_succ.reset(&s.out_deg, &mut s.out_cur);
        self.id_pred.reset(&s.in_deg, &mut s.in_cur);
        let (id_succ, id_pred, out_cur, in_cur) = (
            &mut self.id_succ,
            &mut self.id_pred,
            &mut s.out_cur,
            &mut s.in_cur,
        );
        let mut put = |u: usize, w: usize, e: Ix| {
            id_succ.items[out_cur[u]] = NodeId::new(w);
            id_succ.edges[out_cur[u]] = e;
            out_cur[u] += 1;
            id_pred.items[in_cur[w]] = NodeId::new(u);
            id_pred.edges[in_cur[w]] = e;
            in_cur[w] += 1;
        };
        let none = <Ix as IndexType>::max();
        for &v in &s.sources {
            put(start, v, none);
        }
        for &v in &s.sinks {
            put(v, end, none);
        }
        // Newest edge first, the order petgraph's adjacency lists had. It
        // decides the topological order and which of equal-score alignments
        // the search returns.
        for (e, edge) in self.edges.iter().enumerate().rev() {
            put(edge.source.index(), edge.target.index(), Ix::new(e));
        }

        self.compute_rev_postorder();

        // The view, by rank. Each list keeps its order.
        self.rank_of.clear();
        self.rank_of.resize(n, POANodeIndex::default());
        for (rank, &id) in self.node_of_rank.iter().enumerate() {
            self.rank_of[id.index()] = POANodeIndex::new(rank);
        }
        self.end_node = self.rank_of[end];
        self.rank_symbols.clear();
        self.rank_symbols
            .extend(self.node_of_rank.iter().map(|id| self.symbols[id.index()]));
        let s = &mut self.scratch;
        s.degrees.clear();
        s.degrees
            .extend(self.node_of_rank.iter().map(|id| s.out_deg[id.index()]));
        self.succ.reset(&s.degrees, &mut s.out_cur);
        s.degrees.clear();
        s.degrees
            .extend(self.node_of_rank.iter().map(|id| s.in_deg[id.index()]));
        self.pred.reset(&s.degrees, &mut s.in_cur);
        for (rank, &id) in self.node_of_rank.iter().enumerate() {
            for (from, to) in [
                (&self.id_succ, &mut self.succ),
                (&self.id_pred, &mut self.pred),
            ] {
                let r = from.range(id.index());
                let at = to.start[rank];
                for (k, ix) in r.enumerate() {
                    to.items[at + k] = self.rank_of[from.items[ix].index()];
                    to.edges[at + k] = from.edges[ix];
                }
            }
        }

        // Relax depth_min / depth_max in rank order. Start sentinel stays at
        // (0, 0); each other node takes min/max over its predecessors'
        // depths plus one.
        let zero = Ix::new(0);
        self.depth_min.clear();
        self.depth_min.resize(n, zero);
        self.depth_max.clear();
        self.depth_max.resize(n, zero);
        for v in 1..n {
            let preds = &self.pred.items[self.pred.range(v)];
            let dmin = preds
                .iter()
                .map(|u| self.depth_min[u.index()].index() + 1)
                .min()
                .unwrap_or(0);
            let dmax = preds
                .iter()
                .map(|u| self.depth_max[u.index()].index() + 1)
                .max()
                .unwrap_or(0);
            self.depth_min[v] = Ix::new(dmin);
            self.depth_max[v] = Ix::new(dmax);
        }
    }

    /// Whether the store has no cycle. The rank order is topological exactly
    /// then, so every edge must go to a higher rank.
    pub(crate) fn is_acyclic(&self) -> bool {
        (0..self.node_of_rank.len()).all(|v| {
            self.succ.items[self.succ.range(v)]
                .iter()
                .all(|w| w.index() > v)
        })
    }

    /// Fill `node_of_rank` by iterative DFS from the start sentinel (reverse
    /// postorder) on the lists by NodeId. Nodes the DFS does not reach come
    /// last. Only a cycle without a way in from a source has them.
    fn compute_rev_postorder(&mut self) {
        let n = self.symbols.len();
        let start = self.start_id.index();
        self.visited_scratch.clear();
        self.visited_scratch.resize(n, false);
        self.node_of_rank.clear();
        self.dfs_stack.clear();

        self.visited_scratch[start] = true;
        self.dfs_stack.push((start, self.id_succ.start[start]));
        while let Some(top) = self.dfs_stack.last_mut() {
            let (v, next) = *top;
            if next < self.id_succ.start[v + 1] {
                top.1 += 1;
                let w = self.id_succ.items[next].index();
                if !self.visited_scratch[w] {
                    self.visited_scratch[w] = true;
                    self.dfs_stack.push((w, self.id_succ.start[w]));
                }
            } else {
                self.dfs_stack.pop();
                self.node_of_rank.push(NodeId::new(v));
            }
        }

        self.node_of_rank.reverse();
        let unreached = (0..n)
            .filter(|&v| !self.visited_scratch[v])
            .map(NodeId::new);
        self.node_of_rank.extend(unreached);
    }
}

impl<Ix> GraphBase for POAGraph<Ix>
where
    Ix: IndexType,
{
    type NodeType = POANodeIndex<Ix>;
    type NodeIter<'a> = std::iter::Map<std::ops::Range<usize>, fn(usize) -> POANodeIndex<Ix>>;
    type Successors<'a> = std::iter::Copied<std::slice::Iter<'a, POANodeIndex<Ix>>>;
    type Predecessors<'a> = std::iter::Copied<std::slice::Iter<'a, POANodeIndex<Ix>>>;

    #[inline]
    fn all_nodes_iter(&self) -> Self::NodeIter<'_> {
        (0..self.node_of_rank.len()).map(POANodeIndex::new)
    }

    #[inline]
    fn node_count(&self) -> usize {
        self.node_of_rank.len()
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.node_of_rank.len() <= 2 // Graph with only start and end nodes is considered empty
    }

    #[inline]
    fn successors(&self, node: Self::NodeType) -> Self::Successors<'_> {
        self.succ.items[self.succ.range(node.index())]
            .iter()
            .copied()
    }

    #[inline]
    fn predecessors(&self, node: Self::NodeType) -> Self::Predecessors<'_> {
        self.pred.items[self.pred.range(node.index())]
            .iter()
            .copied()
    }

    #[inline]
    fn out_degree(&self, node: Self::NodeType) -> usize {
        self.succ.range(node.index()).len()
    }

    #[inline]
    fn in_degree(&self, node: Self::NodeType) -> usize {
        self.pred.range(node.index()).len()
    }
}

impl<Ix> GraphWithNodeOrdering for POAGraph<Ix>
where
    Ix: IndexType,
{
    /// Always rank 0.
    #[inline(always)]
    fn start_node(&self) -> Self::NodeType {
        POANodeIndex::new(0)
    }

    #[inline(always)]
    fn end_node(&self) -> Self::NodeType {
        self.end_node
    }

    #[inline(always)]
    fn node_rank(&self, node: Self::NodeType) -> usize {
        node.index()
    }

    #[inline(always)]
    fn rank_to_node(&self, node_rank: usize) -> Self::NodeType {
        POANodeIndex::new(node_rank)
    }
}

impl<Ix> AlignableGraph for POAGraph<Ix>
where
    Ix: IndexType,
{
    type Node = POANodeIndex<Ix>;

    #[inline(always)]
    fn get_node_symbol(&self, p: Self::Node) -> u8 {
        self.rank_symbols[p.index()]
    }
}

/// Update the graph with an aligned sequence to the graph
impl<Ix> AddAlignment<AlignOutput<POAGraph<Ix>>> for POAGraph<Ix>
where
    Ix: IndexType,
{
    type Error = PoastaError<Ix>;

    fn add_alignment(
        &mut self,
        sequence_name: &str,
        sequence: &[u8],
        alignment_opt: Option<&AlignOutput<POAGraph<Ix>>>,
        weights: &[u32],
    ) -> Result<(), PoastaError<Ix>> {
        if sequence.len() != weights.len() {
            return Err(PoastaError::WeightsUnequalSize(
                sequence.len(),
                weights.len(),
            ));
        }

        if alignment_opt.is_none() {
            // No aligned bases, just add unaligned nodes
            let (nfirst, _) = self
                .add_nodes_for_sequence(sequence, weights, 0, sequence.len())
                .unwrap();
            self.sequences
                .push(Sequence(sequence_name.to_owned(), nfirst));
            self.post_process();

            return Ok(());
        }

        let aln_result = alignment_opt.unwrap();

        // Check start and end of alignment
        let valid_ix: Vec<usize> = aln_result
            .alignment
            .iter()
            .filter_map(|e| e.query_pos())
            .filter(|qpos| *qpos < sequence.len())
            .collect();

        if valid_ix.is_empty() {
            return Err(PoastaError::InvalidAlignment);
        }

        // Add unaligned bases
        let first = valid_ix.first().unwrap();
        let last = valid_ix.last().unwrap();

        let mut nodes_unaligned_begin = self.add_nodes_for_sequence(sequence, weights, 0, *first);

        let mut prev = if let Some((_, begin_n2)) = nodes_unaligned_begin {
            Some(begin_n2)
        } else {
            None
        };

        let nodes_unaligned_end =
            self.add_nodes_for_sequence(sequence, weights, last + 1, sequence.len());

        // Add aligned bases. The alignment names nodes by rank, and the view
        // stays as it is until `post_process` below.
        for AlignedPair(rpos, qpos) in &aln_result.alignment {
            if qpos.is_none() {
                continue;
            }

            let q = qpos.unwrap();
            let mut curr: Option<NodeId<Ix>> = None;
            let qsymbol = sequence[q];

            if let Some(r) = rpos.map(|r| self.node_id(r)) {
                // We got an aligned pair
                let rsymbol = self.stored_symbol(r);
                if rsymbol == qsymbol {
                    curr = Some(r);
                } else {
                    // Aligned to a node with a different symbol
                    // Check if that node is already aligned to other nodes in the graph with that symbol
                    for other_ix in self.aligned_nodes(r) {
                        if self.stored_symbol(*other_ix) == qsymbol {
                            curr = Some(*other_ix);
                            break;
                        }
                    }

                    if curr.is_none() {
                        // Even the selected node does not have any matching aligning nodes, create a new node with this symbol
                        let new_node = self.add_node(qsymbol);
                        curr = Some(new_node);

                        // Add this new node to the `aligned_nodes` in the other existing nodes
                        let other_nodes = self.aligned_nodes(r).to_vec();
                        for other_ix in other_nodes {
                            self.link_aligned(other_ix, new_node);
                        }

                        self.link_aligned(r, new_node);
                    }
                }
            } else {
                // It's an insertion
                let new_node = self.add_node(qsymbol);
                curr = Some(new_node);
            }

            if nodes_unaligned_begin.is_none() {
                nodes_unaligned_begin = Some((curr.unwrap(), curr.unwrap()));
            }

            // `curr` should be set by now. Add edge from previous node if exists.
            if let Some(ref p) = prev {
                self.add_edge(
                    *p,
                    curr.unwrap(),
                    self.sequences.len(),
                    weights[q - 1] + weights[q],
                );
            }

            prev = curr;
        }

        if let Some((unaligned_end, _)) = nodes_unaligned_end {
            self.add_edge(
                prev.unwrap(),
                unaligned_end,
                self.sequences.len(),
                weights[*last] + weights[*last + 1],
            );
        }

        self.sequences.push(Sequence(
            sequence_name.to_owned(),
            nodes_unaligned_begin.unwrap().0,
        ));

        self.post_process();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "edge index overflow")]
    fn edge_index_overflow_panics() {
        let mut g = POAGraph::<u8>::new();
        let nodes: Vec<_> = (0..30).map(|_| g.add_node(b'A')).collect();
        for (i, &s) in nodes.iter().enumerate() {
            for &t in &nodes[i + 1..] {
                g.add_edge(s, t, 0, 1);
            }
        }
    }
}
