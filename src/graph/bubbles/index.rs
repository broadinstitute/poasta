//! Bubble index: for every node, the bubbles it lies in with its shortest
//! and longest distance to each exit, plus the shortest and longest
//! distance to the graph end.
//!
//! Built in three passes over ranks, all on flat vectors. The index is
//! rebuilt for every alignment in the POA flow, because the graph changes
//! between alignments, and on graphs with many bubbles the old hash-map
//! and BFS build was 12 to 15 % of the run.

use crate::graph::bubbles::finder::SuperbubbleFinder;
use crate::graph::traits::{GraphNodeId, GraphWithNodeOrdering};

#[derive(Copy, Clone)]
enum BubbleNode<N> {
    /// Indicates that this node is not a bubble entrance or exit
    None,

    /// Represents a bubble entrance, the associated data is the corresponding exit node ID
    Entrance(N),

    /// Represents a bubble exit, the associated data is the corresponding entrance node ID
    Exit(N),
}

impl<N> BubbleNode<N> {
    #[inline]
    pub fn is_entrance(&self) -> bool {
        matches!(self, BubbleNode::Entrance(_))
    }

    #[inline]
    pub fn is_exit(&self) -> bool {
        matches!(self, BubbleNode::Exit(_))
    }
}

/// Node → bubble membership, cached alongside the (min, max) path length from
/// that node to the corresponding bubble exit.
#[derive(Clone)]
pub struct BubbleIndex<N> {
    /// Vector indicating whether a node is a bubble entrance
    bubble_entrance: Vec<BubbleNode<N>>,

    /// Vector indicating whether a node is a bubble exit
    bubble_exit: Vec<BubbleNode<N>>,

    /// CSR offsets by node index: the bubbles of node `i` are
    /// `node_bubbles[node_bubble_offsets[i]..node_bubble_offsets[i + 1]]`.
    node_bubble_offsets: Vec<u32>,

    /// Every node's bubbles, with the distances to each bubble's exit.
    node_bubbles: Vec<NodeBubbleMap<N>>,

    /// For each node, stores the shortest and longest path length to the graph end node
    dist_to_end: Vec<(usize, usize)>,
}

impl<N> BubbleIndex<N>
where
    N: GraphNodeId,
{
    pub fn new<G>(graph: &G) -> Self
    where
        G: GraphWithNodeOrdering<NodeType = N>,
    {
        // Identify the bubbles in the graph and store entrances and exits
        let finder = SuperbubbleFinder::new(graph);

        // Two separate lists because nodes can be both an entrance and an exit
        let mut bubble_entrances = vec![BubbleNode::None; graph.node_count()];
        let mut bubble_exits = vec![BubbleNode::None; graph.node_count()];

        for (entrance, exit) in finder.iter() {
            bubble_entrances[entrance.index()] = BubbleNode::Entrance(exit);
            bubble_exits[exit.index()] = BubbleNode::Exit(entrance);
        }

        // Every path from a node in a bubble to the graph end passes through the
        // bubble's exit, so the node's distance to the exit is its distance to
        // the end minus the exit's. That holds for the shortest and the longest
        // path alike, so one pass over ranks gives both and no BFS is needed.
        // Successors have higher ranks, so their entries are final when a node
        // is processed.
        let n = graph.node_count();
        let mut dist_to_end = vec![(0usize, 0usize); n];
        for rank in (0..n).rev() {
            let node = graph.rank_to_node(rank);
            let (mut lo, mut hi) = (usize::MAX, 0usize);
            for succ in graph.successors(node) {
                let (succ_lo, succ_hi) = dist_to_end[succ.index()];
                lo = lo.min(succ_lo + 1);
                hi = hi.max(succ_hi + 1);
            }
            if lo == usize::MAX {
                lo = 0;
            }
            dist_to_end[node.index()] = (lo, hi);
        }

        let mut node_bubble_offsets = vec![0u32; n + 1];
        for_each_node_bubble(graph, &bubble_entrances, &bubble_exits, |node, _| {
            node_bubble_offsets[node.index() + 1] += 1;
        });
        for i in 0..n {
            node_bubble_offsets[i + 1] += node_bubble_offsets[i];
        }

        let total = node_bubble_offsets[n] as usize;
        let placeholder = NodeBubbleMap::new(graph.end_node(), 0, 0);
        let mut node_bubbles = vec![placeholder; total];
        let mut cursor: Vec<u32> = node_bubble_offsets.clone();
        for_each_node_bubble(graph, &bubble_entrances, &bubble_exits, |node, exit| {
            let slot = &mut cursor[node.index()];
            node_bubbles[*slot as usize] = NodeBubbleMap {
                bubble_exit: exit,
                min_dist_to_exit: (dist_to_end[node.index()].0 - dist_to_end[exit.index()].0) as u32,
                max_dist_to_exit: (dist_to_end[node.index()].1 - dist_to_end[exit.index()].1) as u32,
            };
            *slot += 1;
        });

        Self {
            bubble_entrance: bubble_entrances,
            bubble_exit: bubble_exits,
            node_bubble_offsets,
            node_bubbles,
            dist_to_end,
        }
    }

    #[inline]
    pub fn is_entrance(&self, node: N) -> bool {
        self.bubble_entrance[node.index()].is_entrance()
    }

    #[inline]
    pub fn is_exit(&self, node: N) -> bool {
        self.bubble_exit[node.index()].is_exit()
    }

    #[inline]
    pub fn get_node_bubbles(&self, node: N) -> &[NodeBubbleMap<N>] {
        let i = node.index();
        &self.node_bubbles[self.node_bubble_offsets[i] as usize..self.node_bubble_offsets[i + 1] as usize]
    }

    #[inline]
    pub fn node_is_part_of_bubble(&self, node: N) -> bool {
        let i = node.index();
        self.node_bubble_offsets[i] != self.node_bubble_offsets[i + 1]
    }

    #[inline]
    pub fn num_bubbles(&self) -> usize {
        self.bubble_entrance.iter().filter(|v| v.is_entrance()).count()
    }

    #[inline]
    pub fn get_min_dist_to_end(&self, node: N) -> usize {
        self.dist_to_end[node.index()].0
    }

    #[inline]
    pub fn get_max_dist_to_end(&self, node: N) -> usize {
        self.dist_to_end[node.index()].1
    }

    #[inline]
    pub fn get_dist_to_end(&self) -> &[(usize, usize)] {
        &self.dist_to_end
    }
}

/// Walk the ranks once with a stack of the bubbles that are open, and
/// report every (node, exit) pair.
///
/// This works because a superbubble of a DAG is an interval in a
/// DFS-based topological order (Gärtner et al. 2018, the property the
/// finder itself uses): a node lies in a bubble exactly when its rank is
/// between the entrance's and the exit's. Bubbles nest or are disjoint,
/// so the innermost open bubble is always on top of the stack and the
/// exit that closes at this rank is the top entry.
///
/// A node that is the exit of one bubble and the entrance of the next
/// belongs to both, with distance zero to the first exit. Entries come
/// out outer to inner with the node's own entrance bubble last. Pruning
/// takes `all()` over them, so the order does not matter to it.
fn for_each_node_bubble<N, G>(
    graph: &G,
    entrances: &[BubbleNode<N>],
    exits: &[BubbleNode<N>],
    mut emit: impl FnMut(N, N),
) where
    N: GraphNodeId,
    G: GraphWithNodeOrdering<NodeType = N>,
{
    let mut open: Vec<N> = Vec::new();
    for rank in 0..graph.node_count() {
        let node = graph.rank_to_node(rank);
        for &exit in &open {
            emit(node, exit);
        }
        if exits[node.index()].is_exit() {
            let closed = open.pop();
            debug_assert_eq!(closed, Some(node));
        }
        if let BubbleNode::Entrance(exit) = entrances[node.index()] {
            open.push(exit);
            emit(node, exit);
        }
    }
    debug_assert!(open.is_empty());
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct NodeBubbleMap<N> {
    pub bubble_exit: N,
    pub min_dist_to_exit: u32,
    pub max_dist_to_exit: u32,
}

impl<N> NodeBubbleMap<N>
where
    N: GraphNodeId,
{
    pub fn new(bubble_exit: N, min_dist_to_exit: u32, max_dist_to_exit: u32) -> Self {
        NodeBubbleMap {
            bubble_exit,
            min_dist_to_exit,
            max_dist_to_exit,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::BubbleIndex;
    use crate::graph::bubbles::finder::SuperbubbleFinder;
    use crate::graph::mock::{create_test_graph1, create_test_graph2, NIx};
    use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};
    use petgraph::graph::NodeIndex;

    type Ix = NodeIndex<NIx>;

    /// Reduce a bubble map entry to (exit, min_dist_to_exit) for comparison
    /// against hand-checked truth data, sorted by exit so the per-node order
    /// of the index does not matter.
    fn reduced(index: &BubbleIndex<Ix>, node: Ix, skip_end: Option<Ix>) -> Vec<(Ix, u32)> {
        let mut entries: Vec<_> = index
            .get_node_bubbles(node)
            .iter()
            .filter(|b| Some(b.bubble_exit) != skip_end)
            .map(|b| (b.bubble_exit, b.min_dist_to_exit))
            .collect();
        entries.sort_by_key(|(exit, _)| exit.index());
        entries
    }

    fn sorted(truth: &[(Ix, u32)]) -> Vec<(Ix, u32)> {
        let mut truth = truth.to_vec();
        truth.sort_by_key(|(exit, _)| exit.index());
        truth
    }

    #[test]
    pub fn test_bubble_map_builder() {
        let graph1 = create_test_graph1();
        let index1 = BubbleIndex::new(&graph1);

        let truth1 = [
            vec![(Ix::new(1), 1)],
            vec![(Ix::new(2), 1), (Ix::new(1), 0)],
            vec![(Ix::new(2), 0)],
            vec![(Ix::new(4), 1)],
            vec![(Ix::new(5), 1), (Ix::new(4), 0)],
            vec![(Ix::new(5), 0)],
            vec![(Ix::new(7), 1)],
            vec![(Ix::new(8), 1), (Ix::new(7), 0)],
            vec![(Ix::new(8), 0)],
        ];

        for n in graph1.all_nodes_iter() {
            if n == graph1.end_node() {
                continue;
            }
            assert_eq!(reduced(&index1, n, Some(graph1.end_node())), sorted(&truth1[n.index()]));
        }

        let graph2 = create_test_graph2();
        let index2 = BubbleIndex::new(&graph2);

        let truth2 = [
            vec![(Ix::new(2), 1)],
            vec![(Ix::new(2), 1)],
            vec![(Ix::new(7), 2), (Ix::new(2), 0)],
            vec![(Ix::new(7), 1)],
            vec![(Ix::new(6), 2), (Ix::new(7), 3)],
            vec![(Ix::new(7), 2), (Ix::new(6), 1)],
            vec![(Ix::new(7), 1), (Ix::new(6), 0)],
            vec![(Ix::new(14), 1), (Ix::new(7), 0)],
            vec![(Ix::new(7), 3), (Ix::new(6), 2)],
            vec![(Ix::new(7), 2), (Ix::new(6), 1)],
            vec![(Ix::new(11), 1), (Ix::new(7), 2)],
            vec![(Ix::new(7), 1), (Ix::new(11), 0)],
            vec![(Ix::new(14), 1)],
            vec![(Ix::new(14), 1)],
            vec![(Ix::new(14), 0)],
        ];

        for n in graph2.all_nodes_iter() {
            assert_eq!(reduced(&index2, n, None), sorted(&truth2[n.index()]));
        }
    }

    /// `max_dist_to_exit` is computed from `dist_to_end`, which is itself
    /// checked against hand-verified numbers in `test_dist_to_exit`. Confirm
    /// the two stay consistent: `max_dist_to_exit == max_dist_to_end(node) -
    /// max_dist_to_end(bubble_exit)`.
    #[test]
    pub fn test_max_dist_to_exit_consistent_with_dist_to_end() {
        let graph = create_test_graph2();
        let index = BubbleIndex::new(&graph);

        for n in graph.all_nodes_iter() {
            for bubble in index.get_node_bubbles(n) {
                let expected = (index.get_max_dist_to_end(n) - index.get_max_dist_to_end(bubble.bubble_exit)) as u32;
                assert_eq!(bubble.max_dist_to_exit, expected);
            }
        }
    }

    /// The twin of the test above for the shortest distance.
    #[test]
    pub fn test_min_dist_to_exit_consistent_with_dist_to_end() {
        let graph = create_test_graph2();
        let index = BubbleIndex::new(&graph);

        for n in graph.all_nodes_iter() {
            for bubble in index.get_node_bubbles(n) {
                let expected = (index.get_min_dist_to_end(n) - index.get_min_dist_to_end(bubble.bubble_exit)) as u32;
                assert_eq!(bubble.min_dist_to_exit, expected);
            }
        }
    }

    /// A node lists a bubble exactly when its rank lies between the
    /// bubble's entrance and exit ranks.
    #[test]
    pub fn test_bubble_membership_is_an_interval() {
        let graph = create_test_graph2();
        let index = BubbleIndex::new(&graph);

        for (s, t) in SuperbubbleFinder::new(&graph).iter() {
            let ranks = graph.node_rank(s)..=graph.node_rank(t);
            for n in graph.all_nodes_iter() {
                let listed = index.get_node_bubbles(n).iter().any(|b| b.bubble_exit == t);
                assert_eq!(listed, ranks.contains(&graph.node_rank(n)), "node {n:?}, bubble ({s:?}, {t:?})");
            }
        }
    }

    #[test]
    pub fn test_dist_to_exit() {
        let graph = create_test_graph2();
        let index = BubbleIndex::new(&graph);

        assert_eq!(
            index.get_dist_to_end(),
            &vec![
                (4, 10),
                (4, 9),
                (3, 8),
                (2, 4),
                (4, 7),
                (3, 6),
                (2, 4),
                (1, 3),
                (4, 6),
                (3, 5),
                (3, 5),
                (2, 4),
                (1, 2),
                (1, 1),
                (0, 0)
            ]
        );
    }
}
