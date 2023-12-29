use std::cmp::max;
use std::collections::VecDeque;
use rustc_hash::FxHashSet;
use crate::bubbles::finder::SuperbubbleFinder;
use crate::graphs::{AlignableRefGraph, NodeIndexType};

#[derive(Copy, Clone)]
enum BubbleNode<N> {
    /// Indicates that this node is not a bubble entrance or exit
    None,

    /// Represents a bubble entrance, the associated data is the corresponding exit node ID
    Entrance(N),

    /// Represents a bubble exit, the associated data is the corresponding entrance node ID
    Exit(N)
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


#[derive(Clone)]
pub struct BubbleIndex<N> {
    /// Vector indicating whether a node is a bubble entrance
    bubble_entrance: Vec<BubbleNode<N>>,

    /// Vector indicating whether a node is a bubble exit
    bubble_exit: Vec<BubbleNode<N>>,

    /// For each node, stores which bubbles it is part of, and the distance to the bubble exit
    node_bubble_map: Vec<Vec<NodeBubbleMap<N>>>,

    /// For each node, stores the shortest and longest path length to the POA graph end node
    dist_to_end: Vec<(usize, usize)>
}

impl<N> BubbleIndex<N>
where
    N: NodeIndexType,
{
    pub fn new<G>(graph: &G) -> Self
        where G: AlignableRefGraph<NodeIndex=N>
    {
        // Identify the bubbles in the graph and store entrances and exits
        let finder = SuperbubbleFinder::new(graph);

        // Two separate lists because nodes can be both an entrance and an exit
        let mut bubble_entrances = vec![BubbleNode::None; graph.node_count_with_start_and_end()];
        let mut bubble_exits = vec![BubbleNode::None; graph.node_count_with_start_and_end()];

        for (entrance, exit) in finder.iter() {
            bubble_entrances[entrance.index()] = BubbleNode::Entrance(exit);
            bubble_exits[exit.index()] = BubbleNode::Exit(entrance);
        }

        // To identify which nodes are part of which bubbles, and the path length to the corresponding
        // bubble exit, we run a backward BFS from the POA graph end node.
        // Since the edges in the POA graph do not have weights, we can use BFS also to compute
        // the shortest path length from the graph end node (our BFS starting point).
        // We keep a stack of bubbles that we have entered, and each time a node is popped
        // from the queue, we store which bubbles are currently active for the popped node, and the
        // shortest path distance to the exit.
        let mut node_bubble_map = vec![Vec::default(); graph.node_count_with_start_and_end()];
        let mut dists_to_end = vec![(0, 0); graph.node_count_with_start_and_end()];

        // Check if the end node is a bubble exit, and if yes, initialize the bubble stack with that bubble.
        let end_node_bubblestack = if bubble_exits[graph.end_node().index()].is_exit() {
            vec![(0usize, graph.end_node())]
        } else {
            Vec::default()
        };

        let mut queue: VecDeque<_> = vec![
            (graph.end_node(), 0usize, end_node_bubblestack)
        ].into();

        let mut visited = FxHashSet::default();
        visited.insert(graph.end_node());

        while !queue.is_empty() {
            let (curr, dist_from_end, bubble_stack) = queue.pop_front().unwrap();

            // For each active bubble (as indicated by the presence of a bubble on the bubble stack),
            // store this bubble in the node -> bubble map, accompanied with the minimum distance to the exit.
            for (bubble_dist_from_end, bubble_exit) in bubble_stack.iter() {
                node_bubble_map[curr.index()].push(NodeBubbleMap {
                    bubble_exit: *bubble_exit,
                    min_dist_to_exit: dist_from_end - *bubble_dist_from_end,
                    max_dist_to_exit: 0 // to be computed later
                })
            }

            // BFS can be used to compute the shortest path length to the POA graph end node
            dists_to_end[curr.index()].0 = dist_from_end;

            for pred in graph.predecessors(curr) {
                if !visited.contains(&pred) {
                    let new_dist_from_end = dist_from_end + 1;
                    let mut new_bubble_stack = bubble_stack.clone();

                    // if the predecessor is an entrance, remove it from the stack since we're moving backward,
                    // and thus are exiting this bubble.
                    if bubble_entrances[pred.index()].is_entrance() {
                        let (bubble_dist_from_start, bubble_exit) = new_bubble_stack.pop().unwrap();
                        node_bubble_map[pred.index()].push(NodeBubbleMap {
                            bubble_exit,
                            min_dist_to_exit: new_dist_from_end - bubble_dist_from_start,
                            max_dist_to_exit: 0 // to be computed later
                        });
                    }

                    // We're entering a new bubble
                    if bubble_exits[pred.index()].is_exit() {
                        new_bubble_stack.push((new_dist_from_end, pred));
                    }

                    visited.insert(pred);
                    queue.push_back((pred, new_dist_from_end, new_bubble_stack));
                }
            }
        }

        // While BFS can be found to compute the shortest path length, to compute the
        // longest path length to the POA graph end node we separately process each node in post order.
        for n in finder.inv_rev_postorder().iter().rev() {
            let mut max_dist_to_end = 0;

            for succ in graph.successors(*n) {
                max_dist_to_end = max(max_dist_to_end, dists_to_end[succ.index()].1 + 1);
            }

            dists_to_end[n.index()].1 = max_dist_to_end;

            // Compute max distance to bubble exits
            for bubble in node_bubble_map[n.index()].iter_mut() {
                bubble.max_dist_to_exit = max_dist_to_end - dists_to_end[bubble.bubble_exit.index()].1;
            }
        }

        Self {
            bubble_entrance: bubble_entrances,
            bubble_exit: bubble_exits,
            node_bubble_map,
            dist_to_end: dists_to_end
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
        &self.node_bubble_map[node.index()]
    }

    #[inline]
    pub fn node_is_part_of_bubble(&self, node: N) -> bool {
        !self.node_bubble_map[node.index()].is_empty()
    }

    #[inline]
    pub fn num_bubbles(&self) -> usize {
        self.bubble_entrance.iter()
            .filter(|v| v.is_entrance())
            .count()
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


#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct NodeBubbleMap<N> {
    pub bubble_exit: N,
    pub min_dist_to_exit: usize,
    pub max_dist_to_exit: usize
}

impl<N> NodeBubbleMap<N>
where
    N: NodeIndexType,
{
    pub fn new(bubble_exit: N, min_dist_to_exit: usize, max_dist_to_exit: usize) -> Self {
        NodeBubbleMap {
            bubble_exit,
            min_dist_to_exit,
            max_dist_to_exit
        }
    }
}

#[cfg(test)]
mod tests {
    use petgraph::graph::NodeIndex;
    use crate::bubbles::index::NodeBubbleMap;
    use crate::graphs::AlignableRefGraph;
    use crate::graphs::mock::{create_test_graph1, create_test_graph2};
    use super::BubbleIndex;

    type NIx = NodeIndex<crate::graphs::mock::NIx>;

    #[test]
    pub fn test_bubble_map_builder() {
        let graph1 = create_test_graph1();
        let index1 = BubbleIndex::new(&graph1);

        let truth1 = [
            vec![NodeBubbleMap::new(NIx::new(1), 1),],
            vec![NodeBubbleMap::new(NIx::new(2), 1), NodeBubbleMap::new(NIx::new(1), 0)],
            vec![NodeBubbleMap::new(NIx::new(2), 0)],
            vec![NodeBubbleMap::new(NIx::new(4), 1)],
            vec![NodeBubbleMap::new(NIx::new(5), 1), NodeBubbleMap::new(NIx::new(4), 0)],
            vec![NodeBubbleMap::new(NIx::new(5), 0)],
            vec![NodeBubbleMap::new(NIx::new(7), 1)],
            vec![NodeBubbleMap::new(NIx::new(8), 1), NodeBubbleMap::new(NIx::new(7), 0)],
            vec![NodeBubbleMap::new(NIx::new(8), 0)]
        ];

        for n in graph1.node_indices() {
            if n == graph1.end_node() {
                continue;
            }

            let bubbles_no_end_node: Vec<_> = index1.get_node_bubbles(n)
                .iter()
                .filter(|b| b.bubble_exit != graph1.end_node())
                .copied()
                .collect();
            assert_eq!(&bubbles_no_end_node, &truth1[n.index()])
        }

        let graph2 = create_test_graph2();
        let index2 = BubbleIndex::new(&graph2);

        let truth2 = [
            vec![NodeBubbleMap::new(NIx::new(2), 1)],
            vec![NodeBubbleMap::new(NIx::new(2), 1)],
            vec![NodeBubbleMap::new(NIx::new(7), 2), NodeBubbleMap::new(NIx::new(2), 0)],
            vec![NodeBubbleMap::new(NIx::new(7), 1)],
            vec![NodeBubbleMap::new(NIx::new(6), 2), NodeBubbleMap::new(NIx::new(7), 3)],
            vec![NodeBubbleMap::new(NIx::new(7), 2), NodeBubbleMap::new(NIx::new(6), 1)],
            vec![NodeBubbleMap::new(NIx::new(7), 1), NodeBubbleMap::new(NIx::new(6), 0)],
            vec![NodeBubbleMap::new(NIx::new(14), 1), NodeBubbleMap::new(NIx::new(7), 0)],
            vec![NodeBubbleMap::new(NIx::new(7), 3), NodeBubbleMap::new(NIx::new(6), 2)],
            vec![NodeBubbleMap::new(NIx::new(7), 2), NodeBubbleMap::new(NIx::new(6), 1)],
            vec![NodeBubbleMap::new(NIx::new(11), 1), NodeBubbleMap::new(NIx::new(7), 2)],
            vec![NodeBubbleMap::new(NIx::new(7), 1), NodeBubbleMap::new(NIx::new(11), 0)],
            vec![NodeBubbleMap::new(NIx::new(14), 1)],
            vec![NodeBubbleMap::new(NIx::new(14), 1)],
            vec![NodeBubbleMap::new(NIx::new(14), 0)],
        ];

        for n in graph2.node_indices() {
            assert_eq!(index2.get_node_bubbles(n), &truth2[n.index()])
        }
    }

    #[test]
    pub fn test_dist_to_exit() {
        let graph = create_test_graph2();
        let index = BubbleIndex::new(&graph);

        assert_eq!(index.get_dist_to_end(), &vec![
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
        ]);
    }

}
