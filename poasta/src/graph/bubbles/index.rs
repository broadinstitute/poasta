use crate::graph::traits::{GraphNodeId, GraphWithNodeLengths, GraphWithStartEnd};
use std::collections::BTreeMap;

use super::finder::SuperbubbleFinder;

#[derive(Copy, Clone, Debug)]
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

#[derive(Clone)]
pub struct BubbleIndex<N> {
    /// Vector indicating whether a node is a bubble entrance
    bubble_entrance: Vec<BubbleNode<N>>,

    /// Vector indicating whether a node is a bubble exit
    bubble_exit: Vec<BubbleNode<N>>,

    /// For each node, stores which bubbles it is part of, and the distance to the bubble exit
    node_bubble_map: Vec<BTreeMap<N, BubbleDistToExit>>,

    /// For each node, stores the shortest and longest path length to the POA graph end node
    dist_to_end: Vec<(usize, usize)>,
}

impl<N> BubbleIndex<N>
where
    N: GraphNodeId,
{
    pub fn new<G>(graph: &G) -> Self
    where
        G: GraphWithStartEnd<NodeType = N> + GraphWithNodeLengths<NodeType = N>,
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

        // To identify which nodes are part of which bubbles, we run BFS on the graph
        // and track entering and exiting bubbles while visiting nodes.
        let mut node_bubble_map: Vec<BTreeMap<N, BubbleDistToExit>> =
            vec![BTreeMap::default(); graph.node_count()];

        // First, iterate over nodes in reverse post order, determining which nodes
        // are contained in which bubbles.
        for n in finder.inv_rev_postorder().iter() {
            let mut to_add = vec![];

            for pred in graph.predecessors(*n) {
                // Find bubbles that contain the predecessor; those bubbles also contain the current node
                // unless the corresponding bubble exit is the current node itself.
                for pred_bubble_exit in node_bubble_map[pred.index()].keys() {
                    if *pred_bubble_exit == *n {
                        continue;
                    }

                    to_add.push(*pred_bubble_exit);
                }

                // Add the bubble and its corresponding exit to the current node if it's not already present.
                // Initialize with default values since we will compute the min and max distance to the exit
                // later.
                for bubble_to_add in &to_add {
                    node_bubble_map[n.index()]
                        .entry(*bubble_to_add)
                        .or_default();
                }

                to_add.clear();
            }

            // Check if the current node is a bubble entrance, and if so,
            // add the bubble and its corresponding exit to the current node.
            if let BubbleNode::Entrance(bubble_exit) = bubble_entrances[n.index()] {
                node_bubble_map[n.index()].entry(bubble_exit).or_default();
            }
        }

        // With the node => bubble map completed, compute the min and max distances to the bubble exits and the graph end node.
        // We compute this by visiting nodes in postorder, recursively computing the distances.
        let mut dist_to_end = vec![(0usize, 0usize); graph.node_count()];
        for n in finder.inv_rev_postorder().iter().rev() {
            let node_length = graph.node_length(*n);

            let min_dist = graph
                .successors(*n)
                .map(|succ| dist_to_end[succ.index()].0 + node_length)
                .min()
                .unwrap_or(0);

            let max_dist = graph
                .successors(*n)
                .map(|succ| dist_to_end[succ.index()].1 + node_length)
                .max()
                .unwrap_or(0);

            dist_to_end[n.index()] = (min_dist, max_dist);

            // Compute min/max distance to bubble exits
            for (bubble_exit, bubble_dists) in node_bubble_map[n.index()].iter_mut() {
                bubble_dists.min_dist =
                    dist_to_end[n.index()].0 - dist_to_end[bubble_exit.index()].0;
                bubble_dists.max_dist =
                    dist_to_end[n.index()].1 - dist_to_end[bubble_exit.index()].1;
            }
        }

        Self {
            bubble_entrance: bubble_entrances,
            bubble_exit: bubble_exits,
            node_bubble_map,
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
    pub fn get_node_bubbles(&self, node: N) -> &BTreeMap<N, BubbleDistToExit> {
        &self.node_bubble_map[node.index()]
    }

    #[inline]
    pub fn node_is_part_of_bubble(&self, node: N) -> bool {
        !self.node_bubble_map[node.index()].is_empty()
    }

    #[inline]
    pub fn num_bubbles(&self) -> usize {
        self.bubble_entrance
            .iter()
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

#[derive(Default, Copy, Clone, Debug, PartialEq, Eq)]
pub struct BubbleDistToExit {
    pub min_dist: usize,
    pub max_dist: usize,
}

impl BubbleDistToExit {
    pub fn new(min_dist: usize, max_dist: usize) -> Self {
        BubbleDistToExit { min_dist, max_dist }
    }
}

#[cfg(test)]
mod tests {
    use petgraph::graph::NodeIndex;

    use super::BubbleDistToExit;
    use super::BubbleIndex;
    use crate::graph::mock::{create_test_graph1, create_test_graph2};
    use crate::graph::traits::GraphWithStartEnd;

    type NIx = NodeIndex<crate::graph::mock::NIx>;

    #[test]
    pub fn test_bubble_map_builder() {
        let graph1 = create_test_graph1();
        let index1 = BubbleIndex::new(&graph1);

        let truth1 = [
            vec![(NIx::new(1), BubbleDistToExit::new(1, 1))],
            vec![(NIx::new(2), BubbleDistToExit::new(1, 1))],
            vec![],
            vec![(NIx::new(4), BubbleDistToExit::new(1, 1))],
            vec![(NIx::new(5), BubbleDistToExit::new(1, 1))],
            vec![],
            vec![(NIx::new(7), BubbleDistToExit::new(1, 1))],
            vec![(NIx::new(8), BubbleDistToExit::new(1, 1))],
            vec![],
            vec![],
        ];

        assert_eq!(index1.node_bubble_map.len(), truth1.len());

        for n in graph1.node_indices() {
            let i = n.index();
            let excl_end_node_bubbles = index1.node_bubble_map[i]
                .iter()
                .map(|(k, v)| (*k, *v))
                .filter(|(k, _)| *k != graph1.end_node())
                .collect::<Vec<_>>();

            assert_eq!(excl_end_node_bubbles, truth1[i]);
        }

        let graph2 = create_test_graph2();
        let index2 = BubbleIndex::new(&graph2);

        let truth2 = [
            vec![(
                NIx::new(2),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 2,
                },
            )], // node 0
            vec![(
                NIx::new(2),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 1,
                },
            )], // node 1
            vec![(
                NIx::new(7),
                BubbleDistToExit {
                    min_dist: 2,
                    max_dist: 5,
                },
            )], // node 2
            vec![(
                NIx::new(7),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 1,
                },
            )], // node 3
            vec![
                (
                    NIx::new(6),
                    BubbleDistToExit {
                        min_dist: 2,
                        max_dist: 3,
                    },
                ),
                (
                    NIx::new(7),
                    BubbleDistToExit {
                        min_dist: 3,
                        max_dist: 4,
                    },
                ),
            ], // node 4
            vec![
                (
                    NIx::new(6),
                    BubbleDistToExit {
                        min_dist: 1,
                        max_dist: 2,
                    },
                ),
                (
                    NIx::new(7),
                    BubbleDistToExit {
                        min_dist: 2,
                        max_dist: 3,
                    },
                ),
            ], // node 5
            vec![(
                NIx::new(7),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 1,
                },
            )], // node 6
            vec![(
                NIx::new(14),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 3,
                },
            )], // node 7
            vec![
                (
                    NIx::new(6),
                    BubbleDistToExit {
                        min_dist: 2,
                        max_dist: 2,
                    },
                ),
                (
                    NIx::new(7),
                    BubbleDistToExit {
                        min_dist: 3,
                        max_dist: 3,
                    },
                ),
            ], // node 8
            vec![
                (
                    NIx::new(6),
                    BubbleDistToExit {
                        min_dist: 1,
                        max_dist: 1,
                    },
                ),
                (
                    NIx::new(7),
                    BubbleDistToExit {
                        min_dist: 2,
                        max_dist: 2,
                    },
                ),
            ], // node 9
            vec![
                (
                    NIx::new(7),
                    BubbleDistToExit {
                        min_dist: 2,
                        max_dist: 2,
                    },
                ),
                (
                    NIx::new(11),
                    BubbleDistToExit {
                        min_dist: 1,
                        max_dist: 1,
                    },
                ),
            ], // node 10
            vec![(
                NIx::new(7),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 1,
                },
            )], // node 11
            vec![(
                NIx::new(14),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 2,
                },
            )], // node 12
            vec![(
                NIx::new(14),
                BubbleDistToExit {
                    min_dist: 1,
                    max_dist: 1,
                },
            )], // node 13
            vec![], // node 14
        ];

        for n in graph2.node_indices() {
            let i = n.index();

            let bubble_dists = index2.node_bubble_map[i]
                .iter()
                .map(|(k, v)| (*k, *v))
                .collect::<Vec<_>>();

            assert_eq!(bubble_dists, truth2[i]);
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
