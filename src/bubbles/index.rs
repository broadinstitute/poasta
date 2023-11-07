use std::collections::VecDeque;
use rustc_hash::FxHashSet;
use crate::aligner::offsets::OffsetType;
use crate::bubbles::finder::SuperbubbleFinder;
use crate::graphs::{AlignableGraph, NodeIndexType};

#[derive(Copy, Clone)]
enum BubbleNode<N> {
    /// Represents a bubble entrance, the associated data is the corresponding exit node ID
    Entrance(N),

    /// Represents a bubble exit, the associated data is the corresponding entrance node ID
    Exit(N)
}


pub struct BubbleIndex<N, O> {
    /// Vector indicating whether a node is a bubble entrance
    bubble_entrance: Vec<Option<BubbleNode<N>>>,

    /// Vector indicating whether a node is a bubble exit
    bubble_exit: Vec<Option<BubbleNode<N>>>,

    /// For each node, stores which bubbles it is part of, and the distance to the bubble exit
    node_bubble_map: Vec<Vec<NodeBubbleMap<N, O>>>,
}

impl<N, O> BubbleIndex<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    #[inline]
    pub fn is_entrance(&self, node: N) -> bool {
        self.bubble_entrance[node.index()].is_some()
    }

    #[inline]
    pub fn is_exit(&self, node: N) -> bool {
        self.bubble_exit[node.index()].is_some()
    }

    #[inline]
    pub fn get_node_bubbles(&self, node: N) -> &[NodeBubbleMap<N, O>] {
        &self.node_bubble_map[node.index()]
    }

    pub fn node_is_part_of_bubble(&self, node: N) -> bool {
        !self.node_bubble_map[node.index()].is_empty()
    }

    #[inline]
    pub fn num_bubbles(&self) -> usize {
        self.bubble_entrance.iter()
            .filter(|v| v.is_some())
            .count()
    }
}

pub struct BubbleIndexBuilder<'a, O, G>
where
    G: AlignableGraph,
{
    graph: &'a G,

    /// Superbubble finder object
    finder: SuperbubbleFinder<'a, G>,

    /// Array indicating whether a node is a bubble entrance
    bubble_entrance: Vec<Option<BubbleNode<G::NodeIndex>>>,

    /// Array indicating whether a node is a bubble exit
    bubble_exit: Vec<Option<BubbleNode<G::NodeIndex>>>,

    /// A list of bubbles containing a particular node
    node_bubble_map: Vec<Vec<NodeBubbleMap<G::NodeIndex, O>>>,

    /// Which nodes have we already processed
    visited: FxHashSet<G::NodeIndex>,
}

impl<'a, O, G> BubbleIndexBuilder<'a, O, G>
where
    O: OffsetType,
    G: AlignableGraph,
{
    pub fn new(graph: &'a G) -> Self {
        let finder = SuperbubbleFinder::new(graph);

        // Two separate lists because nodes can be both an entrance and an exit
        let mut bubble_entrances = vec![None; graph.node_count_with_start_and_end()];
        let mut bubble_exits = vec![None; graph.node_count_with_start_and_end()];

        for (entrance, exit) in finder.iter() {
            bubble_entrances[entrance.index()] = Some(BubbleNode::Entrance(exit));
            bubble_exits[exit.index()] = Some(BubbleNode::Exit(entrance));
        }

        Self {
            graph,
            finder,
            bubble_entrance: bubble_entrances,
            bubble_exit: bubble_exits,
            node_bubble_map: vec![Vec::default(); graph.node_count_with_start_and_end()],
            visited: FxHashSet::default(),
        }
    }

    fn bubble_backward_bfs(&mut self, entrance: G::NodeIndex, exit: G::NodeIndex) {
        // BFS queue, containing the next nodes to explore, with per node the distance from start,
        // along with a stack of active bubbles
        let mut queue: VecDeque<_> = vec![
            (exit, 0usize, vec![(0usize, exit)])
        ].into();
        self.visited.insert(exit);

        while !queue.is_empty() {
            let (curr, dist_from_start, bubble_stack) = queue.pop_front().unwrap();

            for (bubble_dist_from_start, bubble_exit) in bubble_stack.iter() {
                self.node_bubble_map[curr.index()].push(NodeBubbleMap {
                    bubble_exit: *bubble_exit,
                    dist_to_exit: O::new(dist_from_start - *bubble_dist_from_start)
                })
            }

            if curr == entrance && self.bubble_exit[curr.index()].is_none() {
                continue;
            }

            for pred in self.graph.predecessors(curr) {
                if !self.visited.contains(&pred) {
                    let new_dist_from_start = dist_from_start + 1;
                    let mut new_bubble_stack = bubble_stack.clone();

                    if self.bubble_entrance[pred.index()].is_some() {
                        let (bubble_dist_from_start, bubble_exit) = new_bubble_stack.pop().unwrap();
                        self.node_bubble_map[pred.index()].push(NodeBubbleMap {
                            bubble_exit,
                            dist_to_exit: O::new(new_dist_from_start - bubble_dist_from_start)
                        });
                    }

                    if self.bubble_exit[pred.index()].is_some() {
                        new_bubble_stack.push((new_dist_from_start, pred));
                    }

                    self.visited.insert(pred);
                    queue.push_back((pred, new_dist_from_start, new_bubble_stack));
                }
            }
        }
    }

    pub fn build(mut self) -> BubbleIndex<G::NodeIndex, O> {
        for rpo in (0..self.graph.node_count_with_start_and_end()).rev() {
            let inv_rev_postorder = self.finder.get_inv_rev_postorder();
            let node_id = inv_rev_postorder[rpo];
            if self.visited.contains(&node_id) {
                continue;
            }

            if let Some(bubble_exit_node) = self.bubble_exit[node_id.index()] {
                let BubbleNode::Exit(entrance) = bubble_exit_node else {
                    panic!("Unexpected value for bubble exit!");
                };

                self.bubble_backward_bfs(entrance, node_id);
            }
        }

        BubbleIndex {
            bubble_entrance: self.bubble_entrance,
            bubble_exit: self.bubble_exit,
            node_bubble_map: self.node_bubble_map
        }
    }
}


#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct NodeBubbleMap<N, O> {
    pub bubble_exit: N,
    pub dist_to_exit: O
}

impl<N, O> NodeBubbleMap<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new(bubble_exit: N, dist_to_exit: O) -> Self {
        NodeBubbleMap {
            bubble_exit,
            dist_to_exit
        }
    }
}

#[cfg(test)]
mod tests {
    use petgraph::graph::NodeIndex;
    use crate::bubbles::index::NodeBubbleMap;
    use crate::graphs::mock::{create_test_graph1, create_test_graph2};
    use super::BubbleIndexBuilder;

    type NIx = NodeIndex<crate::graphs::mock::NIx>;

    #[test]
    pub fn test_bubble_map_builder() {
        let graph1 = create_test_graph1();
        let index1 = BubbleIndexBuilder::<u32, _>::new(&graph1)
            .build();

        let truth1 = [
            vec![NodeBubbleMap::new(NIx::new(1), 1u32),],
            vec![NodeBubbleMap::new(NIx::new(2), 1u32), NodeBubbleMap::new(NIx::new(1), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(2), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(4), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(5), 1u32), NodeBubbleMap::new(NIx::new(4), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(5), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(8), 1u32), NodeBubbleMap::new(NIx::new(7), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(8), 0u32)]
        ];

        for n in graph1.node_indices() {
            assert_eq!(index1.get_node_bubbles(n), &truth1[n.index()])
        }

        let graph2 = create_test_graph2();
        let index2 = BubbleIndexBuilder::<u32, _>::new(&graph2)
            .build();

        let truth2 = [
            vec![NodeBubbleMap::new(NIx::new(2), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(2), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 2u32), NodeBubbleMap::new(NIx::new(2), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(6), 2u32), NodeBubbleMap::new(NIx::new(7), 3u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 2u32), NodeBubbleMap::new(NIx::new(6), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 1u32), NodeBubbleMap::new(NIx::new(6), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(13), 1u32), NodeBubbleMap::new(NIx::new(7), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 3u32), NodeBubbleMap::new(NIx::new(6), 2u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 2u32), NodeBubbleMap::new(NIx::new(6), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(11), 1u32), NodeBubbleMap::new(NIx::new(7), 2u32)],
            vec![NodeBubbleMap::new(NIx::new(7), 1u32), NodeBubbleMap::new(NIx::new(11), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(13), 1u32)],
            vec![NodeBubbleMap::new(NIx::new(13), 0u32)],
            vec![NodeBubbleMap::new(NIx::new(13), 1u32)],
        ];

        for n in graph2.node_indices() {
            assert_eq!(index2.get_node_bubbles(n), &truth2[n.index()])
        }
    }

}
