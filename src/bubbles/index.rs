use std::collections::VecDeque;
use rustc_hash::FxHashSet;
use crate::bubbles::finder::SuperbubbleFinder;
use crate::graphs::{AlignableRefGraph, NodeIndexType};

#[derive(Copy, Clone)]
enum BubbleNode<N> {
    /// Represents a bubble entrance, the associated data is the corresponding exit node ID
    Entrance(N),

    /// Represents a bubble exit, the associated data is the corresponding entrance node ID
    Exit(N)
}


#[derive(Clone)]
pub struct BubbleIndex<N> {
    /// Vector indicating whether a node is a bubble entrance
    bubble_entrance: Vec<Option<BubbleNode<N>>>,

    /// Vector indicating whether a node is a bubble exit
    bubble_exit: Vec<Option<BubbleNode<N>>>,

    /// For each node, stores which bubbles it is part of, and the distance to the bubble exit
    node_bubble_map: Vec<Vec<NodeBubbleMap<N>>>,
}

impl<N> BubbleIndex<N>
where
    N: NodeIndexType,
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
            .filter(|v| v.is_some())
            .count()
    }
}

pub struct BubbleIndexBuilder<'a, G>
where
    G: AlignableRefGraph,
{
    graph: &'a G,

    /// Array indicating whether a node is a bubble entrance
    bubble_entrance: Vec<Option<BubbleNode<G::NodeIndex>>>,

    /// Array indicating whether a node is a bubble exit
    bubble_exit: Vec<Option<BubbleNode<G::NodeIndex>>>,

    /// A list of bubbles containing a particular node
    node_bubble_map: Vec<Vec<NodeBubbleMap<G::NodeIndex>>>,

    /// For each node, the distance to the graph end node
    node_dist_to_end: Vec<usize>,

    /// Which nodes have we already processed
    visited: FxHashSet<G::NodeIndex>,
}

impl<'a, G> BubbleIndexBuilder<'a, G>
where
    G: AlignableRefGraph,
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
            bubble_entrance: bubble_entrances,
            bubble_exit: bubble_exits,
            node_bubble_map: vec![Vec::default(); graph.node_count_with_start_and_end()],
            node_dist_to_end: vec![0; graph.node_count_with_start_and_end()],
            visited: FxHashSet::default(),
        }
    }

    fn backward_bfs(&mut self) {
        // BFS queue, containing the next nodes to explore, with per node the distance from start,
        // along with a stack of active bubbles
        let end_node_bubblestack = if self.bubble_exit[self.graph.end_node().index()].is_some() {
            vec![(0usize, self.graph.end_node())]
        } else {
            Vec::default()
        };
        let mut queue: VecDeque<_> = vec![
            (self.graph.end_node(), 0usize, end_node_bubblestack)
        ].into();

        self.visited.insert(self.graph.end_node());

        while !queue.is_empty() {
            let (curr, dist_from_end, bubble_stack) = queue.pop_front().unwrap();

            for (bubble_dist_from_end, bubble_exit) in bubble_stack.iter() {
                self.node_bubble_map[curr.index()].push(NodeBubbleMap {
                    bubble_exit: *bubble_exit,
                    dist_to_exit: dist_from_end - *bubble_dist_from_end
                })
            }

            self.node_dist_to_end[curr.index()] = dist_from_end;

            for pred in self.graph.predecessors(curr) {
                if !self.visited.contains(&pred) {
                    let new_dist_from_end = dist_from_end + 1;
                    let mut new_bubble_stack = bubble_stack.clone();

                    if self.bubble_entrance[pred.index()].is_some() {
                        let (bubble_dist_from_start, bubble_exit) = new_bubble_stack.pop().unwrap();
                        self.node_bubble_map[pred.index()].push(NodeBubbleMap {
                            bubble_exit,
                            dist_to_exit: new_dist_from_end - bubble_dist_from_start
                        });
                    }

                    if self.bubble_exit[pred.index()].is_some() {
                        new_bubble_stack.push((new_dist_from_end, pred));
                    }

                    self.visited.insert(pred);
                    queue.push_back((pred, new_dist_from_end, new_bubble_stack));
                }
            }
        }
    }

    pub fn build(mut self) -> (BubbleIndex<G::NodeIndex>, Vec<usize>) {
        self.backward_bfs();

        (BubbleIndex {
            bubble_entrance: self.bubble_entrance,
            bubble_exit: self.bubble_exit,
            node_bubble_map: self.node_bubble_map,
        }, self.node_dist_to_end)
    }
}


#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct NodeBubbleMap<N> {
    pub bubble_exit: N,
    pub dist_to_exit: usize
}

impl<N> NodeBubbleMap<N>
where
    N: NodeIndexType,
{
    pub fn new(bubble_exit: N, dist_to_exit: usize) -> Self {
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
        let index1 = BubbleIndexBuilder::new(&graph1)
            .build();

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
            assert_eq!(index1.get_node_bubbles(n), &truth1[n.index()])
        }

        let graph2 = create_test_graph2();
        let index2 = BubbleIndexBuilder::new(&graph2)
            .build();

        let truth2 = [
            vec![NodeBubbleMap::new(NIx::new(2), 1)],
            vec![NodeBubbleMap::new(NIx::new(2), 1)],
            vec![NodeBubbleMap::new(NIx::new(7), 2), NodeBubbleMap::new(NIx::new(2), 0)],
            vec![NodeBubbleMap::new(NIx::new(7), 1)],
            vec![NodeBubbleMap::new(NIx::new(6), 2), NodeBubbleMap::new(NIx::new(7), 3)],
            vec![NodeBubbleMap::new(NIx::new(7), 2), NodeBubbleMap::new(NIx::new(6), 1)],
            vec![NodeBubbleMap::new(NIx::new(7), 1), NodeBubbleMap::new(NIx::new(6), 0)],
            vec![NodeBubbleMap::new(NIx::new(13), 1), NodeBubbleMap::new(NIx::new(7), 0)],
            vec![NodeBubbleMap::new(NIx::new(7), 3), NodeBubbleMap::new(NIx::new(6), 2)],
            vec![NodeBubbleMap::new(NIx::new(7), 2), NodeBubbleMap::new(NIx::new(6), 1)],
            vec![NodeBubbleMap::new(NIx::new(11), 1), NodeBubbleMap::new(NIx::new(7), 2)],
            vec![NodeBubbleMap::new(NIx::new(7), 1), NodeBubbleMap::new(NIx::new(11), 0)],
            vec![NodeBubbleMap::new(NIx::new(13), 1)],
            vec![NodeBubbleMap::new(NIx::new(13), 0)],
            vec![NodeBubbleMap::new(NIx::new(13), 1)],
        ];

        for n in graph2.node_indices() {
            assert_eq!(index2.get_node_bubbles(n), &truth2[n.index()])
        }
    }

}
