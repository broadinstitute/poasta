//! Small graph helpers shared across modules (bubble finder, A* engine).

use crate::graph::traits::GraphWithNodeOrdering;

/// Return all nodes of `graph` in reverse-postorder (topological) order.
///
/// `GraphWithNodeOrdering` implementations already maintain this order via
/// `node_rank`/`rank_to_node` (computed once via DFS in `post_process`), so
/// this is just a materialization of `rank_to_node(0..node_count())` — no
/// second traversal.
pub fn rev_postorder_nodes<G: GraphWithNodeOrdering>(graph: &G) -> Vec<G::NodeType> {
    (0..graph.node_count()).map(|r| graph.rank_to_node(r)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::mock::create_test_graph1;
    use crate::graph::traits::GraphWithNodeOrdering;

    #[test]
    fn rev_postorder_nodes_matches_node_rank() {
        let g = create_test_graph1();
        let order = rev_postorder_nodes(&g);

        assert_eq!(order.len(), g.node_count());
        for (rank, &node) in order.iter().enumerate() {
            assert_eq!(g.node_rank(node), rank);
        }

        // Start node must come first, end node last, in a rank-ordered DAG.
        assert_eq!(order[0], g.start_node());
        assert_eq!(*order.last().unwrap(), g.end_node());
    }
}
