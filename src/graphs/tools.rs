use std::cell::RefCell;
use rustc_hash::FxHashSet;
use crate::graphs::AlignableRefGraph;

pub fn rev_postorder_nodes<G: AlignableRefGraph>(graph: &G) -> Vec<G::NodeIndex> {
    let mut ordered = Vec::with_capacity(graph.node_count_with_start_and_end());

    let mut stack = vec![
        (graph.start_node(), RefCell::new(graph.successors(graph.start_node())))
    ];

    let next_valid_child = |succ_iter: &RefCell<G::SuccessorIterator<'_>>, visited: &FxHashSet<G::NodeIndex>| {
        while let Some(child) = succ_iter.borrow_mut().next() {
            if !visited.contains(&child) {
                return Some(child)
            }
        }

        None
    };

    let mut visited = FxHashSet::default();
    while !stack.is_empty() {
        let (_, succ_iter) = stack.last().unwrap();

        if let Some(child) = next_valid_child(succ_iter, &visited) {
            visited.insert(child);
            stack.push((child, RefCell::new(graph.successors(child))))
        } else {
            let last = stack.pop().unwrap();
            ordered.push(last.0);
        }
    }

    ordered.reverse();
    ordered
}