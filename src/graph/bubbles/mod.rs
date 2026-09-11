//! Superbubbles: [`finder`] detects them, [`index`] records per node which
//! bubbles it lies in and how far the graph end is. The A* heuristic and
//! bubble-exit pruning both read the index.

pub mod finder;
pub mod index;
