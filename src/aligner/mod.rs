pub mod offsets;
pub mod extend;
pub mod scoring;
pub mod alignment;
pub mod visited;
pub mod layers;

use crate::graphs::AlignableGraph;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentCosts;

use crate::debug::DebugOutputWriter;

pub use alignment::{AlignedPair, Alignment};

use scoring::LayerCompute;


pub struct PoastaAligner<'a, C>
where
    C: AlignmentCosts
{
    costs: C,
    debug_output: Option<&'a DebugOutputWriter>,
}

impl<'a, C> PoastaAligner<'a, C>
where
    C: AlignmentCosts
{
    pub fn new(costs: C) -> Self {
        Self {
            costs,
            debug_output: None,
        }
    }

    pub fn new_with_debug_output(costs: C, debug_writer: &'a DebugOutputWriter) -> Self {
        PoastaAligner {
            costs,
            debug_output: Some(debug_writer),
        }
    }

    pub fn align<O, G, S>(&mut self, graph: &G, sequence: &S) -> (usize, Alignment<G::NodeIndex>)
    where
        O: OffsetType,
        G: AlignableGraph,
        S: AsRef<[u8]>,
    {
        let seq = sequence.as_ref();
        let max_offset = O::max_value().as_usize();

        assert!(seq.len() - 1 < max_offset, "Sequence is too long for Offset integer type!");

        let mut layers: <C as AlignmentCosts>::LayerComputeType<'_, O> = self.costs.new_layer_compute(graph);

        if let Some(debug_output) = self.debug_output {
            layers = layers.with_debug(debug_output);
        }

        let mut score = 0;
        let reached_end_state;
        loop {
            eprintln!("-------- EXTEND - Score {}", score);
            if let Some(end_state) = layers.extend(graph, seq) {
                reached_end_state = end_state;
                break;
            }

            score += 1;

            eprintln!("-------- NEXT - Score {}", score);
            if let Some(debug_output) = self.debug_output {
            }
            layers.next_layer(graph, seq, score);
        }

        let alignment = layers.backtrace(graph, seq, reached_end_state.0, &reached_end_state.1);

        // let dp_cells = ((seq.len() + 1) * graph.node_count()) * 3;
        // eprintln!("States explored: {:?}, {:.1}% of DP cells ({})", layers.num_nodes(), (layers.num_nodes() as f64 * 100.0) / dp_cells as f64, dp_cells);

        (score as usize, alignment)
    }
}
