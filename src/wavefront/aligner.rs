use crate::graph::POAGraph;
use crate::alignment::Alignment;
use crate::debug::DebugOutputWriter;

use super::extend::ExtendPaths;
use super::compute::WFCompute;

use num::Bounded;
use crate::debug::messages::DebugOutputMessage;


pub struct WavefrontPOAligner<'a, Compute>
    where
        Compute: WFCompute,
{
    compute: Compute,
    debug_output: Option<&'a DebugOutputWriter>,
    num_seq_aligned: usize,
}

impl<'a, Compute> WavefrontPOAligner<'a, Compute>
    where
        Compute: WFCompute,
{
    pub fn new() -> Self {
        Self {
            compute: Compute::default(),
            debug_output: None,
            num_seq_aligned: 1
        }
    }

    pub fn new_with_debug_output(debug_writer: &'a DebugOutputWriter) -> Self {
        WavefrontPOAligner {
            compute: Compute::default(),
            debug_output: Some(debug_writer),
            num_seq_aligned: 1
        }
    }

    pub fn align<T: AsRef<[u8]>>(&mut self, graph: &POAGraph, sequence: T) -> Alignment {
        let seq = sequence.as_ref();

        let max_offset = match Compute::OffsetType::max_value().try_into() {
            Ok(v) => v,
            Err(_) => panic!("Could not determine maximum value for Offset type!")
        };

        assert!(seq.len() - 1 < max_offset, "Sequence is too long for Offset type!");

        self.compute.reset();
        let reached_end_node;
        let mut score: i64 = 0;
        loop {
            eprintln!("----- EXTEND: score {}", score);
            if let Some(end_node) = self.wf_extend(graph, seq) {
                eprintln!("Reached end node {:?}!", end_node);
                reached_end_node = end_node;
                break;
            }

            score += 1;
            eprintln!("----- NEXT: score {}", score);
            self.compute.next(graph, seq.len(), score);

            if let Some(debug) = self.debug_output {
                self.compute.log_debug_data(debug);
            }
        }

        eprintln!("Sequence length: {:?}", seq.len());
        self.num_seq_aligned += 1;

        self.compute.backtrace(graph, seq, reached_end_node)
    }

    fn wf_extend(&mut self, graph: &POAGraph, seq: &[u8]) -> Option<usize> {
        // Check if we can extend without mismatches along some paths in the graph. Traverse the
        // graph in a depth first manner, and update furthest reaching points accordingly.
        let candidates = self.compute.extend_candidates();
        eprintln!("- Extend start points: {:?}", candidates);

        for (start_node, offset) in candidates {
            let mut extend_paths = ExtendPaths::new(graph, seq, start_node, offset);

            while let Some(path) = extend_paths.next(&self.compute) {
                if let Some(debug) = self.debug_output {
                    debug.log(DebugOutputMessage::new_extended_path(start_node, offset, path.clone()));
                }

                self.compute.update_extended_path(start_node, path);
            }

            if let Some(end_node) = self.compute.reached_end(graph, seq.len()) {
                return Some(end_node);
            }
        }

        None
    }
}

impl<'a, Compute: WFCompute> Default for WavefrontPOAligner<'a, Compute> {
    fn default() -> Self {
        Self::new()
    }
}
