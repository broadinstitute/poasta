use std::fs::File;
use std::io::BufWriter;
use crate::graph::POAGraph;
use crate::alignment::Alignment;

use super::compute::WFCompute;
use super::fr_points::ExtendCandidate;
use super::{DiagIx, FRPoint, WFPoint};

use std::path::PathBuf;

use num::{One, Bounded};


pub struct WavefrontPOAligner<Compute>
    where
        Compute: WFCompute,
{
    compute: Compute,
    debug_output_dir: Option<PathBuf>,
    num_seq_aligned: usize,
}

impl<Compute> WavefrontPOAligner<Compute>
    where
        Compute: WFCompute,
{
    pub fn new() -> Self {
        Self {
            compute: Compute::default(),
            debug_output_dir: None,
            num_seq_aligned: 1
        }
    }

    pub fn new_with_debug_output(debug_output_dir: impl Into<PathBuf>) -> Self {
        WavefrontPOAligner {
            compute: Compute::default(),
            debug_output_dir: Some(debug_output_dir.into()),
            num_seq_aligned: 1
        }
    }

    pub fn align<T: AsRef<[u8]>>(&mut self, graph: &POAGraph, sequence: T) -> Alignment {
        let seq = sequence.as_ref();

        let max_offset = match Compute::OffsetType::max_value().try_into() {
            Ok(v) => v,
            Err(_) => panic!("Could not determine maximum value for Offset type!")
        };

        assert!(seq.len() - 1 < max_offset as usize, "Sequence is too long for Offset type!");

        self.compute.reset();
        let reached_end;
        let mut score: i64 = 0;
        loop {
            eprintln!("----- EXTEND: score {}", score);
            self.wf_extend(graph, seq);

            if let Some(debug_output_dir) = &self.debug_output_dir {
                let f_after_extend = debug_output_dir.join(format!("seq{}.g{}.s{}.score{score}.after_extend.tsv",
                                                                    self.num_seq_aligned, graph.max_rank(), seq.len()));
                let fout = File::create(f_after_extend).unwrap();
                let mut writer = BufWriter::new(fout);
                self.compute.write_csv(&mut writer).unwrap();
            }

            if let Some(end_point) = self.compute.reached_end(graph, seq.len()) {
                eprintln!("Reached end point {:?}!", end_point);
                reached_end = Some(end_point);
                break;
            }

            score += 1;
            eprintln!("----- NEXT: score {}", score);
            self.compute.next(graph, seq.len(), score);

            if let Some(debug_output_dir) = &self.debug_output_dir {
                let f_before_extend = debug_output_dir.join(format!("seq{}.g{}.s{}.score{score}.before_extend.tsv", self.num_seq_aligned, graph.max_rank(), seq.len()));
                let fout = File::create(f_before_extend).unwrap();
                let mut writer = BufWriter::new(fout);
                self.compute.write_csv(&mut writer).unwrap();
            }
        }

        eprintln!("Sequence length: {:?}", seq.len());
        self.num_seq_aligned += 1;

        self.compute.backtrace(graph, reached_end.unwrap())
    }

    fn neighboring_fr_points(&self, graph: &POAGraph, point: FRPoint<Compute::OffsetType>) -> Vec<FRPoint<Compute::OffsetType>> {
        graph.neighbors_for_rank(point.rank())
            .map(|succ| {
                let succ_rank = graph.get_node_rank(succ);
                let succ_k = point.diag() - (succ_rank - point.rank() - 1) as DiagIx;

                (succ_k, point.offset() + Compute::OffsetType::one())
            })
            .collect()
    }

    fn wf_extend(&mut self, graph: &POAGraph, seq: &[u8]) {
        // Check if we can extend without mismatches along our current diagonals. Follow paths
        // in the graph in a depth first manner, and update furthest reaching points accordingly.
        let mut stack: Vec<ExtendCandidate<Compute::OffsetType>> = self.compute.extend_candidates().into_iter()
            .flat_map(|p| {
                self.neighboring_fr_points(graph, p).into_iter()
                    .map(|neighbor| ExtendCandidate(neighbor, Some(p)))
                    .collect::<Vec<ExtendCandidate<Compute::OffsetType>>>()
            })
            .collect();
        eprintln!("Stack at start: {:?}", stack);

        while let Some(candidate) = stack.pop() {
            let rank = candidate.curr().rank();

            let offset = match candidate.curr().offset().try_into() {
                Ok(v) => v,
                Err(_) => panic!("Could not obtain offset!")
            } - 1;

            eprintln!("Popped item {:?}, node_rank: {}, seq offset: {:?}", candidate, rank, offset);

            if (offset as usize) < seq.len() && graph.is_symbol_equal(rank, seq[offset as usize]) {
                eprintln!("Sequence matches, try to extend offset.");

                if self.compute.extend(&candidate) {
                    // Sequence matches, and increased the furthest reaching point, add successors to stack
                    for new_point in self.neighboring_fr_points(graph, candidate.curr()) {
                        // Add neighbor with updated diagonal and one step further along `seq`
                        stack.push(candidate.new_for_successor(new_point));
                        eprintln!("Added neighbor {:?}", new_point);
                    }
                }
            }
        }
    }
}

impl<Compute: WFCompute> Default for WavefrontPOAligner<Compute> {
    fn default() -> Self {
        Self::new()
    }
}
