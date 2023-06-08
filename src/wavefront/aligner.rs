use std::fs::File;
use std::io::BufWriter;
use crate::graph::POAGraph;
use crate::alignment::Alignment;

use super::compute::WFCompute;
use super::fr_points::ExtendCandidate;
use super::{OffsetPrimitive, DiagIx, FRPoint, WFPoint};

use std::marker::PhantomData;
use std::path::PathBuf;


pub struct WavefrontPOAligner<'a, Offset, Compute>
    where
        Offset: OffsetPrimitive,
        Compute: WFCompute<Offset>,
{
    graph: &'a POAGraph,
    compute: Compute,
    debug_output_dir: PathBuf,
    dummy: PhantomData<Offset>
}

impl<'a, Offset, Compute> WavefrontPOAligner<'a, Offset, Compute>
    where
        Offset: OffsetPrimitive,
        Compute: WFCompute<Offset>,
{
    pub fn new(graph: &'a POAGraph, debug_output_dir: impl Into<PathBuf>) -> Self {
        WavefrontPOAligner {
            graph,
            compute: Compute::default(),
            debug_output_dir: debug_output_dir.into(),
            dummy: PhantomData
        }
    }

    pub fn align(&mut self, seq: &[u8]) -> Alignment {
        let max_offset = match Offset::max_value().try_into() {
            Ok(v) => v,
            Err(_) => panic!("Could not determine maximum value for Offset type!")
        };

        assert!(seq.len() - 1 < max_offset as usize, "Sequence is too long for Offset type!");

        let reached_end;
        let mut score: i64 = 0;
        loop {
            eprintln!("----- EXTEND: score {}", score);
            self.wf_extend(seq);

            let f_after_extend = self.debug_output_dir.join(format!("g{}.s{}.score{score}.after_extend.tsv", self.graph.graph.node_count(), seq.len()));
            let mut fout = File::create(f_after_extend).unwrap();
            let mut writer = BufWriter::new(fout);
            self.compute.write_csv(&mut writer).unwrap();
            drop(writer);

            if let Some(end_point) = self.compute.reached_end(self.graph, seq.len()) {
                eprintln!("Reached end point {:?}!", end_point);
                reached_end = Some(end_point);
                break;
            }

            score += 1;
            eprintln!("----- NEXT: score {}", score);
            self.compute.next(self.graph, seq.len(), score);

            let f_before_extend = self.debug_output_dir.join(format!("g{}.s{}.score{score}.before_extend.tsv", self.graph.graph.node_count(), seq.len()));
            fout = File::create(f_before_extend).unwrap();
            writer = BufWriter::new(fout);
            self.compute.write_csv(&mut writer).unwrap();
        }

        eprintln!("Sequence length: {:?}", seq.len());
        self.compute.backtrace(self.graph, reached_end.unwrap())
    }

    fn neighboring_fr_points(&self, point: FRPoint<Offset>) -> Vec<FRPoint<Offset>> {
        self.graph.neighbors_for_rank(point.rank())
            .map(|succ| {
                let succ_rank = self.graph.get_node_rank(succ);
                let succ_k = point.diag() - (succ_rank - point.rank() - 1) as DiagIx;

                (succ_k, point.offset() + Offset::one())
            })
            .collect()
    }

    fn wf_extend(&mut self, seq: &[u8]) {
        // Check if we can extend without mismatches along our current diagonals. Follow paths
        // in the graph in a depth first manner, and update furthest reaching points accordingly.
        let mut stack: Vec<ExtendCandidate<Offset>> = self.compute.extend_candidates().into_iter()
            .flat_map(|p| {
                self.neighboring_fr_points(p).into_iter()
                    .map(|neighbor| ExtendCandidate(neighbor, Some(p)))
                    .collect::<Vec<ExtendCandidate<Offset>>>()
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

            if (offset as usize) < seq.len() && self.graph.is_symbol_equal(rank, seq[offset as usize]) {
                eprintln!("Sequence matches, try to extend offset.");

                if self.compute.extend(&candidate) {
                    // Sequence matches, and increased the furthest reaching point, add successors to stack
                    for new_point in self.neighboring_fr_points(candidate.curr()) {
                        // Add neighbor with updated diagonal and one step further along `seq`
                        stack.push(candidate.new_for_successor(new_point));
                        eprintln!("Added neighbor {:?}", new_point);
                    }
                }
            }
        }
    }
}
