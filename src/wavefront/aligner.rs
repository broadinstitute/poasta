use crate::graph::POAGraph;
use super::compute::WFCompute;
use super::{OffsetPrimitive, DiagIx, FRPoint, WFPoint};

use std::marker::PhantomData;

pub struct WavefrontPOAligner<'a, Offset, Compute, A>
    where
        Offset: OffsetPrimitive,
        Compute: WFCompute<Offset>,
        A: Eq
{
    graph: &'a POAGraph<A>,
    compute: Compute,
    dummy: PhantomData<Offset>
}

impl<'a, Offset, Compute, A> WavefrontPOAligner<'a, Offset, Compute, A>
    where
        Offset: OffsetPrimitive,
        Compute: WFCompute<Offset>,
        A: Eq + Clone,
{
    pub fn new(graph: &'a POAGraph<A>) -> Self {
        WavefrontPOAligner {
            graph,
            compute: Compute::default(),
            dummy: PhantomData
        }
    }

    pub fn align(&mut self, seq: &[A]) {
        let max_offset = match Offset::max_value().try_into() {
            Ok(v) => v,
            Err(_) => panic!("Could not determine maximum value for Offset type!")
        };

        assert!(seq.len() < max_offset as usize, "Sequence is too long for Offset type!");

        let k_end = seq.len() as DiagIx - self.graph.graph.node_count() as i64;
        let offset_end = Offset::from(seq.len()).unwrap();
        let fr_end = (k_end, offset_end) as FRPoint<Offset>;
        eprintln!("Need to reach: {:?}", fr_end);

        let mut score: i64 = 0;
        loop {
            eprintln!("----- EXTEND: score {}", score);
            self.wf_extend(seq);

            if self.compute.reached_point(&fr_end) {
                eprintln!("Reached the end!");
                break;
            }

            score += 1;
            eprintln!("----- NEXT: score {}", score);
            self.compute.next(self.graph, score);
        }
    }

    fn wf_extend(&mut self, seq: &[A]) {
        // Check if we can extend without mismatches along our current diagonals. Follow paths
        // in the graph in a depth first manner, and update furthest reaching points accordingly.
        let mut stack: Vec<FRPoint<Offset>> = self.compute.extend_candidates();
        eprintln!("Stack at start: {:?}", stack);

        while let Some(point) = stack.pop() {
            let rank = point.rank();
            let node = self.graph.get_node_by_rank(rank);

            eprintln!("Popped item {:?}, node_rank: {}", point, rank);
            eprintln!("Remaining stack: {:?}", stack);

            // Sequence matches, add successors to stack
            let offset = match point.offset().try_into() {
                Ok(v) => v,
                Err(_) => panic!("Could not obtain offset!")
            };

            if self.graph.graph[node].code == seq[offset as usize] {
                self.compute.extend(&point);

                for neighbor in self.graph.graph.neighbors(node) {
                    // `succ_rank` is always greater than `rank` because of topological order
                    let succ_rank = self.graph.get_node_rank(neighbor);
                    let new_k = point.diag() - (succ_rank - rank - 1) as DiagIx;

                    // Add neighbor with updated diagonal and one step further along `seq`
                    if (offset as usize) < seq.len() - 1 {
                        stack.push((new_k, point.offset() + Offset::one()));
                        eprintln!("Added neighbor {:?}", (new_k, point.offset() + Offset::one()))
                    }
                }
            }
        }
    }
}
