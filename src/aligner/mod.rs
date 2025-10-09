use std::sync::Arc;
use std::{marker::PhantomData, ops::Bound};

use tracing::{debug, span, Level};

use crate::errors::PoastaError;
use crate::graph::bubbles::index::BubbleIndex;
use astar::{heuristic::AstarHeuristic, AstarResult};
use cost_models::AlignmentCostModel;
use traits::AlignableGraph;

pub mod astar;
pub mod cost_models;
pub(crate) mod extension;
pub(crate) mod fr_points;
pub mod traits;
pub mod utils;

/// Enum representing the kind of alignment to perform
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    /// Perform global alignment of the query sequence to the graph
    Global,

    /// Allow free indels at the beginning or end (optionally up to a given maximum),
    /// e.g., for semi-global alignment.
    EndsFree {
        qry_free_begin: Bound<usize>,
        qry_free_end: Bound<usize>,
        graph_free_begin: Bound<usize>,
        graph_free_end: Bound<usize>,
    },
}

pub struct PoastaAligner<H, C, G> {
    heuristic: H,
    dummy: PhantomData<(C, G)>,
}

impl<H, C, G> PoastaAligner<H, C, G>
where
    H: AstarHeuristic<C, G>,
    C: AlignmentCostModel,
    G: AlignableGraph,
{
    pub fn new(heuristic: H) -> Self {
        Self {
            heuristic,
            dummy: PhantomData,
        }
    }

    pub fn align(
        &self,
        graph: &G,
        seq: impl AsRef<[u8]>,
        alignment_mode: AlignmentMode,
    ) -> Result<AstarResult<G>, PoastaError> {
        let bubble_index = Arc::new(BubbleIndex::new(graph));
        self.align_u8(graph, seq.as_ref(), bubble_index, alignment_mode)
    }

    pub fn align_with_bubble_index(
        &self,
        graph: &G,
        seq: impl AsRef<[u8]>,
        bubble_index: Arc<BubbleIndex>,
        alignment_mode: AlignmentMode,
    ) -> Result<AstarResult<G>, PoastaError> {
        self.align_u8(graph, seq.as_ref(), bubble_index, alignment_mode)
    }

    fn align_u8(
        &self,
        graph: &G,
        seq: &[u8],
        bubble_index: Arc<BubbleIndex>,
        alignment_mode: AlignmentMode,
    ) -> Result<AstarResult<G>, PoastaError> {
        let mut runnable = self
            .heuristic
            .init(graph, seq, bubble_index, alignment_mode);

        let span = span!(Level::INFO, "astar_run");
        let _enter = span.enter();

        let mut result = AstarResult::default();

        let (end_score, end_point) = loop {
            debug!("pop");
            let Some(front) = runnable.pop_front() else {
                panic!("Empty queue before reaching end!")
            };

            {
                let visit_span = span!(Level::DEBUG, "visit");
                let _visit_enter = visit_span.enter();

                debug!(front=?front);

                if runnable.prune(&front) {
                    debug!("State pruned.");
                    continue;
                }

                if runnable.is_end(graph, &front) {
                    let score = runnable.get_score(&front);
                    debug!(score=?score, "Reached end.");

                    break (score, front);
                }

                result.num_visited += 1;
                runnable.relax(graph, seq, &front);
            }
        };

        result.score = end_score;
        result.alignment = runnable.backtrace(graph, &end_point);

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;

    use noodles::fasta;

    use crate::aligner::utils::print_alignment;
    use crate::graph::io::dot::graph_to_dot;
    use crate::graph::poa::POASeqGraph;

    use super::astar::heuristic::Dijkstra;
    use super::cost_models::affine::Affine;
    use super::{AlignmentMode, PoastaAligner};

    #[test]
    fn test_alignment() {
        let cost_model = Affine::<i32, u32>::new(4, 6, 2);
        let heuristic = Dijkstra::new(cost_model);
        let mut graph = POASeqGraph::<u32>::new();

        let aligner = PoastaAligner::new(heuristic);

        let mut reader = File::open("tests/test2_from_abpoa.fa")
            .map(BufReader::new)
            .map(fasta::io::Reader::new)
            .unwrap();

        let sequences: Vec<_> = reader.records().map(|v| v.unwrap()).collect();

        for record in &sequences {
            let name = std::str::from_utf8(record.name()).unwrap();
            if graph.is_empty() {
                graph
                    .add_aligned_sequence(
                        name,
                        record.sequence(),
                        vec![1; record.sequence().len()],
                        None,
                    )
                    .unwrap();
            } else {
                {
                    let mut writer =
                        File::create(format!("tests/output/graph_for_{}.dot", name)).unwrap();
                    graph_to_dot(&mut writer, &graph).unwrap();
                }

                let aln = aligner
                    .align(&graph, record.sequence(), AlignmentMode::Global)
                    .unwrap();

                let aln_str = print_alignment(&graph, record.sequence().as_ref(), &aln.alignment);
                eprintln!("{aln_str}");

                graph
                    .add_aligned_sequence(
                        name,
                        record.sequence(),
                        vec![1; record.sequence().len()],
                        Some(&aln.alignment),
                    )
                    .unwrap();
            }
        }
    }
}
