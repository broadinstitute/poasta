use std::{marker::PhantomData, ops::Bound};

use astar::{heuristic::AstarHeuristic, Astar, AstarResult, AstarState};
use cost_models::AlignmentCostModel;
use fr_points::{DiagType, PosType};
use traits::AlignableGraph;
use crate::errors::PoastaError;

pub mod traits;
pub mod astar;
pub mod cost_models;
pub(crate) mod extension;
pub(crate) mod fr_points;
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

pub trait GraphAligner<G>
where
    G: AlignableGraph,
{
    fn align<S>(
        &self,
        graph: &G,
        seq: S,
        mode: AlignmentMode,
    ) -> Result<AstarResult<G>, PoastaError>
    where
        S: AsRef<[u8]>;
}

pub struct PoastaAligner<H, D, O, C, G> {
    cost_model: C,
    dummy: PhantomData<(H, D, O, G)>,
}

impl<H, D, O, C, G> PoastaAligner<H, D, O, C, G>
where
    C: AlignmentCostModel,
    H: AstarHeuristic<
        G,
        <<C as AlignmentCostModel>::AstarStateType<G, D, O> as AstarState<G>>::AstarItem,
    >,
    G: AlignableGraph,
    D: DiagType,
    O: PosType,
{
    pub fn new(cost_model: C) -> Self {
        Self {
            cost_model,
            dummy: PhantomData,
        }
    }

    pub fn align_with_precomputed_heuristic(
        &self,
        graph: &G,
        seq: impl AsRef<[u8]>,
        alignment_mode: AlignmentMode,
        heuristic: H,
    ) -> Result<AstarResult<G>, PoastaError> {
        self.align_u8(graph, seq.as_ref(), alignment_mode, heuristic)
    }

    fn align_u8(
        &self,
        graph: &G,
        seq: &[u8],
        alignment_mode: AlignmentMode,
        mut heuristic: H,
    ) -> Result<AstarResult<G>, PoastaError> {
        heuristic.init(graph, seq);
        let astar_state = self
            .cost_model
            .initialize(graph, seq, alignment_mode, |item| heuristic.h(item));

        let mut astar = Astar::new(graph, seq, &heuristic, alignment_mode, astar_state);

        astar.run()
    }
}

impl<H, D, O, C, G> GraphAligner<G> for PoastaAligner<H, D, O, C, G>
where
    C: AlignmentCostModel,
    H: AstarHeuristic<
        G,
        <<C as AlignmentCostModel>::AstarStateType<G, D, O> as AstarState<G>>::AstarItem,
    >,
    G: AlignableGraph,
    D: DiagType,
    O: PosType,
{
    fn align<S>(
        &self,
        graph: &G,
        seq: S,
        alignment_mode: AlignmentMode,
    ) -> Result<AstarResult<G>, PoastaError>
    where
        S: AsRef<[u8]>,
    {
        let heuristic = H::default();
        self.align_u8(graph, seq.as_ref(), alignment_mode, heuristic)
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
    use super::{AlignmentMode, GraphAligner, PoastaAligner};

    #[test]
    fn test_alignment() {
        let cost_model = Affine::new(4, 6, 2);
        let mut graph = POASeqGraph::<u32>::new();

        let aligner = PoastaAligner::<Dijkstra, i32, u32, _, _>::new(cost_model);

        let mut reader = File::open("../tests/test2_from_abpoa.fa")
            .map(BufReader::new)
            .map(fasta::Reader::new)
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
                        File::create(format!("../tests/output/graph_for_{}.dot", name)).unwrap();
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
