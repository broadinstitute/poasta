//! Aligning a query sequence to a graph.
//!
//! [`PoastaAligner`] drives the batch flow: align one sequence with an
//! [`AlignmentEngine`], fold the result into the graph with `AddAlignment`,
//! repeat. [`AlignmentMode`] and [`AlignState`] are shared by the engine and
//! its tests. The engine itself lives in [`astar`], the cost models in
//! [`cost_models`].

use std::fmt;
use std::ops::Bound;

use crate::graph::alignment::AddAlignment;

use self::traits::AlignmentStats;
use self::traits::{AlignResult, AlignableGraph, AlignmentEngine};

pub mod cost_models;
pub mod astar;
pub mod offsets;
pub mod score;
pub mod traits;

/// Test-only helper for verifying that a backtrace-produced alignment is
/// self-consistent with its reported score, independent of the DP oracle.
#[cfg(test)]
pub(crate) mod test_support {
    use super::cost_models::AlignmentCostModel;
    use super::traits::AlignedPair;
    use super::score::sat;
    use super::traits::AlignableGraph;

    /// Recompute an alignment's score by replaying its `AlignedPair`s through
    /// the cost model, independent of however the engine computed its DP
    /// scores. Assumes a single-piece gap-affine cost model (`gap_open`/
    /// `gap_extend`), so it applies to `Affine` but not `TwoPieceAffine`.
    pub(crate) fn score_from_alignment<G, C>(
        graph: &G,
        query: &[u8],
        aln: &[AlignedPair<G::Node>],
        costs: &C,
    ) -> u32
    where
        G: AlignableGraph,
        C: AlignmentCostModel,
    {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        let mut score = 0u32;
        let mut in_gap = false;

        for pair in aln {
            match (pair.node(), pair.query_pos()) {
                (Some(node), Some(qpos)) => {
                    let gsym = graph.get_node_symbol(node);
                    let qsym = query[qpos];
                    score = sat(
                        score,
                        if gsym == qsym {
                            0u32
                        } else {
                            costs.mismatch() as u32
                        },
                    );
                    in_gap = false;
                }
                (Some(_), None) | (None, Some(_)) => {
                    score = sat(score, if in_gap { ge } else { go + ge });
                    in_gap = true;
                }
                (None, None) => {}
            }
        }
        score
    }
}

/// Represents the kind of alignment to perform
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

/// Represents the alignment state of a particular cell in the alignment matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum AlignState {
    Match,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model
    Insertion2,
}

/// Error returned by [`PoastaAligner`] when either the engine or the graph fails.
#[derive(Debug)]
pub enum PoastaAlignerError<EErr, GErr> {
    Engine(EErr),
    AddAlignment(GErr),
}

impl<EErr, GErr> fmt::Display for PoastaAlignerError<EErr, GErr>
where
    EErr: fmt::Display,
    GErr: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Engine(e) => write!(f, "alignment engine error: {e}"),
            Self::AddAlignment(e) => write!(f, "graph add_alignment error: {e}"),
        }
    }
}

impl<EErr, GErr> std::error::Error for PoastaAlignerError<EErr, GErr>
where
    EErr: std::error::Error + 'static,
    GErr: std::error::Error + 'static,
{
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Engine(e) => Some(e),
            Self::AddAlignment(e) => Some(e),
        }
    }
}

/// Running totals across a batch of alignments.
///
/// Populated by [`PoastaAligner`] after every successful per-sequence
/// alignment so callers can report averages once the run finishes.
#[derive(Debug, Default, Clone, Copy)]
pub struct RunStats {
    pub n_alignments: usize,
    pub sum_cells_computed: u128,
    pub sum_fraction: f64,
}

impl RunStats {
    fn record(&mut self, stats: &AlignmentStats) {
        self.n_alignments += 1;
        self.sum_cells_computed += stats.cells_computed as u128;
        self.sum_fraction += stats.fraction_of_full_matrix;
    }

    pub fn avg_cells_computed(&self) -> f64 {
        if self.n_alignments == 0 {
            0.0
        } else {
            self.sum_cells_computed as f64 / self.n_alignments as f64
        }
    }

    pub fn avg_fraction(&self) -> f64 {
        if self.n_alignments == 0 {
            0.0
        } else {
            self.sum_fraction / self.n_alignments as f64
        }
    }
}

/// High-level batch aligner that owns an alignment engine and a graph, and
/// incrementally builds up the graph by aligning sequences one at a time.
pub struct PoastaAligner<E, G> {
    engine: E,
    graph: G,
    run_stats: RunStats,
}

impl<E, G> PoastaAligner<E, G>
where
    G: AlignableGraph,
{
    pub fn new(engine: E, graph: G) -> Self {
        Self {
            engine,
            graph,
            run_stats: RunStats::default(),
        }
    }

    pub fn graph(&self) -> &G {
        &self.graph
    }

    pub fn run_stats(&self) -> &RunStats {
        &self.run_stats
    }

    pub fn into_graph(self) -> G {
        self.graph
    }

    /// Align each sequence to the current graph and fold the alignment back in.
    ///
    /// Sequence names are synthesised as `seq_0`, `seq_1`, ... and every base
    /// gets a weight of 1.
    pub fn align_all<'s, I>(
        &mut self,
        sequences: I,
    ) -> Result<(), PoastaAlignerError<<E as AlignmentEngine<&'s [u8]>>::Error, G::Error>>
    where
        I: IntoIterator<Item = &'s [u8]>,
        E: AlignmentEngine<&'s [u8], Graph = G>,
        G: AddAlignment<<E as AlignmentEngine<&'s [u8]>>::Success>,
        <E as AlignmentEngine<&'s [u8]>>::Success: AlignResult,
    {
        for (i, seq) in sequences.into_iter().enumerate() {
            let name = format!("seq_{i}");
            let weights = vec![1usize; seq.len()];
            self.align_one(&name, seq, &weights)?;
        }
        Ok(())
    }

    /// Align each `(name, sequence)` pair, using weights of 1 for every base.
    pub fn align_all_named<'s, I>(
        &mut self,
        sequences: I,
    ) -> Result<(), PoastaAlignerError<<E as AlignmentEngine<&'s [u8]>>::Error, G::Error>>
    where
        I: IntoIterator<Item = (&'s str, &'s [u8])>,
        E: AlignmentEngine<&'s [u8], Graph = G>,
        G: AddAlignment<<E as AlignmentEngine<&'s [u8]>>::Success>,
        <E as AlignmentEngine<&'s [u8]>>::Success: AlignResult,
    {
        for (name, seq) in sequences {
            let weights = vec![1usize; seq.len()];
            self.align_one(name, seq, &weights)?;
        }
        Ok(())
    }

    /// Align each `(name, sequence, weights)` triple. `weights.len()` must equal
    /// `sequence.len()`.
    pub fn align_all_with_weights<'s, I>(
        &mut self,
        sequences: I,
    ) -> Result<(), PoastaAlignerError<<E as AlignmentEngine<&'s [u8]>>::Error, G::Error>>
    where
        I: IntoIterator<Item = (&'s str, &'s [u8], &'s [usize])>,
        E: AlignmentEngine<&'s [u8], Graph = G>,
        G: AddAlignment<<E as AlignmentEngine<&'s [u8]>>::Success>,
        <E as AlignmentEngine<&'s [u8]>>::Success: AlignResult,
    {
        for (name, seq, weights) in sequences {
            self.align_one(name, seq, weights)?;
        }
        Ok(())
    }

    /// Align one named sequence using a weight of 1 for every base.
    pub fn align_one_named<'s>(
        &mut self,
        name: &str,
        seq: &'s [u8],
    ) -> Result<(), PoastaAlignerError<<E as AlignmentEngine<&'s [u8]>>::Error, G::Error>>
    where
        E: AlignmentEngine<&'s [u8], Graph = G>,
        G: AddAlignment<<E as AlignmentEngine<&'s [u8]>>::Success>,
        <E as AlignmentEngine<&'s [u8]>>::Success: AlignResult,
    {
        let weights = vec![1usize; seq.len()];
        self.align_one(name, seq, &weights)
    }

    fn align_one<'s>(
        &mut self,
        name: &str,
        seq: &'s [u8],
        weights: &[usize],
    ) -> Result<(), PoastaAlignerError<<E as AlignmentEngine<&'s [u8]>>::Error, G::Error>>
    where
        E: AlignmentEngine<&'s [u8], Graph = G>,
        G: AddAlignment<<E as AlignmentEngine<&'s [u8]>>::Success>,
        <E as AlignmentEngine<&'s [u8]>>::Success: AlignResult,
    {
        let span = tracing::info_span!("align_seq", seq_id = name);
        let _enter = span.enter();

        if self.graph.is_empty() {
            tracing::debug!(
                "Initialize graph from first sequence {name} (len: {})",
                seq.len()
            );

            self.graph
                .add_alignment(name, seq, None, weights)
                .map_err(PoastaAlignerError::AddAlignment)?;
        } else {
            tracing::debug!("Aligning {name} (len: {})", seq.len());
            let result = self
                .engine
                .align(&self.graph, seq)
                .map_err(PoastaAlignerError::Engine)?;
            self.run_stats.record(&result.alignment_stats());
            tracing::debug!(score = result.alignment_score(), "Alignment complete");

            tracing::debug!("Updating graph...");
            self.graph
                .add_alignment(name, seq, Some(&result), weights)
                .map_err(PoastaAlignerError::AddAlignment)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cost_models::affine::Affine;
    use crate::align::astar::AstarEngine;
    use crate::graph::poa::POAGraph;
    use crate::graph::traits::GraphBase;

    #[test]
    fn align_all_builds_graph_from_empty() {
        let graph = POAGraph::<u32>::new();
        let engine = AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(1, 2, 1));
        let mut aligner = PoastaAligner::new(engine, graph);

        let seqs: &[&[u8]] = &[b"ACGT", b"ACCT", b"ACGG"];
        aligner
            .align_all(seqs.iter().copied())
            .expect("align_all should succeed");

        assert_eq!(aligner.graph().sequences.len(), 3);
        assert!(aligner.graph().node_count() > 2);
    }

    #[test]
    fn align_all_with_weights_accepts_custom_weights() {
        let graph = POAGraph::<u32>::new();
        let engine = AstarEngine::<POAGraph<u32>, Affine>::new(Affine::new(1, 2, 1));
        let mut aligner = PoastaAligner::new(engine, graph);

        let s0: &[u8] = b"ACGT";
        let s1: &[u8] = b"ACCT";
        let w0 = vec![2usize; s0.len()];
        let w1 = vec![3usize; s1.len()];

        let inputs: Vec<(&str, &[u8], &[usize])> = vec![("first", s0, &w0), ("second", s1, &w1)];

        aligner
            .align_all_with_weights(inputs)
            .expect("weighted align_all should succeed");

        assert_eq!(aligner.graph().sequences.len(), 2);
    }
}
