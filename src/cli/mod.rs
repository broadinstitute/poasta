//! Command-line argument definitions for the `poasta` and `lasagna` binaries
//! (and `poasta_vs_spoa` behind the `bench` feature). The drivers themselves
//! live in `src/bin/`.

use clap::{Args, ValueEnum};

use std::ops::Bound;

use crate::align::AlignmentMode;
use crate::align::astar::SeedParams;

pub mod lasagna;
pub mod poasta;

#[cfg(feature = "bench")]
pub mod poasta_vs_spoa;

/// An enum indicating what kind of alignment to perform. Shared between
/// `poasta align` and `lasagna`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum AlignmentSpan {
    /// Perform global alignment
    Global,

    /// Perform semi-global alignment, i.e., globally align query but allow free gaps in the graph
    /// at the beginning and end
    SemiGlobal,
}

impl AlignmentSpan {
    /// The engine-level mode. `SemiGlobal` frees the graph's leading and
    /// trailing nodes and aligns the query end to end.
    pub fn mode(self) -> AlignmentMode {
        match self {
            AlignmentSpan::Global => AlignmentMode::Global,
            AlignmentSpan::SemiGlobal => AlignmentMode::EndsFree {
                qry_free_begin: Bound::Excluded(0),
                qry_free_end: Bound::Excluded(0),
                graph_free_begin: Bound::Unbounded,
                graph_free_end: Bound::Unbounded,
            },
        }
    }
}

/// Log verbosity levels exposed on the CLI.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl LogLevel {
    pub fn as_tracing_str(self) -> &'static str {
        match self {
            LogLevel::Error => "error",
            LogLevel::Warn => "warn",
            LogLevel::Info => "info",
            LogLevel::Debug => "debug",
            LogLevel::Trace => "trace",
        }
    }
}

/// Cost model selection.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum CostModelKind {
    /// Affine gap cost (single piecewise-linear gap penalty).
    Affine,

    /// Linear gap cost (no gap-open penalty).
    Linear,

    /// Two-piece affine gap cost.
    TwoPiece,
}

/// Which A* heuristic to use. Shared between `poasta align` and `lasagna`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum HeuristicKind {
    /// Admissible min-gap-cost lower bound (default): tighter bound, visits
    /// fewer states than Dijkstra.
    MinGapCost,

    /// Plain Dijkstra (`h ≡ 0`): exact, but explores every state within the
    /// optimal-score radius.
    Dijkstra,

    /// Gap-chaining seed heuristic: charges for query seeds that cannot be
    /// matched ahead, and prunes matches the search has passed. Much faster
    /// on long reads and on short reads against a large graph.
    Chain,
}

/// A* search options shared by `poasta align` and `lasagna`.
#[derive(Args, Debug, Clone)]
pub struct SearchArgs {
    /// Which A* heuristic to use.
    #[arg(long, value_enum, default_value = "min-gap-cost")]
    #[clap(help_heading = "Search")]
    pub heuristic: HeuristicKind,

    /// Seed length for `--heuristic chain`. Longer seeds match less by
    /// chance, shorter ones survive more sequencing errors.
    #[arg(long, default_value_t = 8, value_parser = clap::value_parser!(u8).range(4..=32))]
    #[clap(help_heading = "Search")]
    pub seed_length: u8,

    /// Substitutions allowed inside one seed match for `--heuristic chain`,
    /// 0 or 1.
    #[arg(long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=1))]
    #[clap(help_heading = "Search")]
    pub seed_mismatches: u8,

    /// Turn off bubble-exit pruning.
    #[arg(long)]
    #[clap(help_heading = "Search")]
    pub no_pruning: bool,

    /// Turn off seed-match pruning. For diagnosis: without it the chain
    /// heuristic is quadratic on long reads.
    #[arg(long)]
    #[clap(help_heading = "Search")]
    pub no_match_pruning: bool,
}

impl SearchArgs {
    pub fn seed_params(&self) -> SeedParams {
        SeedParams {
            k: self.seed_length as usize,
            r: self.seed_mismatches as usize,
            prune_matches: !self.no_match_pruning,
        }
    }
}

/// Instantiate engines with the search options in `SearchArgs`.
/// `$new_plain_engine` is a closure that creates an engine without them.
/// `$body` runs with `$new_engine` bound to a closure that creates one with
/// them: the chosen heuristic, and pruning unless turned off. A macro,
/// because each heuristic gives the engine another type and a closure cannot
/// be generic over it.
///
/// ```ignore
/// instantiate_engine!(&args.search, || AstarEngine::new(costs), |new_engine| run(new_engine()))
/// ```
#[macro_export]
macro_rules! instantiate_engine {
    ($search:expr, $new_plain_engine:expr, |$new_engine:ident| $body:expr) => {{
        let search: &$crate::cli::SearchArgs = $search;
        let new_plain_engine = $new_plain_engine;
        let new_pruned_engine = || {
            let engine = new_plain_engine();
            if search.no_pruning {
                engine.without_pruning()
            } else {
                engine
            }
        };
        match search.heuristic {
            $crate::cli::HeuristicKind::MinGapCost => {
                let $new_engine = new_pruned_engine;
                $body
            }
            $crate::cli::HeuristicKind::Dijkstra => {
                let $new_engine =
                    || new_pruned_engine().with_heuristic($crate::align::astar::DijkstraHeuristic);
                $body
            }
            $crate::cli::HeuristicKind::Chain => {
                let config = $crate::align::astar::SmoothHeuristicConfig(search.seed_params());
                let $new_engine = || new_pruned_engine().with_heuristic(config);
                $body
            }
        }
    }};
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::HeuristicKind;
    use crate::cli::lasagna::LasagnaCli;
    use crate::cli::poasta::{CliArgs, CliSubcommand};

    fn align(extra: &[&str]) -> Result<crate::cli::SearchArgs, clap::Error> {
        let mut argv = vec!["poasta", "align", "in.fa"];
        argv.extend_from_slice(extra);
        match CliArgs::try_parse_from(argv)?.command {
            Some(CliSubcommand::Align(args)) => Ok(args.search),
            _ => unreachable!(),
        }
    }

    #[test]
    fn search_defaults_and_ranges() {
        let s = align(&[]).unwrap();
        assert_eq!(s.heuristic, HeuristicKind::MinGapCost);
        assert_eq!(
            (
                s.seed_length,
                s.seed_mismatches,
                s.no_pruning,
                s.no_match_pruning
            ),
            (8, 0, false, false)
        );
        assert!(s.seed_params().prune_matches);

        let s = align(&[
            "--heuristic",
            "chain",
            "--seed-length",
            "12",
            "--seed-mismatches",
            "1",
            "--no-match-pruning",
        ])
        .unwrap();
        assert_eq!(s.heuristic, HeuristicKind::Chain);
        let p = s.seed_params();
        assert_eq!((p.k, p.r, p.prune_matches), (12, 1, false));

        assert!(align(&["--seed-length", "3"]).is_err());
        assert!(align(&["--seed-mismatches", "2"]).is_err());
    }

    #[test]
    fn lasagna_accepts_semi_global() {
        let cli =
            LasagnaCli::try_parse_from(["lasagna", "g.gfa", "r.fa", "-m", "semi-global"]).unwrap();
        assert!(matches!(
            cli.alignment_span.mode(),
            crate::align::AlignmentMode::EndsFree { .. }
        ));
    }

    #[test]
    fn lasagna_takes_the_same_search_options() {
        let cli = LasagnaCli::try_parse_from([
            "lasagna",
            "g.gfa",
            "r.fa",
            "--heuristic",
            "chain",
            "--no-pruning",
        ])
        .unwrap();
        assert_eq!(cli.search.heuristic, HeuristicKind::Chain);
        assert!(cli.search.no_pruning);
    }
}
