//! Randomized differential test: every engine configuration against the DP
//! oracle, on small random POA graphs.
//!
//! Each case builds a graph from mutated copies of a random template and
//! aligns a mutated window of one copy, or one time in ten an unrelated
//! string. For every cost model, both modes, pruning on and off, and every
//! heuristic, the score must equal the oracle's, and for single-piece
//! models the backtrace must rescore to it. On sampled states `h` must not
//! exceed the oracle's remaining cost. The run also counts pruned states,
//! reopened states, and consistency violations of `h` over every edge of
//! the alignment graph, and prints them with `--nocapture`. Min-gap-cost
//! and the smooth seed heuristic are both expected to be consistent, so a
//! violation fails the test.
//!
//! `POASTA_DIFF_CASES` sets the number of cases (default 40 in release
//! builds, 8 in debug builds) and `POASTA_DIFF_SEED` the random seed.

use std::ops::Bound;
use std::sync::Arc;

use super::dp_oracle::DpOracle;
use super::heuristic::{AstarHeuristic, MinimumGapCost, SeedParams, SmoothHeuristic, SmoothHeuristicConfig};
use super::*;
use crate::align::PoastaAligner;
use crate::align::cost_models::linear::Linear;
use crate::align::cost_models::two_piece::TwoPieceAffine;
use crate::align::score::INF;
use crate::align::test_support::score_from_alignment;
use crate::graph::poa::{POAGraph, POANodeIndex};
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

type G = POAGraph<u32>;
type N = POANodeIndex<u32>;

/// Xorshift, so a failing case reproduces from its seed.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn chance(&mut self, p: f64) -> bool {
        (self.next() % 1_000_000) as f64 / 1e6 < p
    }
}

/// Substitute, delete or insert after each base with probability `rate`.
fn mutate(rng: &mut Rng, seq: &[u8], alphabet: &[u8], rate: f64) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len() + 4);
    for &b in seq {
        if rng.chance(rate) {
            match rng.below(3) {
                0 => out.push(alphabet[rng.below(alphabet.len())]),
                1 => {}
                _ => {
                    out.push(b);
                    out.push(alphabet[rng.below(alphabet.len())]);
                }
            }
        } else {
            out.push(b);
        }
    }
    if out.is_empty() {
        out.push(alphabet[0]);
    }
    out
}

/// A random POA graph and a query, with a description to paste into a
/// regression test. Templates of 8 to 57 bases over ACGT or AC (repeats),
/// 1 to 6 mutated copies.
fn random_case(rng: &mut Rng) -> (G, Vec<u8>, String) {
    let alphabet: &[u8] = if rng.below(2) == 0 { b"ACGT" } else { b"AC" };
    let tlen = 8 + rng.below(50);
    let template: Vec<u8> = (0..tlen).map(|_| alphabet[rng.below(alphabet.len())]).collect();
    let n_seqs = 1 + rng.below(6);
    let rate = rng.below(26) as f64 / 100.0;
    let seqs: Vec<Vec<u8>> = (0..n_seqs).map(|_| mutate(rng, &template, alphabet, rate)).collect();

    let engine = AstarEngine::<G, Affine>::new(Affine::new(4, 6, 2)).with_heuristic(DijkstraHeuristic);
    let mut aligner = PoastaAligner::new(engine, POAGraph::new());
    for (i, s) in seqs.iter().enumerate() {
        aligner.align_one_named(&format!("s{i}"), s).unwrap();
    }
    let graph = aligner.into_graph();

    let query = if rng.chance(0.1) {
        (0..1 + rng.below(30)).map(|_| alphabet[rng.below(alphabet.len())]).collect()
    } else {
        let src = &seqs[rng.below(n_seqs)];
        let a = rng.below(src.len());
        let b = a + 1 + rng.below(src.len() - a);
        let qrate = rng.below(31) as f64 / 100.0;
        mutate(rng, &src[a..b], alphabet, qrate)
    };
    let desc = format!(
        "seqs={:?} query={:?}",
        seqs.iter().map(|s| String::from_utf8_lossy(s).into_owned()).collect::<Vec<_>>(),
        String::from_utf8_lossy(&query)
    );
    (graph, query, desc)
}

fn semi_global() -> AlignmentMode {
    AlignmentMode::EndsFree {
        qry_free_begin: Bound::Excluded(0),
        qry_free_end: Bound::Excluded(0),
        graph_free_begin: Bound::Unbounded,
        graph_free_end: Bound::Unbounded,
    }
}

fn short_name<T>() -> &'static str {
    std::any::type_name::<T>().rsplit("::").next().unwrap()
}

#[derive(Default)]
struct Tally {
    alignments: usize,
    failures: Vec<String>,
    pruned: usize,
    reopened: usize,
    /// Alignments that both pruned and reopened states.
    both: usize,
    sampled_states: usize,
    inadmissible: Vec<String>,
    edges: usize,
    violations: usize,
    max_violation: usize,
    /// Heuristics with a consistency violation, which the proofs rule out.
    inconsistent: Vec<String>,
}

/// Align with one configuration and compare with the oracle's score.
#[allow(clippy::too_many_arguments)]
fn check<C, H>(t: &mut Tally, g: &G, q: &[u8], costs: C, mode: AlignmentMode, h: H, prune: bool, expected: u32, desc: &str)
where
    C: AstarCostModel<N, Offset> + std::fmt::Debug,
    H: HeuristicConfig,
{
    let mut engine = AstarEngine::<G, C>::new(costs).with_mode(mode).with_heuristic(h).with_reopen_tracking();
    if !prune {
        engine = engine.without_pruning();
    }
    let out = engine.align(g, q).unwrap();
    t.alignments += 1;
    t.pruned += out.stats.states_pruned;
    t.reopened += out.stats.states_reopened;
    if out.stats.states_pruned > 0 && out.stats.states_reopened > 0 {
        t.both += 1;
    }
    let label = format!("{} {costs:?} prune={prune} {mode:?}", short_name::<H>());
    if out.score != expected {
        t.failures.push(format!("score {} want {expected}: {label} {desc}", out.score));
    }
    // score_from_alignment knows one gap piece only.
    if costs.gap_extend2() == 0 {
        let rescored = score_from_alignment(g, q, &out.alignment, &costs);
        if rescored != out.score {
            t.failures.push(format!("backtrace rescores to {rescored}, score {}: {label} {desc}", out.score));
        }
    }
}

/// Successors of state (v, i, s) with their costs, as the DP oracle defines
/// the alignment graph. `pieces` holds (open, extend, deletion slot,
/// insertion slot) per gap piece.
#[allow(clippy::too_many_arguments)]
fn successors(
    g: &G,
    q: &[u8],
    mismatch: u32,
    pieces: &[(u32, u32, usize, usize)],
    states: usize,
    v: N,
    i: usize,
    s: usize,
) -> Vec<(N, usize, usize, u32)> {
    let mut out = Vec::new();
    let linear = states == 1;
    if s == 0 {
        for u in g.successors(v) {
            if u == g.end_node() {
                out.push((u, i, 0, 0));
            } else if i < q.len() {
                let c = if g.is_symbol_equal(u, q[i]) { 0 } else { mismatch };
                out.push((u, i + 1, 0, c));
            }
            for &(go, ge, del, _) in pieces {
                out.push((u, i, if linear { 0 } else { del }, go + ge));
            }
        }
        if i < q.len() {
            for &(go, ge, _, ins) in pieces {
                out.push((v, i + 1, if linear { 0 } else { ins }, go + ge));
            }
        }
    } else {
        let (_, ge, del, _) = *pieces.iter().find(|p| p.2 == s || p.3 == s).unwrap();
        if s == del {
            for u in g.successors(v) {
                out.push((u, i, s, ge));
            }
        } else if i < q.len() {
            out.push((v, i + 1, s, ge));
        }
        out.push((v, i, 0, 0));
    }
    out
}

const STATES: [AlignState; 5] =
    [AlignState::Match, AlignState::Deletion, AlignState::Insertion, AlignState::Deletion2, AlignState::Insertion2];

/// Sample 12 states for `h ≤ cost_from`, and count consistency violations
/// of `h` over every edge of the alignment graph.
#[allow(clippy::too_many_arguments)]
fn check_heuristic<C, Hs>(
    t: &mut Tally,
    g: &G,
    q: &[u8],
    costs: C,
    mode: AlignmentMode,
    h: &mut Hs,
    name: &str,
    rng: &mut Rng,
    desc: &str,
) where
    C: AstarCostModel<N, Offset> + std::fmt::Debug,
    Hs: AstarHeuristic<N, Offset>,
{
    let oracle: DpOracle<C, G> =
        if mode == AlignmentMode::Global { DpOracle::new(costs) } else { DpOracle::new_semi_global(costs) };
    let states = C::STATES;
    let mut pieces = vec![(costs.gap_open() as u32, costs.gap_extend() as u32, 1usize, 2usize)];
    if costs.gap_extend2() > 0 {
        pieces.push((costs.gap_open2() as u32, costs.gap_extend2() as u32, 3, 4));
    }
    let nodes: Vec<N> = g.all_nodes_iter().collect();
    let mut hv = |v: N, i: usize, s: usize| h.h(g, &AlignmentGraphNode::new(v, i as Offset), STATES[s]);
    let violations_before = t.violations;

    for _ in 0..12 {
        let (v, i, s) = (nodes[rng.below(nodes.len())], rng.below(q.len() + 1), rng.below(states));
        let truth = oracle.cost_from(g, q, v, i, s);
        if truth >= INF {
            continue;
        }
        t.sampled_states += 1;
        let est = hv(v, i, s);
        if est as u32 > truth {
            t.inadmissible.push(format!("{name} {costs:?} {mode:?} state=({v:?},{i},{s}) h={est} true={truth} {desc}"));
        }
    }
    for &v in &nodes {
        for i in 0..=q.len() {
            for s in 0..states {
                let hx = hv(v, i, s);
                for (u, j, s2, c) in successors(g, q, costs.mismatch() as u32, &pieces, states, v, i, s) {
                    t.edges += 1;
                    let hy = hv(u, j, s2);
                    if hx > c as usize + hy {
                        t.violations += 1;
                        t.max_violation = t.max_violation.max(hx - c as usize - hy);
                    }
                }
            }
        }
    }
    if t.violations > violations_before {
        t.inconsistent.push(format!("{name} {costs:?} {mode:?}: {} violations {desc}", t.violations - violations_before));
    }
}

fn run_case<C>(t: &mut Tally, g: &G, q: &[u8], costs: C, rng: &mut Rng, desc: &str)
where
    C: AstarCostModel<N, Offset> + std::fmt::Debug,
{
    for mode in [AlignmentMode::Global, semi_global()] {
        let mut oracle: DpOracle<C, G> =
            if mode == AlignmentMode::Global { DpOracle::new(costs) } else { DpOracle::new_semi_global(costs) };
        let expected = oracle.align(g, q).unwrap().score;
        for prune in [false, true] {
            check(t, g, q, costs, mode, DijkstraHeuristic, prune, expected, desc);
            check(t, g, q, costs, mode, MinGapCostHeuristic, prune, expected, desc);
            for k in 2..=5 {
                for r in 0..=1 {
                    for prune_matches in [false, true] {
                        let params = SeedParams { k, r, prune_matches };
                        let label = format!("{params:?} {desc}");
                        check(t, g, q, costs, mode, SmoothHeuristicConfig(params), prune, expected, &label);
                    }
                }
            }
        }
        let bubbles = Arc::new(BubbleIndex::new(g));
        let mut mingap = MinimumGapCost::<N, Offset, C>::new_with_mode(costs, bubbles.clone(), q.len(), mode);
        check_heuristic(t, g, q, costs, mode, &mut mingap, "min-gap-cost", rng, desc);
        let params = SeedParams { k: 2 + rng.below(4), r: rng.below(2), prune_matches: false };
        let mut smooth = SmoothHeuristic::<N, Offset, C>::new(costs, params);
        smooth.prepare(g, q, &bubbles, mode);
        check_heuristic(t, g, q, costs, mode, &mut smooth, &format!("smooth {params:?}"), rng, desc);
    }
}

#[test]
fn differential_against_oracle() {
    let default_cases = if cfg!(debug_assertions) { 8 } else { 40 };
    let cases: usize = std::env::var("POASTA_DIFF_CASES").ok().and_then(|v| v.parse().ok()).unwrap_or(default_cases);
    let seed: u64 = std::env::var("POASTA_DIFF_SEED").ok().and_then(|v| v.parse().ok()).unwrap_or(0x9E37_79B9_7F4A_7C15);
    let mut rng = Rng(seed);
    let mut t = Tally::default();
    for _ in 0..cases {
        let (g, q, desc) = random_case(&mut rng);
        run_case(&mut t, &g, &q, Affine::new(4, 6, 2), &mut rng, &desc);
        run_case(&mut t, &g, &q, Affine::new(1, 2, 1), &mut rng, &desc);
        run_case(&mut t, &g, &q, Linear::new(4, 2), &mut rng, &desc);
        run_case(&mut t, &g, &q, TwoPieceAffine::new(4, 6, 2, 24, 1), &mut rng, &desc);
        run_case(&mut t, &g, &q, TwoPieceAffine::new(1, 2, 2, 6, 1), &mut rng, &desc);
    }
    eprintln!(
        "cases={cases} alignments={} failures={} pruned={} reopened={} both={} \
         sampled={} inadmissible={} edges={} consistency violations={} (max {})",
        t.alignments,
        t.failures.len(),
        t.pruned,
        t.reopened,
        t.both,
        t.sampled_states,
        t.inadmissible.len(),
        t.edges,
        t.violations,
        t.max_violation,
    );
    for f in t.inadmissible.iter().chain(&t.inconsistent).chain(&t.failures).take(10) {
        eprintln!("{f}");
    }
    assert!(t.failures.is_empty() && t.inadmissible.is_empty() && t.inconsistent.is_empty());
}
