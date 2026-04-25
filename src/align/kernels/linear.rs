//! # LinearKernel — 1 state, linear gap cost (`k·ge`)
//!
//! All three DP states (match, deletion, insertion) collapse into a single
//! per-cell score because a linear gap cost has no notion of "opening": every
//! additional step costs exactly `ge` regardless of history.
//!
//! ## Data layout (STATES = 1, w = qhi - qlo + 1)
//!
//! ```text
//!                    ← w query positions →
//!          qi:    0     1     2   ...  w-1
//!          q:   qlo  qlo+1 qlo+2  ... qhi
//!                ┌─────┬─────┬─────┬─────┐
//!     data[.] =  │  ·  │  ·  │  ·  │  ·  │   single state 0 (combined M/D/I)
//!                └─────┴─────┴─────┴─────┘
//! ```
//!
//! ## DP recurrence
//!
//! For graph node `v` with symbol `s`, predecessor `u`, query position `q`:
//!
//! ```text
//!   C[v][q] = min(
//!     C[u][q-1] + subst(s, query[q-1]),   (a) diagonal (match/mismatch)
//!     C[u][q]   + ge,                     (b) deletion (consume graph only)
//!     C[v][q-1] + ge,                     (c) insertion (consume query only)
//!   )
//! ```
//!
//! ## Computation pipeline (band v)
//!
//! ```text
//!   step 1  accumulate_gap(u → v)     — (b) deletion from predecessor u at q
//!   step 2  accumulate_diag(u → v)    — (a) best predecessor score at q-1
//!   step 3  finalize                  — apply subst, then left→right (c) sweep
//! ```
//!
//! ### step 1: `accumulate_gap` — one call per predecessor band `u`
//!
//! ```text
//!     pred band u   . X . . .    (state 0 at query position q)
//!                     │
//!                     │  + ge
//!                     ▼
//!     curr band v   . Y . . .    Y = min(Y, X + ge)
//!                     ▲
//!                     └─── same q, graph step v ← u
//! ```
//!
//! ### step 2: `accumulate_diag` — best predecessor at `q-1`
//!
//! ```text
//!     pred band u   X . . . .    (at q-1)
//!                    ╲
//!                     ▼ diagonal (no subst yet)
//!     diag_buf      . Y . .    Y = min(Y, X)
//! ```
//!
//! ### step 3: `finalize` — apply subst, then left→right insertion sweep
//!
//! ```text
//!     after accumulate_gap:   data    = [ g0  g1  g2  g3 ]   (deletion cost per cell)
//!     after accumulate_diag
//!     + substitution:         diag_buf = [ d0  d1  d2  d3 ]   (diag + subst)
//!
//!     qi = 0:     data[0] = min(g0, d0)                       (no left neighbour)
//!
//!     qi = 1..w:  data[qi] = min( data[qi] ,   (already has (b))
//!                                 diag_buf[qi],  ((a) diagonal)
//!                                 data[qi-1] + ge )  ((c) insertion)
//!
//!          [●]───ge──▶[●]───ge──▶[●]───ge──▶[●]
//!          qi=0       qi=1       qi=2       qi=3
//!     (insertion sweeps left→right within the band)
//! ```
//!
//! ## Backtrace note
//!
//! Because M and D share the same cell value, the current cell alone cannot always
//! distinguish "diagonal with subst = 0" from "deletion + ge"; when insertion is
//! ruled out we return `Diagonal` and let the engine's generic predecessor search
//! resolve ties — see [`LinearKernel::backtrace_op`].

use super::{BacktraceOp, DPKernel};
use crate::align::{
    cost_models::{AlignmentCostModel, linear::Linear},
    engine::dp::sat,
};

/// 1-state linear gap kernel: single state (index 0 = M/D/I combined).
/// Data layout: `data[qi]` (STATES=1, so no state stride).
pub struct LinearKernel;

impl DPKernel for LinearKernel {
    const STATES: usize = 1;
    const STATE_NAMES: &'static [&'static str] = &["M"];
    type Costs = Linear;

    /// Seed the start-sentinel band with linear gap costs along the query axis.
    ///
    /// The start node has implicit score 0 at `q = 0`, and every query position
    /// requires `q` insertions (each costing `ge`) to reach:
    ///
    /// ```text
    ///     qi:    0    1    2    3   ...   w-1
    ///     q :  qlo   ...                  qhi
    ///          ┌────┬────┬────┬────┬─  ─┬──────┐
    ///   data = │  0 │ ge │ 2ge│ 3ge│ ...│ q·ge │
    ///          └────┴────┴────┴────┴─  ─┴──────┘
    ///            ^ only qi=0 && qlo==0 is zero
    /// ```
    fn init_start(data: &mut [u32], w: usize, qlo: usize, costs: &Linear) {
        let ge = costs.gap_extend() as u32;
        for qi in 0..w {
            let q = qlo + qi;
            data[qi] = if q == 0 { 0 } else { ge * q as u32 };
        }
    }

    /// Fold deletion (graph-only step) contributions from predecessor band `src`
    /// into `dst` at the same query position.
    ///
    /// ```text
    ///     pred band u   . . X . . .     (src[pqi], same q = dst_qlo + qi)
    ///                       │
    ///                       │ + ge        (linear deletion cost)
    ///                       ▼
    ///     curr band v   . . Y . . .     Y ← min(Y, X + ge)
    /// ```
    ///
    /// Iterates over `q ∈ [max(dst_qlo, src_qlo), min(dst_qhi, src_qhi)]`, the
    /// query-range overlap between the two bands. No-op when they don't overlap.
    fn accumulate_gap(
        dst: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Linear,
    ) {
        let ge = costs.gap_extend() as u32;
        let q_lo = dst_qlo.max(src_qlo);
        let q_hi = (dst_qlo + w_dst - 1).min(src_qlo + w_src - 1);
        if q_lo > q_hi {
            return;
        }
        for q in q_lo..=q_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            // D[q] = predecessor M[q] + ge (one deletion step)
            dst[qi] = dst[qi].min(sat(src[pqi], ge));
        }
    }

    /// Fold diagonal (predecessor at `q-1`) contributions into `diag_buf`.
    ///
    /// Substitution cost is NOT applied here — it is applied once in `finalize`
    /// after every predecessor has contributed, so a single symbol comparison
    /// suffices regardless of predecessor count.
    ///
    /// ```text
    ///     pred band u   . X . . .     (src[pqi] at q-1)
    ///                     ╲
    ///                      ╲ diagonal, no cost yet
    ///                       ▼
    ///     diag_buf      . . Y . .     Y ← min(Y, X)
    ///
    ///     index shift:  pqi = q - 1 - src_qlo
    ///                   qi  = q     - dst_qlo
    /// ```
    fn accumulate_diag(
        diag_buf: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
    ) {
        let q_lo = dst_qlo.max(src_qlo + 1);
        let q_hi = (dst_qlo + w_dst - 1).min(src_qlo + w_src);
        if q_lo > q_hi {
            return;
        }
        for q in q_lo..=q_hi {
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            diag_buf[qi] = diag_buf[qi].min(src[pqi]);
        }
    }

    /// Close the band by applying substitution to `diag_buf` and sweeping
    /// insertion cost left → right through `data`.
    ///
    /// Prerequisites: every predecessor has already contributed via
    /// [`Self::accumulate_gap`] (into `data`) and [`Self::accumulate_diag`]
    /// (into `diag_buf`). This function is called exactly once per band.
    ///
    /// ```text
    ///   input:
    ///     data      = [ g0  g1  g2  g3 ]       (deletion contributions so far)
    ///     diag_buf  = [ x0  x1  x2  x3 ]       (best pred score at q-1)
    ///
    ///   phase 1 — apply subst to diag_buf (skip qi=0 when qlo==0, no query char):
    ///     diag_buf[qi] += subst(symbol, query[qlo+qi-1])
    ///
    ///   phase 2 — left→right sweep: insertion uses data[qi-1] just written
    ///
    ///     qi=0:   data[0] = min(g0, diag_buf[0])            (no left neighbour)
    ///     qi≥1:   data[qi] = min(g_qi, diag_buf[qi], data[qi-1] + ge)
    ///
    ///            data[0]──ge──▶data[1]──ge──▶data[2]──ge──▶data[3]
    ///                           ▲              ▲              ▲
    ///                           │  insertion chain propagates │
    /// ```
    fn finalize(
        data: &mut [u32],
        w: usize,
        qlo: usize,
        symbol: u8,
        query: &[u8],
        diag_buf: &mut [u32],
        costs: &Linear,
    ) {
        let ge = costs.gap_extend() as u32;
        let eq_cost = 0u32;
        let mm_cost = costs.mismatch() as u32;

        // Apply substitution to diag_buf in-place
        let subst_start = usize::from(qlo == 0);
        for qi in subst_start..w {
            let subst = if symbol == query[qlo + qi - 1] {
                eq_cost
            } else {
                mm_cost
            };
            diag_buf[qi] = sat(diag_buf[qi], subst);
        }

        // Sequential pass: M[qi] = min(deletion_accumulated, diagonal, insertion_from_left)
        // data[qi] already has the deletion cost from accumulate_gap.
        // qi=0: no insertion from left
        data[0] = data[0].min(diag_buf[0]);
        for qi in 1..w {
            let ins = sat(data[qi - 1], ge);
            data[qi] = data[qi].min(diag_buf[qi]).min(ins);
        }
    }

    #[inline]
    fn m_at(data: &[u32], _w: usize, qi: usize) -> u32 {
        data[qi]
    }

    #[inline]
    fn d_at(data: &[u32], _w: usize, qi: usize) -> u32 {
        data[qi]
    }

    #[inline]
    fn i_at(data: &[u32], _w: usize, qi: usize) -> u32 {
        data[qi]
    }

    fn best_state(col: &[u32]) -> (u32, u8) {
        (col[0], 0)
    }

    fn backtrace_op(
        col: &[u32],
        col_prev: Option<&[u32]>,
        _state: u8,
        costs: &Linear,
    ) -> BacktraceOp {
        // NOTE: The linear model folds deletions and diagonals into the same state,
        // making it impossible to distinguish them from the current cell's data alone.
        // When neither insertion nor the q==0 condition is detected, we return Diagonal
        // which causes backtrace_generic to decrement q and search predecessors at q-1.
        // For cells whose optimal path was a deletion (predecessor at same q), this
        // produces an incorrect alignment pair (match/mismatch instead of deletion).
        // Scores are always correct; only the alignment *representation* may have a wrong
        // op at ties. A fix would require backtrace_op to access predecessor data.
        match col_prev {
            None => BacktraceOp::Delete { next_state: 0 },
            Some(prev) => {
                let ge = costs.gap_extend() as u32;
                if col[0] == sat(prev[0], ge) {
                    BacktraceOp::Insert { next_state: 0 }
                } else {
                    BacktraceOp::Diagonal { next_state: 0 }
                }
            }
        }
    }

    fn del_transition_cost(col: &[u32], costs: &Linear) -> (u32, u8) {
        (sat(col[0], costs.gap_extend() as u32), 0)
    }
}

#[cfg(test)]
mod tests {
    use crate::align::{
        cost_models::linear::Linear,
        engine::{AlignOutput, band_doubling::BandDoublingEngineScalar, dp::CanonicalDP},
        traits::AlignmentEngine,
    };
    use crate::graph::{alignment::AddAlignment, poa::POAGraph};

    fn costs() -> Linear {
        Linear::new(1, 1)
    }

    fn linear_graph(seq: &[u8]) -> POAGraph<u32> {
        let mut g = POAGraph::new();
        g.add_alignment("s0", seq, None, &vec![1; seq.len()])
            .unwrap();
        g
    }

    fn run_canonical(graph: &POAGraph<u32>, query: &[u8]) -> AlignOutput<POAGraph<u32>> {
        CanonicalDP::new(costs()).align(graph, query).unwrap()
    }

    fn run_banded(graph: &POAGraph<u32>, query: &[u8]) -> AlignOutput<POAGraph<u32>> {
        BandDoublingEngineScalar::<Linear, u32>::new(costs())
            .align(graph, query)
            .unwrap()
    }

    #[test]
    fn linear_exact_match_scores_zero() {
        assert_eq!(run_canonical(&linear_graph(b"ACGT"), b"ACGT").score, 0);
    }

    #[test]
    fn linear_one_deletion_costs_one_extend() {
        // ACGT vs ACT: 1 deletion → 1 * gap_extend(1) = 1
        assert_eq!(run_canonical(&linear_graph(b"ACGT"), b"ACT").score, 1);
    }

    #[test]
    fn linear_gap_of_3_costs_three_extends() {
        // ACCCGT vs AGT: 3 deletions → 3 * 1 = 3
        assert_eq!(run_canonical(&linear_graph(b"ACCCGT"), b"AGT").score, 3);
    }

    #[test]
    fn linear_banded_matches_canonical() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGT", b"ACGT"),
            (b"ACGT", b"ACT"),
            (b"ACCCGT", b"AGT"),
            (b"TTTT", b"TTTTTTTT"),
        ];
        for (seq, query) in cases {
            let g = linear_graph(seq);
            let canon = run_canonical(&g, query);
            let banded = run_banded(&g, query);
            assert_eq!(
                banded.score,
                canon.score,
                "linear banded/canonical mismatch seq={} query={}",
                std::str::from_utf8(seq).unwrap(),
                std::str::from_utf8(query).unwrap()
            );
        }
    }
}
