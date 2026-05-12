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
//!   step 1  accumulate_pred(u → v)    per pred u — fold (b) into data and
//!                                                   (a) into diag_buf in one pass
//!   step 2  finalize                  once       — apply subst, then left→right (c) sweep
//! ```
//!
//! ### step 1: `accumulate_pred` — one call per predecessor band `u`
//!
//! Single pass over `src[u]` updates both reductions:
//!
//! ```text
//!     pred band u   . X . . .          (state 0 at query position q)
//!                     │
//!                     │  + ge          (b) deletion: same q, v ← u
//!                     ▼
//!     curr band v   . Y . . .          Y = min(Y, X + ge)
//!
//!     pred band u   X . . . .          (at q-1)
//!                    ╲
//!                     ▼ diagonal (no subst yet)   (a) diagonal: q ← q-1, v ← u
//!     diag_buf      . Y . .            Y = min(Y, X)
//! ```
//!
//! ### step 2: `finalize` — apply subst, then left→right insertion sweep
//!
//! ```text
//!     after accumulate_pred:  data     = [ g0  g1  g2  g3 ]   (deletion cost per cell)
//!     after substitution:     diag_buf = [ d0  d1  d2  d3 ]   (diag + subst)
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

    /// Fold one predecessor band into both `dst` (deletion cost per cell) and
    /// `diag_buf` (best predecessor score at `q-1`), in a single pass over `src`.
    ///
    /// The loop is split into three sections so each iteration body is
    /// straight-line: a leading gap-only cell (at most one), the bulk where both
    /// updates apply, and a trailing diag-only cell (at most one). Iteration
    /// covers `q ∈ [max(dst_qlo, src_qlo), min(dst_qhi, src_qhi)+1]` — the
    /// query-range overlap between the two bands, extended by one for the
    /// diagonal lookup at `q-1`.
    ///
    /// ```text
    ///     gap update at q    pred src[q]      → dst[q]  += ge
    ///     diag update at q   pred src[q-1]    → diag_buf[q] = min(., src[q-1])
    /// ```
    ///
    /// Substitution cost is NOT applied here — `finalize` applies it once after
    /// every predecessor has contributed.
    fn accumulate_pred(
        dst: &mut [u32],
        diag_buf: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Linear,
    ) {
        let ge = costs.gap_extend() as u32;

        let gap_q_lo = dst_qlo.max(src_qlo);
        let gap_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src);
        let diag_q_lo = dst_qlo.max(src_qlo + 1);
        let diag_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src + 1);

        // Section A: gap-only prefix (0 or 1 cell).
        if gap_q_lo < diag_q_lo && gap_q_lo < gap_q_hi {
            let q = gap_q_lo;
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            dst[qi] = dst[qi].min(sat(src[pqi], ge));
        }

        // Section B: both updates apply. Split into two passes so each loop
        // is a clean linear stream (better autovectorization).
        let both_lo = diag_q_lo;
        let both_hi = gap_q_hi.min(diag_q_hi);

        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            dst[qi] = dst[qi].min(sat(src[pqi], ge));
        }
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            diag_buf[qi] = diag_buf[qi].min(src[pqi]);
        }

        // Section C: diag-only tail (0 or 1 cell).
        if gap_q_hi < diag_q_hi {
            let q = gap_q_hi;
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            diag_buf[qi] = diag_buf[qi].min(src[pqi]);
        }
    }

    /// Like [`Self::accumulate_pred`] but writes destinations with `=` instead
    /// of `min`-folding. Caller-side preconditions match
    /// [`AffineKernel::init_pred`]: `dst` is pre-initialized to INF, and the
    /// destination band's range matches the predecessor's shifted range
    /// (true for single-pred bands).
    fn init_pred(
        dst: &mut [u32],
        diag_buf: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Linear,
    ) {
        let ge = costs.gap_extend() as u32;

        let gap_q_lo = dst_qlo.max(src_qlo);
        let gap_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src);
        let diag_q_lo = dst_qlo.max(src_qlo + 1);
        let diag_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src + 1);

        // Section A: gap-only prefix.
        if gap_q_lo < diag_q_lo && gap_q_lo < gap_q_hi {
            let q = gap_q_lo;
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            dst[qi] = sat(src[pqi], ge);
        }

        // Section B: both updates apply (write directly).
        let both_lo = diag_q_lo;
        let both_hi = gap_q_hi.min(diag_q_hi);
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            dst[qi] = sat(src[pqi], ge);
        }
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            diag_buf[qi] = src[pqi];
        }

        // Section C: diag-only tail.
        if gap_q_hi < diag_q_hi {
            let q = gap_q_hi;
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            diag_buf[qi] = src[pqi];
        }
    }

    /// Close the band by applying substitution to `diag_buf` and sweeping
    /// insertion cost left → right through `data`.
    ///
    /// Prerequisites: every predecessor has already contributed via
    /// [`Self::accumulate_pred`] (folding deletion cost into `data` and the
    /// diagonal minimum into `diag_buf`). This function is called exactly once
    /// per band.
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
        let mm_cost = costs.mismatch() as u32;

        // Boundary qi=0: no insertion from left, no query char to compare at
        // qlo == 0 (only reachable for the start sentinel's column).
        let diag0 = if qlo == 0 {
            diag_buf[0]
        } else {
            let mismatch = (symbol != query[qlo - 1]) as u32;
            sat(diag_buf[0], mm_cost & mismatch.wrapping_neg())
        };
        data[0] = data[0].min(diag0);

        // Sequential left→right sweep with per-cell substitution mask folded
        // in. The RAW dependency on data[qi-1] makes this loop scalar; folding
        // saves one streaming pass over diag_buf.
        for qi in 1..w {
            let q_byte = query[qlo + qi - 1];
            let mismatch = (symbol != q_byte) as u32;
            let diag = sat(diag_buf[qi], mm_cost & mismatch.wrapping_neg());
            let ins = sat(data[qi - 1], ge);
            data[qi] = data[qi].min(diag).min(ins);
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
