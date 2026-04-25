//! # TwoPieceAffineKernel — 5 states, piecewise affine gap cost
//!
//! The two-piece gap cost for a gap of length `k`:
//!
//! ```text
//!     gap(k) = min( go1 + k·ge1 ,   piece 1 — typically "long gap, cheap extend"
//!                   go2 + k·ge2 )   piece 2 — typically "short gap, cheap open"
//! ```
//!
//! is realised by running two fully independent affine trackers in parallel.
//! Each piece has its own deletion and insertion state; the overall match state
//! takes the minimum across both tracks on every step.
//!
//! ## States
//!
//! | id | state | meaning                           |
//! |----|-------|-----------------------------------|
//! | 0  | M     | match / mismatch                  |
//! | 1  | D1    | deletion under piece 1 (go1, ge1) |
//! | 2  | D2    | deletion under piece 2 (go2, ge2) |
//! | 3  | I1    | insertion under piece 1           |
//! | 4  | I2    | insertion under piece 2           |
//!
//! ## Data layout (STATES = 5)
//!
//! ```text
//!   offset:  0    w     2w    3w    4w      5w
//!            │    │     │     │     │       │
//!            ▼    ▼     ▼     ▼     ▼
//!   data = [ M ][ D1 ][ D2 ][ I1 ][ I2 ]
//!           w u32 each, qi = 0..w-1
//! ```
//!
//! ## DP recurrences
//!
//! ```text
//!   D1[v][q] = min( M[u][q] + go1+ge1, I1[u][q] + go1+ge1, D1[u][q] + ge1 )
//!   D2[v][q] = min( M[u][q] + go2+ge2, I2[u][q] + go2+ge2, D2[u][q] + ge2 )
//!
//!   I1[v][q] = min( M[v][q-1] + go1+ge1, I1[v][q-1] + ge1 )
//!   I2[v][q] = min( M[v][q-1] + go2+ge2, I2[v][q-1] + ge2 )
//!
//!   M [v][q] = min( diag + subst(s, q.char),
//!                   D1[v][q], D2[v][q],
//!                   I1[v][q], I2[v][q] )
//! ```
//!
//! The two pieces NEVER cross-talk (no D1 ↔ D2, no I1 ↔ I2). Each track can
//! only be entered from M, so a gap is evaluated entirely under piece 1 or
//! entirely under piece 2, and the cheaper of the two wins at cell M.
//!
//! ## Per-band computation pipeline
//!
//! ```text
//!   1. accumulate_gap(u → v)    per pred u — updates D1 AND D2 columns
//!   2. accumulate_diag(u → v)   per pred u — diag_buf = min over all 5 states
//!   3. finalize                 once       — subst + (I1, I2, M) sweep
//! ```

use super::{BacktraceOp, DPKernel};
use crate::align::{
    cost_models::{AlignmentCostModel, two_piece::TwoPieceAffine},
    engine::dp::{INF, min3, sat},
};

/// 5-state two-piece affine kernel: states 0=M, 1=D1, 2=D2, 3=I1, 4=I2.
///
/// Gap cost of length k: min(go1 + k*ge1, go2 + k*ge2).
/// D1/I1 track the first (go1, ge1) piece; D2/I2 track the second (go2, ge2) piece.
pub struct TwoPieceAffineKernel;

impl DPKernel for TwoPieceAffineKernel {
    const STATES: usize = 5;
    const STATE_NAMES: &'static [&'static str] = &["M", "D1", "D2", "I1", "I2"];
    type Costs = TwoPieceAffine;

    /// Seed the start-sentinel band with both insertion tracks.
    ///
    /// `M[0] = 0` at the origin; D1/D2 stay at INF (no graph predecessor); I1
    /// and I2 at query offset `q ≥ 1` carry the cost of `q` query insertions
    /// under each piece.
    ///
    /// ```text
    ///   M  : [ 0    INF   INF   INF   ... ]
    ///   D1 : [ INF  INF   INF   INF   ... ]
    ///   D2 : [ INF  INF   INF   INF   ... ]
    ///   I1 : [ INF  go1+ge1  go1+2ge1  go1+3ge1  ... ]
    ///   I2 : [ INF  go2+ge2  go2+2ge2  go2+3ge2  ... ]
    /// ```
    fn init_start(data: &mut [u32], w: usize, qlo: usize, costs: &TwoPieceAffine) {
        let go1 = costs.gap_open() as u32;
        let ge1 = costs.gap_extend() as u32;
        let go2 = costs.gap_open2() as u32;
        let ge2 = costs.gap_extend2() as u32;
        // M[0] = 0
        data[0] = 0;
        // I1[q] = go1 + q*ge1, I2[q] = go2 + q*ge2 for q >= 1
        // State 3 = I1 at offset 3*w, state 4 = I2 at offset 4*w
        for qi in 1..w {
            let q = (qlo + qi) as u32;
            data[3 * w + qi] = sat(go1, q * ge1); // I1
            data[4 * w + qi] = sat(go2, q * ge2); // I2
        }
    }

    /// Fold D1 and D2 contributions from one predecessor band `src` into `dst`.
    ///
    /// Both tracks are updated independently — D1 pulls from M/I1/D1; D2 pulls
    /// from M/I2/D2. The two tracks do not cross-talk.
    ///
    /// ```text
    ///                    pred band u (same q)
    ///              ┌────┬─────┬─────┬─────┬─────┐
    ///     src[.] = │ M  │ D1  │ D2  │ I1  │ I2  │
    ///              └────┴─────┴─────┴─────┴─────┘
    ///
    ///
    ///             piece 1 ─────       ───── piece 2
    ///              M + go1+ge1         M + go2+ge2
    ///              I1+ go1+ge1         I2+ go2+ge2
    ///              D1+ ge1             D2+ ge2
    ///                  ▼                   ▼
    ///           dst[D1][qi]         dst[D2][qi]      (both updated, independently)
    /// ```
    fn accumulate_gap(
        dst: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &TwoPieceAffine,
    ) {
        let go1 = costs.gap_open() as u32;
        let ge1 = costs.gap_extend() as u32;
        let go_ge1 = go1 + ge1;
        let go2 = costs.gap_open2() as u32;
        let ge2 = costs.gap_extend2() as u32;
        let go_ge2 = go2 + ge2;

        let q_lo = dst_qlo.max(src_qlo);
        let q_hi = (dst_qlo + w_dst - 1).min(src_qlo + w_src - 1);
        if q_lo > q_hi {
            return;
        }

        for q in q_lo..=q_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            // Read predecessor states: 0=M, 1=D1, 2=D2, 3=I1, 4=I2
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];

            let m_or_i1 = pm.min(pi1); // M or I1 → D1 transition
            let m_or_i2 = pm.min(pi2); // M or I2 → D2 transition

            // D1[q] = min(M_pred[q] + go_ge1, I1_pred[q] + go_ge1, D1_pred[q] + ge1)
            dst[w_dst + qi] = dst[w_dst + qi].min(sat(m_or_i1, go_ge1)).min(sat(pd1, ge1));

            // D2[q] = min(M_pred[q] + go_ge2, I2_pred[q] + go_ge2, D2_pred[q] + ge2)
            dst[2 * w_dst + qi] = dst[2 * w_dst + qi]
                .min(sat(m_or_i2, go_ge2))
                .min(sat(pd2, ge2));
        }
    }

    /// Fold the diagonal-predecessor minimum across ALL five states into `diag_buf`.
    ///
    /// Substitution cost is applied once in `finalize`, so the buffer here
    /// carries a raw best-score across states at `q-1` in the predecessor.
    ///
    /// ```text
    ///     pred band u at q-1
    ///     ┌────┬─────┬─────┬─────┬─────┐
    ///     │ M  │ D1  │ D2  │ I1  │ I2  │
    ///     └────┴─────┴─────┴─────┴─────┘
    ///         \   |     |     |    /
    ///          ╲  │     │     │   ╱
    ///           ╲ └──▶ min ◀──┘  ╱
    ///            ╲      │       ╱
    ///             ╲     ▼      ╱
    ///              ▶ candidate ◀
    ///                    │
    ///                    ▼
    ///     diag_buf[qi] ← min(diag_buf[qi], candidate)
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
            // min over all 5 predecessor states
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            let best = pm.min(pd1).min(pd2).min(pi1).min(pi2);
            diag_buf[qi] = diag_buf[qi].min(best);
        }
    }

    /// Close the band: apply substitution to `diag_buf`, then compute I1, I2
    /// and M left → right using the already-filled D1 and D2 columns.
    ///
    /// ```text
    ///   phase 1: diag_buf[qi] += subst(symbol, query[qlo+qi-1])   (skip qi=0 when qlo=0)
    ///
    ///   phase 2: qi=0 boundary
    ///               I1[0] = I2[0] = INF
    ///               M[0]  = min( diag_buf[0], D1[0], D2[0] )
    ///
    ///   phase 3: for qi in 1..w  (left→right sweep, both pieces in parallel)
    ///
    ///       ┌────────── piece 1 ──────────┐ ┌────────── piece 2 ──────────┐
    ///       │                             │ │                             │
    ///       │ M[qi-1]  ──+go1+ge1──┐      │ │ M[qi-1]  ──+go2+ge2──┐      │
    ///       │                      ▼      │ │                      ▼      │
    ///       │ I1[qi-1] ───+ge1────▶ I1[qi]│ │ I2[qi-1] ───+ge2────▶ I2[qi]│
    ///       └─────────────────────────────┘ └─────────────────────────────┘
    ///
    ///                        combine:
    ///          M[qi] = min( diag_buf[qi],
    ///                       min(D1[qi], D2[qi]),     ← better of 2 deletion tracks
    ///                       min(I1[qi], I2[qi]) )    ← better of 2 insertion tracks
    /// ```
    ///
    /// Because `M[qi]` only looks at `D[qi]` / `I[qi]` (just computed above)
    /// and at `M[qi-1]` / `I[qi-1]` (written in the previous iteration), the
    /// sweep is single-pass.
    fn finalize(
        data: &mut [u32],
        w: usize,
        qlo: usize,
        symbol: u8,
        query: &[u8],
        diag_buf: &mut [u32],
        costs: &TwoPieceAffine,
    ) {
        let go1 = costs.gap_open() as u32;
        let ge1 = costs.gap_extend() as u32;
        let go_ge1 = go1 + ge1;
        let go2 = costs.gap_open2() as u32;
        let ge2 = costs.gap_extend2() as u32;
        let go_ge2 = go2 + ge2;
        let eq_cost = 0u32;
        let mm_cost = costs.mismatch() as u32;

        // Apply substitution cost in-place to diag_buf
        let subst_start = usize::from(qlo == 0);
        for qi in subst_start..w {
            let subst = if symbol == query[qlo + qi - 1] {
                eq_cost
            } else {
                mm_cost
            };
            diag_buf[qi] = sat(diag_buf[qi], subst);
        }

        // Boundary: I1[0] = I2[0] = INF, M[0] = min(diag, D1[0], D2[0])
        data[3 * w] = INF; // I1[0]
        data[4 * w] = INF; // I2[0]
        let d_best_0 = data[w].min(data[2 * w]);
        data[0] = diag_buf[0].min(d_best_0);

        for qi in 1..w {
            let m_prev = data[qi - 1];
            let i1_prev = data[3 * w + qi - 1];
            let i2_prev = data[4 * w + qi - 1];
            // I1[qi] = min(M[qi-1] + go_ge1, I1[qi-1] + ge1)
            data[3 * w + qi] = sat(m_prev, go_ge1).min(sat(i1_prev, ge1));
            // I2[qi] = min(M[qi-1] + go_ge2, I2[qi-1] + ge2)
            data[4 * w + qi] = sat(m_prev, go_ge2).min(sat(i2_prev, ge2));
            let d_best = data[w + qi].min(data[2 * w + qi]);
            let i_best = data[3 * w + qi].min(data[4 * w + qi]);
            data[qi] = min3(diag_buf[qi], d_best, i_best);
        }
    }

    fn m_at(data: &[u32], _w: usize, qi: usize) -> u32 {
        data[qi]
    }
    fn d_at(data: &[u32], w: usize, qi: usize) -> u32 {
        data[w + qi].min(data[2 * w + qi])
    }
    fn i_at(data: &[u32], w: usize, qi: usize) -> u32 {
        data[3 * w + qi].min(data[4 * w + qi])
    }

    fn best_state(col: &[u32]) -> (u32, u8) {
        let m = col[0];
        let d1 = col[1];
        let d2 = col[2];
        let i1 = col[3];
        let i2 = col[4];
        let best = m.min(d1).min(d2).min(i1).min(i2);
        let state = if m == best {
            0
        } else if d1 == best {
            1
        } else if d2 == best {
            2
        } else if i1 == best {
            3
        } else {
            4
        };
        (best, state)
    }

    fn backtrace_op(
        col: &[u32],
        col_prev: Option<&[u32]>,
        state: u8,
        costs: &TwoPieceAffine,
    ) -> BacktraceOp {
        let go1 = costs.gap_open() as u32;
        let ge1 = costs.gap_extend() as u32;
        let go_ge1 = go1 + ge1;
        let go2 = costs.gap_open2() as u32;
        let ge2 = costs.gap_extend2() as u32;
        let go_ge2 = go2 + ge2;

        match state {
            0 => {
                // M can arrive from diagonal, or from closing a D1/D2/I1/I2 gap
                // at the same cell. Closing emits SwitchState (no move) so the
                // next iteration walks the gap backwards via Delete/Insert.
                // At q=0 only D closures are reachable.
                let d_best = col[1].min(col[2]);
                let i_best = col[3].min(col[4]);
                if col_prev.is_none() {
                    let d_state = if col[0] == col[1] || col[1] <= col[2] {
                        1
                    } else {
                        2
                    };
                    BacktraceOp::SwitchState {
                        next_state: d_state,
                    }
                } else if col[0] == d_best {
                    let d_state = if col[1] <= col[2] { 1 } else { 2 };
                    BacktraceOp::SwitchState {
                        next_state: d_state,
                    }
                } else if col[0] == i_best {
                    let i_state = if col[3] <= col[4] { 3 } else { 4 };
                    BacktraceOp::SwitchState {
                        next_state: i_state,
                    }
                } else {
                    BacktraceOp::Diagonal { next_state: 0 }
                }
            }
            1 => BacktraceOp::Delete { next_state: 0xff }, // D1: resolve pred in caller
            2 => BacktraceOp::Delete { next_state: 0xff }, // D2
            3 => {
                // I1
                let prev = col_prev.expect("q=0 in I1 state");
                let m_prev = prev[0];
                let i1_prev = prev[3];
                let next = if sat(m_prev, go_ge1) <= sat(i1_prev, ge1) {
                    0
                } else {
                    3
                };
                BacktraceOp::Insert { next_state: next }
            }
            4 => {
                // I2
                let prev = col_prev.expect("q=0 in I2 state");
                let m_prev = prev[0];
                let i2_prev = prev[4];
                let next = if sat(m_prev, go_ge2) <= sat(i2_prev, ge2) {
                    0
                } else {
                    4
                };
                BacktraceOp::Insert { next_state: next }
            }
            _ => unreachable!(),
        }
    }

    fn del_transition_cost(col: &[u32], costs: &TwoPieceAffine) -> (u32, u8) {
        let go1 = costs.gap_open() as u32;
        let ge1 = costs.gap_extend() as u32;
        let go_ge1 = go1 + ge1;
        let go2 = costs.gap_open2() as u32;
        let ge2 = costs.gap_extend2() as u32;
        let go_ge2 = go2 + ge2;

        let pm = col[0];
        let pd1 = col[1];
        let pd2 = col[2];
        let pi1 = col[3];
        let pi2 = col[4];

        let cost1_open = sat(pm.min(pi1), go_ge1); // open D1 from M or I1
        let cost1_ext = sat(pd1, ge1); // extend D1
        let cost1 = cost1_open.min(cost1_ext);

        let cost2_open = sat(pm.min(pi2), go_ge2); // open D2 from M or I2
        let cost2_ext = sat(pd2, ge2); // extend D2
        let cost2 = cost2_open.min(cost2_ext);

        if cost1 <= cost2 {
            (cost1, 1)
        } else {
            (cost2, 2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{
        cost_models::two_piece::TwoPieceAffine,
        engine::{AlignOutput, band_doubling::BandDoublingEngineScalar, dp::CanonicalDP},
        traits::AlignmentEngine,
    };
    use crate::graph::{alignment::AddAlignment, poa::POAGraph};

    /// Two-piece affine: (go1=2, ge1=1) and (go2=0, ge2=3).
    /// For a gap of length k, cost = min(2+k, 3k).
    /// k=1: min(3, 3)=3, k=2: min(4, 6)=4, k=3: min(5, 9)=5
    fn costs() -> TwoPieceAffine {
        TwoPieceAffine::new(1, 2, 1, 0, 3)
    }

    fn linear_graph(seq: &[u8]) -> POAGraph<u32> {
        let mut g = POAGraph::new();
        g.add_alignment("s0", seq, None, &vec![1; seq.len()])
            .unwrap();
        g
    }

    fn run_canonical(graph: &POAGraph<u32>, query: &[u8]) -> AlignOutput<POAGraph<u32>> {
        let engine: CanonicalDP<TwoPieceAffine, POAGraph<u32>> = CanonicalDP::new(costs());
        engine.align(graph, query).unwrap()
    }

    #[test]
    fn two_piece_exact_match() {
        let g = linear_graph(b"ACGT");
        assert_eq!(run_canonical(&g, b"ACGT").score, 0);
    }

    #[test]
    fn two_piece_gap_of_1() {
        // ACGT vs ACT: gap length 1 → min(2+1, 3*1) = 3
        let g = linear_graph(b"ACGT");
        assert_eq!(run_canonical(&g, b"ACT").score, 3);
    }

    #[test]
    fn two_piece_gap_of_2() {
        // ACGT vs AT: gap length 2 (delete CG) → min(2+2, 3*2) = min(4, 6) = 4
        let g = linear_graph(b"ACGT");
        assert_eq!(run_canonical(&g, b"AT").score, 4);
    }

    #[test]
    fn two_piece_banded_matches_canonical() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGT", b"ACGT"),
            (b"ACGT", b"ACT"),
            (b"ACGT", b"AT"),
            (b"TTTT", b"TTTTTTTT"),
        ];
        for (seq, query) in cases {
            let g = linear_graph(seq);
            let canon = run_canonical(&g, query);
            let banded = BandDoublingEngineScalar::<TwoPieceAffine, u32>::new(costs())
                .align(&g, query)
                .unwrap();
            assert_eq!(
                banded.score,
                canon.score,
                "two-piece banded/canonical mismatch seq={} query={}",
                std::str::from_utf8(seq).unwrap(),
                std::str::from_utf8(query).unwrap()
            );
        }
    }
}
