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
//!   1. accumulate_pred(u → v)   per pred u — single pass over src:
//!                                              updates D1 AND D2 columns AND
//!                                              folds min over all 5 states into diag_buf
//!   2. finalize                 once       — subst + (I1, I2, M) sweep
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

    /// Fold one predecessor band into `dst` (D1 and D2 gap states) and
    /// `diag_buf` (min across all 5 states at `q-1`), in a single pass over `src`.
    ///
    /// The two deletion tracks are updated independently — D1 pulls from
    /// M/I1/D1 with piece-1 costs, D2 pulls from M/I2/D2 with piece-2 costs;
    /// no D1↔D2 cross-talk. Substitution is NOT added here — `finalize`
    /// applies it once after every predecessor has contributed.
    ///
    /// **Gap update** (same `q`, graph step `v ← u`):
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
    ///
    /// **Diag update** (graph step `v ← u`, query `q ← q-1`):
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
    ///
    /// The implementation splits the loop into three sections (gap-only
    /// prefix, both-updates bulk, diag-only tail) so each iteration body is
    /// straight-line.
    fn accumulate_pred(
        dst: &mut [u32],
        diag_buf: &mut [u32],
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

        let gap_q_lo = dst_qlo.max(src_qlo);
        let gap_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src);
        let diag_q_lo = dst_qlo.max(src_qlo + 1);
        let diag_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src + 1);

        // Section A: gap-only prefix (0 or 1 cell).
        if gap_q_lo < diag_q_lo && gap_q_lo < gap_q_hi {
            let q = gap_q_lo;
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            dst[w_dst + qi] = dst[w_dst + qi]
                .min(sat(pm.min(pi1), go_ge1))
                .min(sat(pd1, ge1));
            dst[2 * w_dst + qi] = dst[2 * w_dst + qi]
                .min(sat(pm.min(pi2), go_ge2))
                .min(sat(pd2, ge2));
        }

        // Section B: both updates apply. Split into independent passes so
        // each loop is a clean SoA stream over `src` (better autovectorization
        // across the 5 contiguous state arrays).
        let both_lo = diag_q_lo;
        let both_hi = gap_q_hi.min(diag_q_hi);

        // Pass 1: D1 / D2 gap fold (writes dst[D1] and dst[D2]).
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            dst[w_dst + qi] = dst[w_dst + qi]
                .min(sat(pm.min(pi1), go_ge1))
                .min(sat(pd1, ge1));
            dst[2 * w_dst + qi] = dst[2 * w_dst + qi]
                .min(sat(pm.min(pi2), go_ge2))
                .min(sat(pd2, ge2));
        }
        // Pass 2: diag fold (writes diag_buf[qi]).
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            diag_buf[qi] = diag_buf[qi].min(pm.min(pd1).min(pd2).min(pi1).min(pi2));
        }

        // Section C: diag-only tail (0 or 1 cell).
        if gap_q_hi < diag_q_hi {
            let q = gap_q_hi;
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            let best = pm.min(pd1).min(pd2).min(pi1).min(pi2);
            diag_buf[qi] = diag_buf[qi].min(best);
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
        costs: &TwoPieceAffine,
    ) {
        let go1 = costs.gap_open() as u32;
        let ge1 = costs.gap_extend() as u32;
        let go_ge1 = go1 + ge1;
        let go2 = costs.gap_open2() as u32;
        let ge2 = costs.gap_extend2() as u32;
        let go_ge2 = go2 + ge2;

        let gap_q_lo = dst_qlo.max(src_qlo);
        let gap_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src);
        let diag_q_lo = dst_qlo.max(src_qlo + 1);
        let diag_q_hi = (dst_qlo + w_dst).min(src_qlo + w_src + 1);

        // Section A: gap-only prefix.
        if gap_q_lo < diag_q_lo && gap_q_lo < gap_q_hi {
            let q = gap_q_lo;
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            dst[w_dst + qi] = sat(pm.min(pi1), go_ge1).min(sat(pd1, ge1));
            dst[2 * w_dst + qi] = sat(pm.min(pi2), go_ge2).min(sat(pd2, ge2));
        }

        // Section B: both updates apply (write directly).
        let both_lo = diag_q_lo;
        let both_hi = gap_q_hi.min(diag_q_hi);
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            dst[w_dst + qi] = sat(pm.min(pi1), go_ge1).min(sat(pd1, ge1));
            dst[2 * w_dst + qi] = sat(pm.min(pi2), go_ge2).min(sat(pd2, ge2));
        }
        for q in both_lo..both_hi {
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            diag_buf[qi] = pm.min(pd1).min(pd2).min(pi1).min(pi2);
        }

        // Section C: diag-only tail.
        if gap_q_hi < diag_q_hi {
            let q = gap_q_hi;
            let qi = q - dst_qlo;
            let pqi = q - 1 - src_qlo;
            let pm = src[pqi];
            let pd1 = src[w_src + pqi];
            let pd2 = src[2 * w_src + pqi];
            let pi1 = src[3 * w_src + pqi];
            let pi2 = src[4 * w_src + pqi];
            diag_buf[qi] = pm.min(pd1).min(pd2).min(pi1).min(pi2);
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
        let mm_cost = costs.mismatch() as u32;

        // Boundary qi=0: I1[0] = I2[0] = INF, M[0] = min(diag + subst, D1[0], D2[0]).
        data[3 * w] = INF; // I1[0]
        data[4 * w] = INF; // I2[0]
        let diag0 = if qlo == 0 {
            diag_buf[0]
        } else {
            let mismatch = (symbol != query[qlo - 1]) as u32;
            sat(diag_buf[0], mm_cost & mismatch.wrapping_neg())
        };
        let d_best_0 = data[w].min(data[2 * w]);
        data[0] = diag0.min(d_best_0);

        // Sequential left→right sweep with per-cell substitution mask folded
        // in. The loop is scalar regardless (RAW deps on M[qi-1], I1[qi-1],
        // I2[qi-1]); folding saves one streaming pass over diag_buf.
        for qi in 1..w {
            let q_byte = query[qlo + qi - 1];
            let mismatch = (symbol != q_byte) as u32;
            let diag = sat(diag_buf[qi], mm_cost & mismatch.wrapping_neg());

            let m_prev = data[qi - 1];
            let i1_prev = data[3 * w + qi - 1];
            let i2_prev = data[4 * w + qi - 1];
            // I1[qi] = min(M[qi-1] + go_ge1, I1[qi-1] + ge1)
            data[3 * w + qi] = sat(m_prev, go_ge1).min(sat(i1_prev, ge1));
            // I2[qi] = min(M[qi-1] + go_ge2, I2[qi-1] + ge2)
            data[4 * w + qi] = sat(m_prev, go_ge2).min(sat(i2_prev, ge2));
            let d_best = data[w + qi].min(data[2 * w + qi]);
            let i_best = data[3 * w + qi].min(data[4 * w + qi]);
            data[qi] = min3(diag, d_best, i_best);
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
