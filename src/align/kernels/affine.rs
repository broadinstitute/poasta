//! # AffineKernel — 3 states, affine gap cost (`go + k·ge`)
//!
//! Tracks three distinct states per cell:
//!
//! - **M** — match/mismatch (diagonal arrival)
//! - **D** — deletion in query = graph-only step (gap open or extend)
//! - **I** — insertion in query = query-only step (gap open or extend)
//!
//! Splitting D and I out of M is what makes affine costs possible: the "open"
//! penalty `go+ge` is only paid on entering a gap, and subsequent extensions
//! inside D or I only pay `ge`.
//!
//! ## Data layout (STATES = 3, w = qhi - qlo + 1)
//!
//! Scores are packed by state then query index:
//!
//! ```text
//!   offset:    0           w          2w         3w
//!              │           │          │          │
//!              ▼           ▼          ▼
//!   data = [ M[0..w) ][ D[0..w) ][ I[0..w) ]
//!           \───────/  \───────/  \───────/
//!             w u32      w u32      w u32
//! ```
//!
//! ## DP recurrences
//!
//! ```text
//!   D[v][q] = min( M[u][q] + go + ge,     gap open from M (in pred u, same q)
//!                  I[u][q] + go + ge,     gap open from I (in pred u, same q)
//!                  D[u][q] + ge )          extend D (from pred u, same q)
//!
//!   I[v][q] = min( M[v][q-1] + go + ge,   gap open from M (this band, q-1)
//!                  I[v][q-1] + ge )       extend I (this band, q-1)
//!
//!   M[v][q] = min( diag + subst(s, q.char),  diagonal (arrives from any state)
//!                  D[v][q],                  close a deletion into a match cell
//!                  I[v][q] )                 close an insertion into a match cell
//! ```
//!
//! Crucially `D ↛ I` and `I ↛ D` directly — every D→I (or I→D) switch must pass
//! through M, which forces a new gap-open penalty. This is what prevents
//! double-charging when two gaps of different orientations abut.
//!
//! ## Per-band computation pipeline
//!
//! ```text
//!   1. accumulate_gap(u → v)    (per predecessor u)  — fills D column
//!   2. accumulate_diag(u → v)   (per predecessor u)  — fills diag_buf
//!   3. finalize                  (once)              — subst + I/M sweep
//! ```

use super::{BacktraceOp, DPKernel};
use crate::align::{
    cost_models::{affine::Affine, AlignmentCostModel},
    engine::dp::{min3, sat, INF},
};

/// 3-state affine gap kernel: states 0=M, 1=D, 2=I.
/// Data layout: `data[s * w + qi]`.
pub struct AffineKernel;

impl DPKernel for AffineKernel {
    const STATES: usize = 3;
    type Costs = Affine;

    /// Seed the start-sentinel band.
    ///
    /// `M[0] = 0` anchors the origin of the alignment; `D` is left at INF
    /// because there is no graph position before the start sentinel; `I[q]` for
    /// `q ≥ 1` is the cost of `q` query insertions from the origin: one
    /// gap-open + q extensions.
    ///
    /// ```text
    ///   state 0 (M):  [ 0   INF INF INF ...   ]
    ///   state 1 (D):  [ INF INF INF INF ...   ]       (untouched here)
    ///   state 2 (I):  [ INF go+ge go+2ge go+3ge ...]
    /// ```
    fn init_start(data: &mut [u32], w: usize, qlo: usize, costs: &Affine) {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        data[0] = 0; // M at qi=0
        for qi in 1..w {
            let q = qlo + qi;
            data[2 * w + qi] = sat(go, q as u32 * ge); // I[qi]
        }
    }

    /// Fold D-state contributions from one predecessor band `src` into `dst`.
    ///
    /// Called once per predecessor; successive calls shrink `dst[D]` via `min`.
    ///
    /// ```text
    ///                     pred band u (same q)
    ///                 ┌─────────┬──────────┬──────────┐
    ///     src[pqi] =  │  M_u    │   D_u    │   I_u    │
    ///                 └─────────┴──────────┴──────────┘
    ///                      │         │          │
    ///                      │ +go+ge  │ +ge      │ +go+ge
    ///                      │         │          │
    ///                      └──────▶  ▼  ◀───────┘
    ///                              min ⇒ candidate for D_v[q]
    ///                                 │
    ///                                 ▼
    ///                    dst[D][qi] ← min( dst[D][qi], candidate )
    ///
    ///     indices:  pqi = q - src_qlo      qi = q - dst_qlo
    /// ```
    fn accumulate_gap(
        dst: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Affine,
    ) {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        let go_ge = go + ge;

        let q_lo = dst_qlo.max(src_qlo);
        let q_hi = (dst_qlo + w_dst - 1).min(src_qlo + w_src - 1);
        if q_lo > q_hi {
            return;
        }
        for q in q_lo..=q_hi {
            let qi = q - dst_qlo;
            let pqi = q - src_qlo;
            let pm = src[pqi];
            let pd = src[w_src + pqi];
            let pi = src[2 * w_src + pqi];
            dst[w_dst + qi] = dst[w_dst + qi].min(sat(pm.min(pi), go_ge)).min(sat(pd, ge));
        }
    }

    /// Fold the diagonal-predecessor minimum into `diag_buf`.
    ///
    /// Substitution is NOT added here — it is applied once in `finalize` after
    /// every predecessor has contributed (so the symbol compare runs once per
    /// cell, not once per predecessor).
    ///
    /// ```text
    ///     pred band u at q-1:
    ///                                ┌──────┬──────┬──────┐
    ///                                │ M_u  │ D_u  │ I_u  │
    ///                                └──────┴──────┴──────┘
    ///                                     \   |   /
    ///                                      ╲  │  ╱   min3(M_u, D_u, I_u)
    ///                                       ▼ ▼ ▼
    ///     diag_buf[qi] ← min(diag_buf[qi], min3(·))
    ///
    ///     (arrives diagonally: graph v ← u, query q ← q-1)
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
            let pm = src[pqi];
            let pd = src[w_src + pqi];
            let pi = src[2 * w_src + pqi];
            diag_buf[qi] = diag_buf[qi].min(min3(pm, pd, pi));
        }
    }

    /// Close the band: apply substitution to `diag_buf`, then compute `I` and
    /// `M` left → right using the already-filled `D` column.
    ///
    /// The sweep is sequential because `I[qi]` depends on `M[qi-1]` and
    /// `I[qi-1]` that we have just written in this loop body.
    ///
    /// ```text
    ///   phase 1: diag_buf[qi] += subst(symbol, query[qlo+qi-1])   (skip qi=0 when qlo=0)
    ///
    ///   phase 2: qi=0 boundary
    ///               I[0] = INF                     (no left neighbour in band)
    ///               M[0] = min(diag_buf[0], D[0])
    ///
    ///   phase 3: for qi in 1..w  (left→right sweep)
    ///
    ///       ┌─────────────── within this band ───────────────┐
    ///       │                                                │
    ///       │   M[qi-1]  ──+go+ge──┐                         │
    ///       │                      ▼                         │
    ///       │   I[qi-1]  ───+ge───▶ I[qi]                    │
    ///       │                                                │
    ///       │   D[qi]  (already filled by accumulate_gap)    │
    ///       │                                                │
    ///       │                          \       │             │
    ///       │                           ▼      ▼             │
    ///       │     min3( diag_buf[qi], D[qi], I[qi] ) → M[qi] │
    ///       └────────────────────────────────────────────────┘
    /// ```
    fn finalize(
        data: &mut [u32],
        w: usize,
        qlo: usize,
        symbol: u8,
        query: &[u8],
        diag_buf: &mut [u32],
        costs: &Affine,
    ) {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        let go_ge = go + ge;
        let eq_cost = costs.equal() as u32;
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

        // I[0] = INF (no left neighbor within band at qi=0)
        data[2 * w] = INF;
        // M[0] = min(diag_buf[0], D[0])
        data[0] = diag_buf[0].min(data[w]);

        for qi in 1..w {
            let m_prev = data[qi - 1];
            let i_prev = data[2 * w + qi - 1];
            data[2 * w + qi] = sat(m_prev, go_ge).min(sat(i_prev, ge));
            data[qi] = min3(diag_buf[qi], data[w + qi], data[2 * w + qi]);
        }
    }

    #[inline]
    fn m_at(data: &[u32], _w: usize, qi: usize) -> u32 {
        data[qi]
    }

    #[inline]
    fn d_at(data: &[u32], w: usize, qi: usize) -> u32 {
        data[w + qi]
    }

    #[inline]
    fn i_at(data: &[u32], w: usize, qi: usize) -> u32 {
        data[2 * w + qi]
    }

    fn best_state(col: &[u32]) -> (u32, u8) {
        let m = col[0];
        let d = col[1];
        let i = col[2];
        let best = min3(m, d, i);
        let state = if m <= d && m <= i {
            0
        } else if d <= i {
            1
        } else {
            2
        };
        (best, state)
    }

    fn del_transition_cost(col: &[u32], costs: &Affine) -> (u32, u8) {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        let go_ge = go + ge;
        let pm = col[0];
        let pd = col[1];
        let pi = col[2];
        let open_cost = sat(pm.min(pi), go_ge);
        let extend_cost = sat(pd, ge);
        if extend_cost <= open_cost {
            (extend_cost, 1) // stay in D state (extend)
        } else if pm <= pi {
            (open_cost, 0) // gap opened from M state
        } else {
            (open_cost, 2) // gap opened from I state
        }
    }

    fn backtrace_op(
        col: &[u32],
        col_prev: Option<&[u32]>,
        state: u8,
        costs: &Affine,
    ) -> BacktraceOp {
        let go = costs.gap_open() as u32;
        let ge = costs.gap_extend() as u32;
        let go_ge = go + ge;

        match state {
            2 => {
                let prev = col_prev.expect("q=0 in I state");
                let m_prev = prev[0];
                let i_prev = prev[2];
                let next = if sat(m_prev, go_ge) <= sat(i_prev, ge) {
                    0
                } else {
                    2
                };
                BacktraceOp::Insert { next_state: next }
            }
            0 => {
                if col_prev.is_none() {
                    BacktraceOp::Diagonal { next_state: 1 }
                } else {
                    let m_val = col[0];
                    let d_val = col[1];
                    let i_val = col[2];
                    if m_val == d_val && d_val <= i_val {
                        BacktraceOp::Diagonal { next_state: 1 }
                    } else if m_val == i_val && i_val < d_val {
                        BacktraceOp::Diagonal { next_state: 2 }
                    } else {
                        BacktraceOp::Diagonal { next_state: 0 }
                    }
                }
            }
            1 => BacktraceOp::Delete { next_state: 0xff },
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cost_models::affine::Affine;
    use crate::align::engine::dp::INF;

    #[test]
    fn affine_kernel_m_at_reads_state0() {
        // data layout: [M_0, M_1, D_0, D_1, I_0, I_1]  for w=2
        let w = 2;
        let mut data = vec![INF; AffineKernel::STATES * w];
        data[0] = 7; // M at qi=0
        data[1] = 8; // M at qi=1
        assert_eq!(AffineKernel::m_at(&data, w, 0), 7);
        assert_eq!(AffineKernel::m_at(&data, w, 1), 8);
    }

    #[test]
    fn affine_kernel_init_start_sets_m_and_i() {
        // w=4, qlo=0: M[0]=0, I[1]=go+ge, I[2]=go+2ge, I[3]=go+3ge; D all INF
        let costs = Affine::new(0, 1, 2, 1);
        let w = 4;
        let mut data = vec![INF; AffineKernel::STATES * w];
        AffineKernel::init_start(&mut data, w, 0, &costs);
        // M state (state 0)
        assert_eq!(data[0], 0, "M[0] should be 0");
        assert_eq!(data[1], INF, "M[1] should be INF");
        // D state (state 1) — all INF
        for qi in 0..w {
            assert_eq!(data[w + qi], INF, "D[{qi}] should be INF");
        }
        // I state (state 2)
        assert_eq!(data[2 * w], INF, "I[0] should be INF");
        assert_eq!(data[2 * w + 1], 3, "I[1] = go(2)+ge(1) = 3");
        assert_eq!(data[2 * w + 2], 4, "I[2] = go(2)+2*ge(1) = 4");
        assert_eq!(data[2 * w + 3], 5, "I[3] = go(2)+3*ge(1) = 5");
    }

    #[test]
    fn affine_kernel_accumulate_gap_basic() {
        // Predecessor band: w_src=3, src_qlo=0, M=[0,INF,INF], D=[INF,INF,INF], I=[INF,INF,INF]
        // Query band: w_dst=3, dst_qlo=0
        // D[q] = min(sat(M[q].min(I[q]), go_ge), sat(D[q], ge))
        //       = sat(0, 3) = 3 at qi=0
        let costs = Affine::new(0, 1, 2, 1); // go=2, ge=1, go_ge=3
        let w = 3;
        let mut src = vec![INF; AffineKernel::STATES * w];
        src[0] = 0; // M[0] = 0
        let mut dst = vec![INF; AffineKernel::STATES * w];
        AffineKernel::accumulate_gap(&mut dst, w, 0, &src, w, 0, &costs);
        assert_eq!(dst[w], 3, "D[qi=0] = sat(M=0, go_ge=3) = 3");
        assert_eq!(dst[w + 1], INF, "D[qi=1] should remain INF (pred M[1]=INF)");
    }

    #[test]
    fn affine_kernel_finalize_exact_match() {
        // Single node, w=1, qlo=1, symbol='A', query=b"A"
        // diag_buf=[0] (predecessor score at q=0 was 0)
        // Expect M[0] = 0+0 = 0 (exact match, equal cost = 0)
        let costs = Affine::new(0, 1, 2, 1);
        let w = 1;
        let mut data = vec![INF; AffineKernel::STATES * w];
        data[w] = INF; // D[0] = INF
        let mut diag_buf = [0u32]; // best predecessor at diagonal position
        AffineKernel::finalize(&mut data, w, 1, b'A', b"A", &mut diag_buf, &costs);
        assert_eq!(data[0], 0, "M[0] should be 0 for exact match");
    }

    #[test]
    fn affine_kernel_finalize_mismatch() {
        // Same as above but query='T' != symbol='A': mismatch cost = 1
        let costs = Affine::new(0, 1, 2, 1);
        let w = 1;
        let mut data = vec![INF; AffineKernel::STATES * w];
        data[w] = INF;
        let mut diag_buf = [0u32];
        AffineKernel::finalize(&mut data, w, 1, b'A', b"T", &mut diag_buf, &costs);
        assert_eq!(data[0], 1, "M[0] should be 1 for mismatch");
    }
}
