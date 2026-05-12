pub mod affine;
pub mod linear;
pub mod two_piece;

use crate::align::cost_models::AlignmentCostModel;

/// Describes what action a backtrace step takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BacktraceOp {
    /// Consume one graph node and one query position (match/mismatch).
    /// Switch to the given state after moving to the predecessor.
    Diagonal { next_state: u8 },
    /// Consume one graph node, stay at same query position (deletion).
    /// Switch to the given state after moving to the predecessor.
    Delete { next_state: u8 },
    /// Stay at same graph node, consume one query position (insertion).
    /// Transition within the same node to the given state.
    Insert { next_state: u8 },
    /// Stay at the same cell (same graph node, same query position); only
    /// switch to the given state. Used when the current state's value came
    /// from another state at the same cell — e.g. `M[v][q]` came from
    /// `D[v][q]` (closing a deletion) or `I[v][q]` (closing an insertion).
    SwitchState { next_state: u8 },
    /// Alignment is complete (reached origin).
    Done,
}

/// A DP kernel encapsulates how many states to track per (node, query-position)
/// cell and how to compute transitions for a specific gap cost model.
///
/// All data slices follow the layout: `data[state * w + qi]`
/// where `w` is the band width and `qi = q - qlo`.
pub trait DPKernel: Sized + 'static {
    /// Number of DP state arrays allocated per band.
    /// 1 for linear, 3 for affine, 5 for two-piece affine.
    const STATES: usize;

    /// Human-readable names for each state slot, in data-layout order.
    /// Length must equal `STATES`. Used to write debug cell TSV headers.
    const STATE_NAMES: &'static [&'static str];

    /// The cost model type this kernel is paired with.
    type Costs: AlignmentCostModel;

    // ── Initialization ──────────────────────────────────────────────────────

    /// Populate the start-sentinel band data slice (length `STATES * w`).
    /// `qlo` is the first query position in this band (usually 0).
    /// Called once before the forward pass.
    fn init_start(data: &mut [u32], w: usize, qlo: usize, costs: &Self::Costs);

    // ── Forward pass ────────────────────────────────────────────────────────

    /// Fold ONE predecessor band into both the current band's gap states and
    /// the shared diagonal-minimum buffer in a single pass over `src`.
    ///
    /// Called once per predecessor band before [`Self::finalize`]. Subsequent
    /// calls shrink `dst` (gap states) and `diag_buf` (diagonal minimum) via
    /// `min`, so the ordering of predecessor calls is irrelevant.
    ///
    /// - `dst`: current band data, len = `STATES * w_dst`. The gap states
    ///   (e.g. `D` for affine, `D1`/`D2` for two-piece) are updated in place.
    ///   M and I-like states are NOT touched here — they are written by
    ///   `finalize`.
    /// - `diag_buf`: shared buffer of length `w_dst` accumulating
    ///   `min(all-states at q-1 in predecessor)`. Substitution cost is NOT
    ///   applied here — it is applied once in [`Self::finalize`] after every
    ///   predecessor has contributed, so a single symbol comparison suffices
    ///   regardless of predecessor count.
    /// - `src`: predecessor band data, len = `STATES * w_src`.
    fn accumulate_pred(
        dst: &mut [u32],
        diag_buf: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Self::Costs,
    );

    /// Like [`Self::accumulate_pred`] but writes destinations with `=` instead
    /// of `min`-folding, so the caller can skip the per-band
    /// `diag_buf.fill(INF)` for **single-predecessor bands**.
    ///
    /// Preconditions enforced by the caller:
    ///
    /// - `dst[D_state][..]` is pre-initialized to [`crate::align::engine::dp::INF`]
    ///   (this is true after `BandMatrix::reset_for_bands`).
    /// - Every cell in `diag_buf[..w_dst]` is covered by Section B + C of this
    ///   call, i.e. the destination band's range equals the predecessor's
    ///   shifted range. This holds for single-pred bands because the
    ///   destination is constructed as `[pred.qlo+1, pred.qhi+1)`.
    ///
    /// For multi-predecessor bands the caller must keep the `fill(INF)` +
    /// `accumulate_pred` path: subsequent `min`-folds against an uninitialized
    /// `diag_buf` are unsound.
    fn init_pred(
        dst: &mut [u32],
        diag_buf: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Self::Costs,
    );

    /// Apply substitution cost to `diag_buf`, then compute sequential states
    /// (M and I-like) left-to-right. The gap states in `data` must already be
    /// filled (via `accumulate_pred` calls) before calling this.
    fn finalize(
        data: &mut [u32],
        w: usize,
        qlo: usize,
        symbol: u8,
        query: &[u8],
        diag_buf: &mut [u32],
        costs: &Self::Costs,
    );

    // ── DPGet compatibility ──────────────────────────────────────────────────

    /// M (match/mismatch) score at query index `qi`.
    fn m_at(data: &[u32], w: usize, qi: usize) -> u32;

    /// D (deletion from graph) score at query index `qi`.
    /// Returns the minimum across all deletion states for multi-state models.
    /// For models without an explicit deletion state (e.g. `LinearKernel`), return `INF`.
    fn d_at(data: &[u32], w: usize, qi: usize) -> u32;

    /// I (insertion into query) score at query index `qi`.
    /// Returns the minimum across all insertion states for multi-state models.
    /// For models without an explicit insertion state (e.g. `LinearKernel`), return `INF`.
    fn i_at(data: &[u32], w: usize, qi: usize) -> u32;

    // ── Backtrace ────────────────────────────────────────────────────────────
    //
    // The backtrace methods take a single-column view `col` of length `STATES`,
    // indexed directly by state id (`col[0]` = M, `col[1]` = D or D1, ...).
    // This avoids materialising a full `n_states * mp1` per-node array on every
    // backtrace step — the forward pass keeps its dense band layout, but
    // backtrace only ever needs the column at `q` (and occasionally `q-1`).

    /// Return `(best_score, best_state_id)` from a single column.
    /// State IDs correspond to `crate::align::AlignState as u8` by convention.
    fn best_state(col: &[u32]) -> (u32, u8);

    /// Cost of transitioning into a deletion (D) state from a predecessor column
    /// at query index `q`, and the state in the predecessor the deletion
    /// originated from.
    ///
    /// Returns `(transition_cost, pred_state_id)`.
    fn del_transition_cost(col: &[u32], costs: &Self::Costs) -> (u32, u8);

    /// Per-cell backtrace decision.
    ///
    /// `col` is the column at the current query position `q`; `col_prev` is the
    /// column at `q-1` or `None` when `q == 0`.
    fn backtrace_op(
        col: &[u32],
        col_prev: Option<&[u32]>,
        state: u8,
        costs: &Self::Costs,
    ) -> BacktraceOp;
}
