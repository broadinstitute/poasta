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

    /// The cost model type this kernel is paired with.
    type Costs: AlignmentCostModel;

    // ── Initialization ──────────────────────────────────────────────────────

    /// Populate the start-sentinel band data slice (length `STATES * w`).
    /// `qlo` is the first query position in this band (usually 0).
    /// Called once before the forward pass.
    fn init_start(data: &mut [u32], w: usize, qlo: usize, costs: &Self::Costs);

    // ── Forward pass ────────────────────────────────────────────────────────

    /// Accumulate gap-state (D-like) contributions from ONE predecessor band
    /// into the current band's data. Called once per predecessor band.
    ///
    /// `dst`: current band data, len = `STATES * w_dst`
    /// `src`: predecessor band data, len = `STATES * w_src`
    fn accumulate_gap(
        dst: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
        costs: &Self::Costs,
    );

    /// Accumulate the diagonal predecessor minimum into `diag_buf` (length `w_dst`).
    /// `diag_buf[qi]` accumulates `min(all-states at q-1 in predecessor)`.
    /// Called once per predecessor band before `finalize`.
    /// No `costs` parameter — substitution cost is applied later in [`Self::finalize`].
    fn accumulate_diag(
        diag_buf: &mut [u32],
        w_dst: usize,
        dst_qlo: usize,
        src: &[u32],
        w_src: usize,
        src_qlo: usize,
    );

    /// Apply substitution cost to `diag_buf`, then compute sequential states
    /// (M and I-like) left-to-right. The gap states in `data` must already be
    /// filled (via `accumulate_gap` calls) before calling this.
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
