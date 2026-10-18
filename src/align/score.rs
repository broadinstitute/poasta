//! Score representation for the aligner: a bare `u32` with [`INF`] as the
//! "unreachable, not yet visited" sentinel instead of an `Option` or enum
//! wrapper. A visited cell defaults to `INF`, `update_score_if_lower` is a
//! plain `new < *slot`, and backtrace equalities are exact `u32` comparisons.

/// Sentinel value for "unreachable / infinite cost".
///
/// Chosen so that `INF + any_realistic_cost` stays below `u32::MAX`, which lets
/// `sat` use a plain wrapping add + min instead of the more expensive
/// `saturating_add` (which requires overflow detection on every u32 lane and
/// maps to multiple instructions instead of `vpaddd` + `vpminud`).
///
/// Real scores must stay below `INF`. A score is at most
/// `sequence_len × max_cost_per_step`, so with per-step costs up to about 25
/// that holds for queries up to about 10 Mbp. Raise this constant (and
/// check `u32` headroom) before aligning longer queries or using much
/// larger costs.
pub const INF: u32 = 1 << 28;

/// Saturating-style cost addition, capped at [`INF`].
///
/// Uses `wrapping_add` rather than `saturating_add` so that the compiler can
/// lower the vectorised form to one `vpaddd` + `vpminud` pair instead of
/// the multi-instruction overflow-detection sequence that `u32::saturating_add`
/// requires.  Correctness relies on the invariant that real scores never exceed
/// `u32::MAX - max_single_step_cost`. See [`INF`] for the bound.
#[inline(always)]
pub fn sat(a: u32, b: u32) -> u32 {
    a.wrapping_add(b).min(INF)
}
