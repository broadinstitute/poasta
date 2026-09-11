//! Shared score representation used by every alignment engine (canonical DP,
//! band-doubling, and A*): a bare `u32` with [`INF`] as the "unreachable /
//! not yet visited" sentinel, instead of an `Option`/enum wrapper.
//!
//! Promoted out of `engine::dp` (which re-exports these for compatibility)
//! so the A* engine can share the exact same arithmetic convention: a visited
//! cell defaults to `INF`, `update_score_if_lower` is a plain `new < *slot`,
//! and backtrace equalities are exact `u32` comparisons.

/// Sentinel value for "unreachable / infinite cost".
///
/// Chosen so that `INF + any_realistic_cost` stays below `u32::MAX`, which lets
/// `sat` use a plain wrapping add + min instead of the more expensive
/// `saturating_add` (which requires overflow detection on every u32 lane and
/// maps to multiple instructions instead of `vpaddd` + `vpminud`).
///
/// Upper bound on a real score: `sequence_len × max_cost_per_step`.  For
/// sequences up to ~10 million bp and costs up to ~1 000, that is ≈10^10 —
/// far below `1<<28 ≈ 2.7×10^8`.  If you align sequences of length > 100 M or
/// use very large cost values, raise this constant accordingly.
pub const INF: u32 = 1 << 28;

/// Saturating-style cost addition, capped at [`INF`].
///
/// Uses `wrapping_add` rather than `saturating_add` so that the compiler can
/// lower the vectorised form to a single `vpaddd` + `vpminud` pair instead of
/// the multi-instruction overflow-detection sequence that `u32::saturating_add`
/// requires.  Correctness relies on the invariant that real scores never exceed
/// `u32::MAX - max_single_step_cost`; see [`INF`] for the bound.
#[inline(always)]
pub fn sat(a: u32, b: u32) -> u32 {
    a.wrapping_add(b).min(INF)
}

#[inline(always)]
pub fn min3(a: u32, b: u32, c: u32) -> u32 {
    a.min(b).min(c)
}
