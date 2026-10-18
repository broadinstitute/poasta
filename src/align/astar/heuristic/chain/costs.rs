//! Seed arithmetic for one cost model: what an unmatched seed costs at
//! least, what the heuristic charges per seed, and the per-state
//! adjustments that keep it admissible and consistent. The proofs are in
//! `docs/chain-heuristic-proofs.md`, section 2.

use crate::align::AlignState;
use crate::align::cost_models::AlignmentCostModel;

/// Seed costs for one cost model, seed length and seed count.
#[derive(Clone, Debug)]
pub(crate) struct SeedCosts {
    /// Least cost any alignment pays per seed it does not match.
    pub miss: usize,
    /// What `h` charges per seed. 0 turns the heuristic off.
    pub charge: i64,
    /// Smallest gap extend cost: the price of one position of shift.
    pub e_min: i64,
    /// Charge for a half-aligned seed without a rest path, by `AlignState`:
    /// the charge in `Match`, at most one extend in a deletion state, and 0
    /// in an insertion state.
    pub partial: [i64; 5],
    /// Insertion-state discount by position inside a seed, by `AlignState`.
    /// Empty for the other states.
    pub ins_discount: [Vec<i64>; 5],
    /// Grid step of every potential and T coordinate.
    pub lattice: i64,
}

impl SeedCosts {
    /// `states` is the alignment graph's state count: 1 (linear), 3
    /// (affine) or 5 (two-piece).
    pub fn new<C: AlignmentCostModel>(
        costs: &C,
        states: usize,
        k: usize,
        r: usize,
        n_seeds: usize,
    ) -> Self {
        let miss = miss(costs, k, r, n_seeds);
        let e_min = gap_pieces(costs).iter().map(|&(_, e)| e).min().unwrap_or(0) as i64;

        let ins_pieces = insertion_pieces(costs, states);

        // The largest charge that keeps every edge consistent: at most
        // `miss` and one mismatch, and for each insertion piece at most k·e,
        // so the discount exists, and at most e + ⌊o/2⌋, so opening an
        // insertion costs at least the drop in `h`: o + e ≥ c + D(0).
        let mut charge = if e_min == 0 {
            0
        } else {
            (miss as i64).min(costs.mismatch() as i64)
        };
        for &(_, open, extend) in &ins_pieces {
            charge = charge.min(k as i64 * extend).min(extend + open / 2);
        }
        let mut ins_discount: [Vec<i64>; 5] = Default::default();
        if charge > 0 {
            for &(state, _, extend) in &ins_pieces {
                ins_discount[state] = insertion_discount(charge, extend, k);
            }
        }

        let mut partial = [charge, 0, 0, 0, 0];
        if states >= 3 {
            partial[AlignState::Deletion as usize] = charge.min(costs.gap_extend() as i64);
        }
        if states == 5 {
            partial[AlignState::Deletion2 as usize] = charge.min(costs.gap_extend2() as i64);
        }

        let mut lattice = gcd(charge, e_min);
        for &p in &partial {
            lattice = gcd(lattice, p);
        }
        for d in &ins_discount {
            for &x in d {
                lattice = gcd(lattice, x);
            }
        }

        Self {
            miss,
            charge,
            e_min,
            partial,
            ins_discount,
            lattice: lattice.max(1),
        }
    }
}

/// (state, open, extend) of each insertion state the alignment graph has.
fn insertion_pieces<C: AlignmentCostModel>(costs: &C, states: usize) -> Vec<(usize, i64, i64)> {
    let mut pieces = Vec::new();
    if states >= 3 {
        pieces.push((
            AlignState::Insertion as usize,
            costs.gap_open() as i64,
            costs.gap_extend() as i64,
        ));
    }
    if states == 5 {
        pieces.push((
            AlignState::Insertion2 as usize,
            costs.gap_open2() as i64,
            costs.gap_extend2() as i64,
        ));
    }
    pieces
}

/// (open, extend) of every gap piece with a nonzero extend cost.
pub(crate) fn gap_pieces<C: AlignmentCostModel>(costs: &C) -> Vec<(usize, usize)> {
    [
        (costs.gap_open() as usize, costs.gap_extend() as usize),
        (costs.gap_open2() as usize, costs.gap_extend2() as usize),
    ]
    .into_iter()
    .filter(|&(_, e)| e > 0)
    .collect()
}

/// Least cost per seed of any edit event that makes a seed unmatched. A
/// seed without a match holds r + 1 mismatches, a gap inside it, or part
/// of an insertion run. A run over s seeds has at least (s − 2)k + 2
/// bases, and its cost is shared by the s seeds it touches.
pub(crate) fn miss<C: AlignmentCostModel>(costs: &C, k: usize, r: usize, n_seeds: usize) -> usize {
    let pieces = gap_pieces(costs);
    let g = |len: usize| {
        pieces
            .iter()
            .map(|&(o, e)| o + e * len)
            .min()
            .unwrap_or(usize::MAX)
    };
    let mut miss = ((r + 1) * costs.mismatch() as usize).min(g(1));
    for s in 2..=n_seeds {
        miss = miss.min(g((s - 2) * k + 2) / s);
    }
    miss
}

/// The insertion-state discount D(t) for one gap piece with extend `e`, at
/// each position t inside a seed of length k, for charge c ≤ k·e. It is the
/// smallest non-negative D with
///
/// - D(0) ≥ c − e and D(t) ≥ c − e(k − t + 1), what an insertion run that is
///   already open can save against the charges of the seeds it enters, so
///   `h` stays admissible,
/// - D(0) − D(1) ≥ c − e and D(t + 1) − D(t) ≤ e, so an insertion step costs
///   at least the drop in `h`.
///
/// The lower bounds of the first row already satisfy the second when
/// c ≤ k·e, so they are the table: c − e·(d + 1), clamped at 0, with d the
/// distance from t to the next seed start.
pub(crate) fn insertion_discount(c: i64, e: i64, k: usize) -> Vec<i64> {
    (0..k)
        .map(|t| (c - e * ((k - t) % k + 1) as i64).max(0))
        .collect()
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cost_models::affine::Affine;
    use crate::align::cost_models::linear::Linear;
    use crate::align::cost_models::two_piece::TwoPieceAffine;

    #[test]
    fn miss_values() {
        assert_eq!(miss(&Affine::new(4, 6, 2), 12, 0, 80), 4);
        assert_eq!(miss(&Affine::new(4, 6, 2), 12, 1, 80), 5);
        assert_eq!(miss(&Linear::new(4, 2), 12, 0, 80), 2);
        assert_eq!(miss(&TwoPieceAffine::new(4, 6, 2, 24, 1), 12, 0, 80), 4);
        assert_eq!(miss(&TwoPieceAffine::new(4, 6, 2, 24, 1), 12, 1, 80), 5);
        assert_eq!(miss(&Affine::new(0, 6, 2), 12, 0, 80), 0);
    }

    #[test]
    fn charges_and_discounts() {
        let c = SeedCosts::new(&Affine::new(4, 6, 2), 3, 8, 0, 100);
        assert_eq!((c.charge, c.e_min, c.lattice), (4, 2, 2));
        assert_eq!(
            c.ins_discount[AlignState::Insertion as usize],
            vec![2, 0, 0, 0, 0, 0, 0, 0]
        );
        assert_eq!(c.partial[AlignState::Deletion as usize], 2);

        // The full charge also for two-piece: no cap.
        let c = SeedCosts::new(&TwoPieceAffine::new(4, 6, 2, 24, 1), 5, 8, 0, 100);
        assert_eq!((c.charge, c.e_min, c.lattice), (4, 1, 1));

        // Mismatch 1 below the piece-one extend 2: D(0) is clamped at 0.
        let c = SeedCosts::new(&TwoPieceAffine::new(1, 2, 2, 6, 1), 5, 4, 0, 100);
        assert_eq!(c.charge, 1);
        assert!(
            c.ins_discount[AlignState::Insertion as usize]
                .iter()
                .all(|&d| d == 0)
        );

        assert_eq!(
            SeedCosts::new(&Affine::new(0, 6, 2), 3, 8, 0, 100).charge,
            0
        );
    }

    /// The charge and discounts by their definition: lower the charge until
    /// every insertion piece has a discount table, found by relaxing the
    /// constraints, with c + max D ≤ o + e.
    fn by_definition<C: AlignmentCostModel>(
        costs: &C,
        states: usize,
        k: usize,
        r: usize,
        n: usize,
    ) -> (i64, [Vec<i64>; 5]) {
        let relax = |c: i64, e: i64| -> Option<Vec<i64>> {
            if c > e * k as i64 {
                return None;
            }
            let mut d: Vec<i64> = (0..k)
                .map(|t| {
                    if t == 0 {
                        (c - e).max(0)
                    } else {
                        (c - e * (k - t + 1) as i64).max(0)
                    }
                })
                .collect();
            loop {
                let mut changed = false;
                if k > 1 && d[0] < d[1] + c - e {
                    d[0] = d[1] + c - e;
                    changed = true;
                }
                for t in (1..k).rev() {
                    if d[t] < d[(t + 1) % k] - e {
                        d[t] = d[(t + 1) % k] - e;
                        changed = true;
                    }
                }
                if !changed {
                    return Some(d);
                }
            }
        };
        let e_min = gap_pieces(costs).iter().map(|&(_, e)| e).min().unwrap_or(0);
        let mut c = (miss(costs, k, r, n) as i64).min(costs.mismatch() as i64);
        while c > 0 {
            let mut tables: [Vec<i64>; 5] = Default::default();
            let fits =
                insertion_pieces(costs, states)
                    .into_iter()
                    .all(|(state, o, e)| match relax(c, e) {
                        Some(d) if c + d.iter().copied().max().unwrap_or(0) <= o + e => {
                            tables[state] = d;
                            true
                        }
                        _ => false,
                    });
            if fits && e_min > 0 {
                return (c, tables);
            }
            c -= 1;
        }
        (0, Default::default())
    }

    /// A charge and its five insertion-discount tables.
    type Tables = (i64, [Vec<i64>; 5]);

    #[test]
    fn closed_forms_match_the_definition() {
        for k in 1..=9 {
            for mismatch in 0..=6 {
                for open in 0..=9 {
                    for extend in 0..=4 {
                        let check = |name: &str, (got, want): (Tables, Tables)| {
                            assert_eq!(got, want, "{name} k={k} x={mismatch} o={open} e={extend}");
                        };
                        let linear = Linear::new(mismatch, extend);
                        let affine = Affine::new(mismatch, open, extend);
                        let two_piece = TwoPieceAffine::new(
                            mismatch,
                            open,
                            extend,
                            open * 3,
                            extend.saturating_sub(1).max(1),
                        );
                        for n in [1, 2, 50] {
                            let new = |c: SeedCosts| (c.charge, c.ins_discount);
                            check(
                                "linear",
                                (
                                    new(SeedCosts::new(&linear, 1, k, 0, n)),
                                    by_definition(&linear, 1, k, 0, n),
                                ),
                            );
                            check(
                                "affine",
                                (
                                    new(SeedCosts::new(&affine, 3, k, 0, n)),
                                    by_definition(&affine, 3, k, 0, n),
                                ),
                            );
                            check(
                                "affine r=1",
                                (
                                    new(SeedCosts::new(&affine, 3, k, 1, n)),
                                    by_definition(&affine, 3, k, 1, n),
                                ),
                            );
                            check(
                                "two-piece",
                                (
                                    new(SeedCosts::new(&two_piece, 5, k, 0, n)),
                                    by_definition(&two_piece, 5, k, 0, n),
                                ),
                            );
                        }
                    }
                }
            }
        }
    }

    /// Every insertion run, at every start inside a stretch of seeds, costs
    /// at least `miss` per seed it touches. That is what `miss` promises for
    /// the runs a completion opens itself.
    #[test]
    fn insertion_runs_cost_at_least_miss_per_seed() {
        for (k, costs) in [
            (4, Affine::new(4, 6, 2)),
            (8, Affine::new(4, 6, 2)),
            (3, Affine::new(1, 2, 1)),
        ] {
            let n = 6;
            let m = miss(&costs, k, 0, n);
            for start in 0..n * k {
                for len in 1..=(n * k - start) {
                    let touched = (start + len - 1) / k - start / k + 1;
                    let cost = costs.gap_open() as usize + len * costs.gap_extend() as usize;
                    assert!(cost >= m * touched, "k {k} start {start} len {len}");
                }
            }
        }
    }
}
