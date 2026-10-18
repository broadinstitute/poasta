//! The seed matches as points in the (T₁, T₂) plane, and the savings of the
//! best chain that starts at each. `docs/chain-heuristic-proofs.md`,
//! section 6, proves the queries exact.
//!
//! Each match m has savings σ(m): its own charge c plus the best
//! continuation. It becomes c / g sub-points with credits σ, σ − g, …,
//! σ − c + g at T(m) − (c − s·g)·(1, 1), where g is the grid step of all
//! T coordinates. The savings at a point q are the largest credit among
//! sub-points at or above q in both coordinates. Credit levels are nested,
//! so a search over levels finds it, galloping from the level of the
//! previous query.

use std::collections::BinaryHeap;

/// A seed match, by its index in `SeedMatches`.
pub type MatchId = u32;

const NONE: MatchId = MatchId::MAX;

/// A point in the (T₁, T₂) plane.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Point {
    pub t1: i64,
    pub t2: i64,
}

impl Point {
    /// This point moved by d along the diagonal.
    #[inline]
    pub fn shifted(self, d: i64) -> Point {
        Point {
            t1: self.t1 + d,
            t2: self.t2 + d,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SubPoint {
    pub p: Point,
    pub m: MatchId,
}

/// Sub-points with one credit, sorted by T₁ descending, and for every
/// prefix the largest T₂ and where it is. One binary search then answers a
/// dominance test and names a dominating match.
#[derive(Default, Clone)]
struct Level {
    pts: Vec<SubPoint>,
    /// `best[ix]`: the largest T₂ in `pts[..=ix]`, and its index.
    best: Vec<(i64, u32)>,
}

impl Level {
    fn fix_best(&mut self, from: usize) {
        let mut acc = if from == 0 {
            (i64::MIN, 0)
        } else {
            self.best[from - 1]
        };
        for ix in from..self.pts.len() {
            if self.pts[ix].p.t2 > acc.0 {
                acc = (self.pts[ix].p.t2, ix as u32);
            }
            self.best[ix] = acc;
        }
    }

    fn insert(&mut self, sp: SubPoint) {
        let at = self.pts.partition_point(|q| q.p.t1 > sp.p.t1);
        self.pts.insert(at, sp);
        self.best.insert(at, (0, 0));
        self.fix_best(at);
    }

    fn remove(&mut self, t1: i64, m: MatchId) {
        let from = self.pts.partition_point(|q| q.p.t1 > t1);
        let at = from
            + self.pts[from..]
                .iter()
                .position(|q| q.m == m)
                .expect("sub-point present");
        self.pts.remove(at);
        self.best.remove(at);
        self.fix_best(at);
    }

    /// A match with a sub-point at or above q in both coordinates, if any.
    #[inline]
    fn dominator(&self, q: Point) -> Option<MatchId> {
        let n = self.pts.partition_point(|sp| sp.p.t1 >= q.t1);
        (n > 0 && self.best[n - 1].0 >= q.t2).then(|| self.pts[self.best[n - 1].1 as usize].m)
    }
}

/// The savings at a point, and the match whose sub-point gives them.
#[cfg(any(test, feature = "viz"))]
#[derive(Clone, Copy, Debug)]
pub struct Credit {
    pub credit: i64,
    pub owner: Option<MatchId>,
}

/// A credit level with its sub-points, T₁ descending.
#[cfg(any(test, feature = "viz"))]
#[derive(Clone, Debug)]
pub struct Staircase {
    pub credit: i64,
    pub points: Vec<SubPoint>,
}

#[derive(Clone, Debug)]
struct MatchState {
    t: Point,
    /// Savings of the best chain that starts with this match. 0 once pruned.
    sigma: i64,
    /// The match whose sub-point gave the savings, or `NONE`.
    parent: MatchId,
    /// Whether the savings may be too high after a prune.
    dirty: bool,
}

#[derive(Default, Clone)]
pub(crate) struct Contours {
    charge: i64,
    grid: i64,
    matches: Vec<MatchState>,
    /// `levels[v]` holds the sub-points with credit v × grid. Level 0 is
    /// implicit: credit 0 is always available.
    levels: Vec<Level>,
    /// Reverse of `parent`. A list may hold stale entries: check `parent`.
    children: Vec<Vec<MatchId>>,
    /// Dirty matches by T₁, largest first. Entries for matches that are no
    /// longer dirty are skipped when popped.
    dirty_heap: BinaryHeap<(i64, MatchId)>,
    /// The level of the last query, where the next one starts searching.
    hint: usize,
}

impl Contours {
    pub fn new(points: Vec<Point>, charge: i64, grid: i64) -> Self {
        let n = points.len();
        let mut s = Self {
            charge,
            grid,
            matches: points
                .into_iter()
                .map(|t| MatchState {
                    t,
                    sigma: 0,
                    parent: NONE,
                    dirty: false,
                })
                .collect(),
            levels: vec![Level::default()],
            children: vec![Vec::new(); n],
            ..Default::default()
        };
        // A sub-point of m lies strictly below T(m), so every sub-point that
        // can dominate T(m) belongs to a match with a larger T₁. Processing
        // by T₁ descending sees all of them.
        let mut order: Vec<MatchId> = (0..n as MatchId).collect();
        order.sort_unstable_by_key(|&m| std::cmp::Reverse(s.matches[m as usize].t.t1));
        for m in order {
            s.compute(m, 0);
        }
        s
    }

    /// Credit levels a match spans: its savings minus 0, g, …, c − g.
    fn levels_of(&self, sigma: i64) -> impl Iterator<Item = (usize, i64)> + use<> {
        let (charge, grid) = (self.charge, self.grid);
        (0..charge / grid).map(move |s| (((sigma - s * grid) / grid) as usize, -charge + s * grid))
    }

    /// Set the savings of m from the current sub-points.
    fn compute(&mut self, m: MatchId, hint: usize) {
        let t = self.matches[m as usize].t;
        let level = self.search(t, hint);
        let owner = self.owner(level, t);
        let sigma = self.charge + level as i64 * self.grid;
        let ms = &mut self.matches[m as usize];
        ms.sigma = sigma;
        ms.parent = owner.unwrap_or(NONE);
        if let Some(o) = owner {
            self.children[o as usize].push(m);
        }
        for (v, d) in self.levels_of(sigma) {
            if self.levels.len() <= v {
                self.levels.resize(v + 1, Level::default());
            }
            self.levels[v].insert(SubPoint { p: t.shifted(d), m });
        }
    }

    fn clear(&mut self, m: MatchId) {
        let ms = &self.matches[m as usize];
        let t1 = ms.t.t1;
        for (v, d) in self.levels_of(ms.sigma) {
            self.levels[v].remove(t1 + d, m);
        }
        self.matches[m as usize].sigma = 0;
    }

    /// The highest level with a sub-point that dominates q, searching from
    /// `hint`. Exact only if no dirty sub-point can dominate q.
    fn search(&self, q: Point, hint: usize) -> usize {
        last_true(hint, self.levels.len() - 1, |v| {
            v == 0 || self.levels[v].dominator(q).is_some()
        })
    }

    fn owner(&self, level: usize, q: Point) -> Option<MatchId> {
        if level == 0 {
            None
        } else {
            self.levels[level].dominator(q)
        }
    }

    /// Recompute the dirty matches that can dominate a point with T₁ ≥ t1. A
    /// dirty sub-point lies at or below T(m) − g, so only dirty matches with
    /// T₁ − g ≥ t1 qualify. They go largest T₁ first, which is the build
    /// order: their savings depend only on matches with a larger T₁, which
    /// are repaired before them.
    fn repair_from(&mut self, t1: i64) {
        while let Some(&(mt1, m)) = self.dirty_heap.peek() {
            if mt1 - self.grid < t1 {
                break;
            }
            self.dirty_heap.pop();
            let ms = &self.matches[m as usize];
            if !ms.dirty || ms.sigma == 0 {
                continue;
            }
            let hint = ((ms.sigma - self.charge).max(0) / self.grid) as usize;
            self.clear(m);
            self.matches[m as usize].dirty = false;
            self.compute(m, hint);
        }
    }

    /// Exact savings at q.
    #[inline]
    pub fn credit(&mut self, q: Point) -> i64 {
        self.repair_from(q.t1);
        self.hint = self.search(q, self.hint);
        self.hint as i64 * self.grid
    }

    /// Exact savings at q, and the match that gives them.
    #[cfg(any(test, feature = "viz"))]
    pub fn best(&mut self, q: Point) -> Credit {
        let credit = self.credit(q);
        Credit {
            credit,
            owner: self.owner(self.hint, q),
        }
    }

    /// Recompute every dirty match.
    #[cfg(any(test, feature = "viz"))]
    pub fn repair(&mut self) {
        self.repair_from(i64::MIN / 2);
    }

    /// Match m followed by the matches its savings came through.
    #[cfg(any(test, feature = "viz"))]
    pub fn chain(&self, mut m: MatchId) -> Vec<MatchId> {
        let mut out = Vec::new();
        while m != NONE {
            out.push(m);
            m = self.matches[m as usize].parent;
        }
        out
    }

    #[cfg(any(test, feature = "viz"))]
    pub fn t(&self, m: MatchId) -> Point {
        self.matches[m as usize].t
    }

    /// Savings of match m. Exact after `repair`.
    #[cfg(any(test, feature = "viz"))]
    pub fn sigma(&self, m: MatchId) -> i64 {
        self.matches[m as usize].sigma
    }

    /// Sub-points by credit level. Level 0 is implicit and left out.
    #[cfg(any(test, feature = "viz"))]
    pub fn staircases(&self) -> Vec<Staircase> {
        self.levels
            .iter()
            .enumerate()
            .skip(1)
            .map(|(v, l)| Staircase {
                credit: v as i64 * self.grid,
                points: l.pts.clone(),
            })
            .collect()
    }

    /// Whether match m is still in the contours.
    pub fn is_live(&self, m: MatchId) -> bool {
        self.matches[m as usize].sigma > 0
    }

    /// Remove match m. Every match whose savings came through m, directly
    /// or not, becomes dirty. Its sub-points stay where they were until a
    /// query needs them. They only lie below and left of their matches, so
    /// a query outside their bounding box cannot see them, and every other
    /// sub-point is exact. Returns the number of matches marked.
    pub fn prune(&mut self, m: MatchId) -> usize {
        if !self.is_live(m) {
            return 0;
        }
        self.clear(m);
        self.matches[m as usize].dirty = false;
        let mut marked = 0;
        let mut stack = vec![m];
        while let Some(p) = stack.pop() {
            let kids = std::mem::take(&mut self.children[p as usize]);
            for ch in kids {
                let cs = &mut self.matches[ch as usize];
                if cs.parent != p || cs.dirty || cs.sigma == 0 {
                    continue;
                }
                cs.dirty = true;
                self.dirty_heap.push((cs.t.t1, ch));
                marked += 1;
                stack.push(ch);
            }
        }
        while self.levels.len() > 1 && self.levels.last().is_some_and(|l| l.pts.is_empty()) {
            self.levels.pop();
        }
        marked
    }
}

/// The largest v ≤ top with ok(v), where ok(0) holds and ok is true up to
/// the answer and false after it. Gallops from `hint` to bracket the
/// answer, then bisects.
fn last_true(hint: usize, top: usize, ok: impl Fn(usize) -> bool) -> usize {
    let hint = hint.min(top);
    let (mut lo, mut hi);
    if ok(hint) {
        lo = hint;
        let mut step = 1;
        loop {
            let probe = lo + step;
            if probe > top {
                hi = top;
                break;
            }
            if ok(probe) {
                lo = probe;
                step *= 2;
            } else {
                hi = probe - 1;
                break;
            }
        }
    } else {
        hi = hint - 1;
        let mut step = 1;
        loop {
            let probe = hi.saturating_sub(step);
            if ok(probe) {
                lo = probe;
                break;
            }
            hi = probe - 1;
            step *= 2;
        }
    }
    while lo < hi {
        let mid = (lo + hi).div_ceil(2);
        if ok(mid) {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    lo
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force savings: σ(m) = c + max over matches m' of
    /// min(σ(m'), U(m') − T(m)), with U(m') = T(m') − c + σ(m').
    fn brute_sigma(t: &[(i64, i64)], alive: &[bool], c: i64) -> Vec<i64> {
        // Fixed point of σ(m) = c + max(0, max over m' ≠ m of savings),
        // over every live match, with no ordering assumed.
        let n = t.len();
        let mut s: Vec<i64> = (0..n).map(|m| if alive[m] { c } else { 0 }).collect();
        loop {
            let mut changed = false;
            for m in (0..n).filter(|&m| alive[m]) {
                let best = (0..n)
                    .filter(|&x| alive[x] && x != m)
                    .map(|x| savings(t[x], s[x], c, t[m]))
                    .max()
                    .unwrap_or(0)
                    .max(0);
                if c + best != s[m] {
                    s[m] = c + best;
                    changed = true;
                }
            }
            if !changed {
                return s;
            }
        }
    }

    fn savings(tm: (i64, i64), sig: i64, c: i64, q: (i64, i64)) -> i64 {
        sig.min(tm.0 - c + sig - q.0).min(tm.1 - c + sig - q.1)
    }

    #[test]
    fn smooth_contours_match_brute_force_with_pruning() {
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut rnd = |n: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x % n) as i64
        };
        for case in 0..2000 {
            let (c, g) = [(4, 2), (4, 1), (2, 2), (6, 3)][case % 4];
            let n = 1 + rnd(30) as usize;
            let t: Vec<(i64, i64)> = (0..n).map(|_| (rnd(24) * g, rnd(24) * g)).collect();
            let mut cs = Contours::new(t.iter().map(|&(t1, t2)| Point { t1, t2 }).collect(), c, g);
            let mut alive = vec![true; n];
            for _ in 0..=n / 2 {
                let want = brute_sigma(&t, &alive, c);
                for _ in 0..10 {
                    let q = ((rnd(26) - 1) * g, (rnd(26) - 1) * g);
                    let expect = (0..n)
                        .filter(|&m| alive[m])
                        .map(|m| savings(t[m], want[m], c, q))
                        .max()
                        .unwrap_or(0)
                        .max(0);
                    cs.hint = rnd(8) as usize;
                    let got = cs.best(Point { t1: q.0, t2: q.1 });
                    assert_eq!(got.credit, expect, "case {case}");
                    assert_eq!(got.owner.is_some(), expect > 0, "case {case}");
                }
                cs.repair();
                assert_eq!(
                    (0..n as MatchId).map(|m| cs.sigma(m)).collect::<Vec<_>>(),
                    want,
                    "case {case}"
                );
                let m = rnd(n as u64) as usize;
                alive[m] = false;
                cs.prune(m as u32);
            }
        }
    }

    #[test]
    fn last_true_from_every_hint() {
        for top in 0..20 {
            for answer in 0..=top {
                for hint in 0..=top + 3 {
                    assert_eq!(
                        last_true(hint, top, |v| v <= answer),
                        answer,
                        "top {top} answer {answer} hint {hint}"
                    );
                }
            }
        }
    }
}
