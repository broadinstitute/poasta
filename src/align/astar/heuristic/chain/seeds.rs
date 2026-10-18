//! The query cut into seeds, and the offset arithmetic around them.

/// Seeds `query[j·k..(j + 1)·k]` for j < n, matched with at most r
/// substitutions. Seeds in a free query tail are not among them.
#[derive(Default, Clone, Debug)]
pub(crate) struct Seeds {
    query: Vec<u8>,
    pub k: usize,
    pub r: usize,
    pub n: usize,
}

impl Seeds {
    pub fn reset(&mut self, query: &[u8], k: usize, r: usize, n: usize) {
        self.query.clear();
        self.query.extend_from_slice(query);
        self.k = k;
        self.r = r;
        self.n = n;
    }

    pub fn seed(&self, j: usize) -> &[u8] {
        &self.query[j * self.k..(j + 1) * self.k]
    }

    /// Seeds that start at or after offset i.
    #[inline]
    pub fn ahead(&self, i: usize) -> usize {
        self.n - i.div_ceil(self.k).min(self.n)
    }

    /// Whether offset i lies strictly inside a seed.
    #[inline]
    pub fn inside(&self, i: usize) -> bool {
        !i.is_multiple_of(self.k) && i / self.k < self.n
    }

    /// The rest of the seed that offset i is inside, from i to its end.
    #[inline]
    pub fn rest(&self, i: usize) -> &[u8] {
        &self.query[i..(i / self.k + 1) * self.k]
    }

    /// Where the match points of seed j sit: just after its first base.
    #[inline]
    pub fn match_offset(&self, j: usize) -> usize {
        j * self.k + 1
    }

    /// The seed whose match points sit at offset i, if any.
    #[inline]
    pub fn seed_matched_at(&self, i: usize) -> Option<usize> {
        let j = i.checked_sub(1)? / self.k;
        (i == self.match_offset(j) && j < self.n).then_some(j)
    }
}

/// One step of a depth-first walk that spells part of a seed: the node
/// reached, the position in the seed, and the substitutions spent.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Walk<N> {
    pub node: N,
    pub pos: usize,
    pub mismatches: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn offsets() {
        let mut s = Seeds::default();
        // Seeds ACG, TAC, and a tail of 2 that is not a seed.
        s.reset(b"ACGTACGT", 3, 0, 2);
        assert_eq!(s.seed(1), b"TAC");
        assert_eq!(
            (0..=8).map(|i| s.ahead(i)).collect::<Vec<_>>(),
            [2, 1, 1, 1, 0, 0, 0, 0, 0]
        );
        assert_eq!(
            (0..=8).filter(|&i| s.inside(i)).collect::<Vec<_>>(),
            [1, 2, 4, 5]
        );
        assert_eq!(s.rest(4), b"AC");
        assert_eq!(
            (0..=8)
                .filter_map(|i| s.seed_matched_at(i))
                .collect::<Vec<_>>(),
            [0, 1]
        );
        assert_eq!(s.seed_matched_at(4), Some(1));

        s.reset(b"ACGT", 1, 0, 4);
        assert_eq!(
            (0..=5).map(|i| s.seed_matched_at(i)).collect::<Vec<_>>(),
            [None, Some(0), Some(1), Some(2), Some(3), None]
        );
    }
}
