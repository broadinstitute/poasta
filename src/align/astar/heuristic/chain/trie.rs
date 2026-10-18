//! The query's seeds as a trie over query bytes, so the finder can test
//! every step with `AlignableGraph::is_symbol_equal` instead of comparing
//! bytes. An overriding `is_symbol_equal` (IUPAC codes) then still finds
//! every match, which admissibility depends on.

use super::seeds::Seeds;

pub(crate) type TrieNode = u32;

#[derive(Clone, Copy, Debug)]
pub(crate) struct TrieEdge {
    pub byte: u8,
    pub child: TrieNode,
}

/// Node 0 is the root. Every seed has length k, so the leaves are the nodes
/// without children, all at depth k. Identical seeds share their leaf.
#[derive(Default, Clone, Debug)]
pub(crate) struct SeedTrie {
    /// Children of node t: `edges[child_start[t]..child_start[t + 1]]`,
    /// sorted by byte.
    child_start: Vec<u32>,
    edges: Vec<TrieEdge>,
    /// Seed ids sorted by seed. Seeds ending at leaf t:
    /// `seeds[leaf_start[t]..leaf_start[t + 1]]`.
    seeds: Vec<u32>,
    leaf_start: Vec<u32>,
    /// Scratch, kept for its allocations: each new edge with its parent, and
    /// the path to the last seed.
    new_edges: Vec<(TrieNode, TrieEdge)>,
    path: Vec<TrieNode>,
}

impl SeedTrie {
    pub fn rebuild(&mut self, seeds: &Seeds) {
        // Sorted seeds share each prefix with the previous one, so one pass
        // creates each trie node once, children in byte order, and leaves in
        // the order of `self.seeds`. So when the seed at position `pos` of
        // that order creates a node, the leaves before it hold exactly `pos`
        // seeds, which is its `leaf_start`.
        self.seeds.clear();
        self.seeds.extend(0..seeds.n as u32);
        self.seeds.sort_unstable_by(|&a, &b| seeds.seed(a as usize).cmp(seeds.seed(b as usize)).then(a.cmp(&b)));
        self.new_edges.clear();
        self.leaf_start.clear();
        self.leaf_start.push(0);
        self.path.clear();
        self.path.push(0);
        let mut prev: Option<&[u8]> = None;
        for (pos, &j) in self.seeds.iter().enumerate() {
            let s = seeds.seed(j as usize);
            let common = prev.map_or(0, |p| p.iter().zip(s).take_while(|(a, b)| a == b).count());
            self.path.truncate(common + 1);
            for (d, &byte) in s.iter().enumerate().skip(common) {
                let child = self.leaf_start.len() as TrieNode;
                self.leaf_start.push(pos as u32);
                self.new_edges.push((self.path[d], TrieEdge { byte, child }));
                self.path.push(child);
            }
            prev = Some(s);
        }
        let n = self.leaf_start.len();
        self.leaf_start.push(seeds.n as u32);

        // Stable, so each parent's children stay in byte order.
        self.new_edges.sort_by_key(|&(parent, _)| parent);
        self.child_start.clear();
        self.child_start.resize(n + 1, 0);
        for &(parent, _) in &self.new_edges {
            self.child_start[parent as usize + 1] += 1;
        }
        for t in 0..n {
            self.child_start[t + 1] += self.child_start[t];
        }
        self.edges.clear();
        self.edges.extend(self.new_edges.iter().map(|&(_, e)| e));
    }

    pub fn children(&self, t: TrieNode) -> &[TrieEdge] {
        &self.edges[self.child_start[t as usize] as usize..self.child_start[t as usize + 1] as usize]
    }

    pub fn is_leaf(&self, t: TrieNode) -> bool {
        self.child_start[t as usize] == self.child_start[t as usize + 1]
    }

    pub fn seeds_at(&self, t: TrieNode) -> &[u32] {
        &self.seeds[self.leaf_start[t as usize] as usize..self.leaf_start[t as usize + 1] as usize]
    }

    #[cfg(test)]
    pub fn node_count(&self) -> usize {
        self.child_start.len() - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trie(q: &[u8], k: usize, n: usize) -> SeedTrie {
        let mut seeds = Seeds::default();
        seeds.reset(q, k, 0, n);
        let mut trie = SeedTrie::default();
        trie.rebuild(&seeds);
        trie
    }

    #[test]
    fn identical_seeds_share_a_leaf() {
        // Seeds ACG, TAC, ACG.
        let trie = trie(b"ACGTACACG", 3, 3);
        let leaves: Vec<u32> = (1..trie.node_count() as u32).filter(|&t| trie.is_leaf(t)).collect();
        assert_eq!(leaves.len(), 2);
        let mut all: Vec<Vec<u32>> = leaves.iter().map(|&t| trie.seeds_at(t).to_vec()).collect();
        all.sort();
        assert_eq!(all, vec![vec![0, 2], vec![1]]);
        let bytes: Vec<u8> = trie.children(0).iter().map(|e| e.byte).collect();
        assert_eq!(bytes, b"AT");
        // Inner nodes hold no seeds.
        for t in (0..trie.node_count() as u32).filter(|&t| !trie.is_leaf(t)) {
            assert!(trie.seeds_at(t).is_empty());
        }
    }

    #[test]
    fn walking_the_trie_spells_every_seed() {
        let q = b"ACGTTGCAACGTGGGA";
        let mut seeds = Seeds::default();
        let mut t = SeedTrie::default();
        for k in 1..=5 {
            let n = q.len() / k;
            seeds.reset(q, k, 0, n);
            t.rebuild(&seeds);
            for j in 0..n {
                let mut node = 0;
                for &b in &q[j * k..(j + 1) * k] {
                    node = t.children(node).iter().find(|e| e.byte == b).unwrap().child;
                }
                assert!(t.is_leaf(node) && t.seeds_at(node).contains(&(j as u32)));
            }
            // A rebuild with fewer seeds leaves nothing behind.
            seeds.reset(q, k, 0, 1);
            t.rebuild(&seeds);
            assert_eq!(t.node_count(), k + 1);
        }
    }
}
