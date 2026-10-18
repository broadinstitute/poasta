//! Finding every seed match in the graph, per alignment.
//!
//! A match of seed j is a path of k graph nodes that spells the seed with
//! at most r substitutions. The index must be complete: a missed match
//! makes `h` charge for a seed the alignment can match for free, and `h`
//! overestimates. A full graph k-mer index does not fit, because the number
//! of k-node paths grows with bubble density (50 per node on a dense graph,
//! 1,487 at the worst node). So the finder is driven by the query instead:
//! one sweep over the ranks keeps, per node, the trie nodes of the seed
//! prefixes that end there. A join merges its predecessors' lists and drops
//! duplicates, so a bubble costs the union of its prefixes, not the product
//! of its paths.

use std::ops::Range;

use super::seeds::{Seeds, Walk};
use super::trie::{SeedTrie, TrieNode};
use crate::align::traits::AlignableGraph;
use crate::graph::traits::GraphNodeId;

/// A partial seed match ending at a graph node: a trie node and the
/// substitutions spent so far.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Partial {
    trie_node: TrieNode,
    mismatches: u8,
}

/// The k-node path from `start` to `end` that spells a seed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct SeedMatch<N> {
    pub start: N,
    pub end: N,
}

/// A match of a seed, as the finder collects them.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Hit<N> {
    seed: u32,
    m: SeedMatch<N>,
}

/// Every match, grouped by seed, without duplicates. A match's id is its
/// index: ids run seed by seed, and within a seed by start node, then end
/// node.
#[derive(Clone, Debug)]
pub(crate) struct SeedMatches<N> {
    matches: Vec<SeedMatch<N>>,
    /// Seed j has `matches[seed_start[j]..seed_start[j + 1]]`.
    seed_start: Vec<u32>,
}

impl<N> Default for SeedMatches<N> {
    fn default() -> Self {
        Self {
            matches: Vec::new(),
            seed_start: Vec::new(),
        }
    }
}

impl<N: GraphNodeId> SeedMatches<N> {
    pub fn of_seed(&self, j: usize) -> &[SeedMatch<N>] {
        &self.matches[self.seed_start[j] as usize..self.seed_start[j + 1] as usize]
    }

    /// Ids of the matches of seed j that start at `node`.
    pub fn starting_at(&self, j: usize, node: N) -> Range<u32> {
        let of_seed = self.of_seed(j);
        let lo = of_seed.partition_point(|m| m.start.index() < node.index());
        let hi = lo + of_seed[lo..].partition_point(|m| m.start == node);
        self.seed_start[j] + lo as u32..self.seed_start[j] + hi as u32
    }

    pub fn len(&self) -> usize {
        self.matches.len()
    }
}

/// Keeps its buffers between alignments.
#[derive(Clone, Debug)]
pub(crate) struct SeedFinder<N> {
    /// Partial matches ending at each rank: `live[live_start[rank]..live_start[rank + 1]]`.
    live: Vec<Partial>,
    live_start: Vec<u32>,
    scratch: Vec<Partial>,
    hits: Vec<Hit<N>>,
    walk: Vec<Walk<N>>,
}

impl<N> Default for SeedFinder<N> {
    fn default() -> Self {
        Self {
            live: Vec::new(),
            live_start: Vec::new(),
            scratch: Vec::new(),
            hits: Vec::new(),
            walk: Vec::new(),
        }
    }
}

impl<N: GraphNodeId> SeedFinder<N> {
    /// Find every match of every seed in `trie`, which holds `seeds`.
    pub fn find<G>(&mut self, graph: &G, seeds: &Seeds, trie: &SeedTrie, out: &mut SeedMatches<N>)
    where
        G: AlignableGraph<Node = N>,
    {
        let (start, end) = (graph.start_node(), graph.end_node());
        let n = graph.node_count();
        self.live.clear();
        self.live_start.clear();
        self.live_start.push(0);
        self.hits.clear();

        for rank in 0..n {
            let v = graph.rank_to_node(rank);
            if v != start && v != end {
                // A fresh prefix at v, and every prefix that ends at a
                // predecessor, each extended by v's symbol.
                self.scratch.clear();
                extend(
                    graph,
                    trie,
                    v,
                    Partial {
                        trie_node: 0,
                        mismatches: 0,
                    },
                    seeds.r,
                    &mut self.scratch,
                );
                for p in graph.predecessors(v) {
                    if p == start {
                        continue;
                    }
                    let pr = graph.node_rank(p);
                    for ix in self.live_start[pr] as usize..self.live_start[pr + 1] as usize {
                        extend(graph, trie, v, self.live[ix], seeds.r, &mut self.scratch);
                    }
                }
                // One entry per trie node, with the fewest substitutions.
                self.scratch.sort_unstable();
                self.scratch.dedup_by_key(|p| p.trie_node);
                for &p in &self.scratch {
                    if trie.is_leaf(p.trie_node) {
                        for &j in trie.seeds_at(p.trie_node) {
                            walk_back(graph, seeds, j, v, &mut self.walk, &mut self.hits);
                        }
                    } else {
                        self.live.push(p);
                    }
                }
            }
            self.live_start.push(self.live.len() as u32);
        }

        self.hits
            .sort_unstable_by_key(|h| (h.seed, h.m.start.index(), h.m.end.index()));
        self.hits.dedup();

        out.matches.clear();
        out.matches.extend(self.hits.iter().map(|h| h.m));
        out.seed_start.clear();
        out.seed_start.resize(seeds.n + 1, 0);
        for h in &self.hits {
            out.seed_start[h.seed as usize + 1] += 1;
        }
        for j in 0..seeds.n {
            out.seed_start[j + 1] += out.seed_start[j];
        }
    }
}

/// Every trie child of `p` that node v can extend, within the budget.
fn extend<G: AlignableGraph>(
    graph: &G,
    trie: &SeedTrie,
    v: G::Node,
    p: Partial,
    r: usize,
    out: &mut Vec<Partial>,
) {
    for e in trie.children(p.trie_node) {
        if graph.is_symbol_equal(v, e.byte) {
            out.push(Partial {
                trie_node: e.child,
                mismatches: p.mismatches,
            });
        } else if (p.mismatches as usize) < r {
            out.push(Partial {
                trie_node: e.child,
                mismatches: p.mismatches + 1,
            });
        }
    }
}

/// Seed j ends at `end`. Walk back k − 1 steps to the start nodes of the
/// paths that spell it, one hit each.
fn walk_back<G: AlignableGraph>(
    graph: &G,
    seeds: &Seeds,
    j: u32,
    end: G::Node,
    walk: &mut Vec<Walk<G::Node>>,
    hits: &mut Vec<Hit<G::Node>>,
) {
    let seed = seeds.seed(j as usize);
    let k = seed.len();
    walk.clear();
    walk.push(Walk {
        node: end,
        pos: k - 1,
        mismatches: usize::from(!graph.is_symbol_equal(end, seed[k - 1])),
    });
    while let Some(w) = walk.pop() {
        if w.mismatches > seeds.r {
            continue;
        }
        if w.pos == 0 {
            hits.push(Hit {
                seed: j,
                m: SeedMatch { start: w.node, end },
            });
            continue;
        }
        for p in graph.predecessors(w.node) {
            if p != graph.start_node() {
                let mismatches =
                    w.mismatches + usize::from(!graph.is_symbol_equal(p, seed[w.pos - 1]));
                walk.push(Walk {
                    node: p,
                    pos: w.pos - 1,
                    mismatches,
                });
            }
        }
    }
}

/// Whether some path from a successor of `node` spells the rest of the
/// seed that offset i is inside, with at most r substitutions. `walk` is
/// scratch space.
pub(crate) fn rest_matches<G: AlignableGraph>(
    graph: &G,
    seeds: &Seeds,
    node: G::Node,
    i: usize,
    walk: &mut Vec<Walk<G::Node>>,
) -> bool {
    let rest = seeds.rest(i);
    walk.clear();
    walk.push(Walk {
        node,
        pos: 0,
        mismatches: 0,
    });
    while let Some(w) = walk.pop() {
        for x in graph.successors(w.node) {
            if x == graph.end_node() {
                continue;
            }
            let mismatches = w.mismatches + usize::from(!graph.is_symbol_equal(x, rest[w.pos]));
            if mismatches <= seeds.r {
                if w.pos + 1 == rest.len() {
                    return true;
                }
                walk.push(Walk {
                    node: x,
                    pos: w.pos + 1,
                    mismatches,
                });
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::PoastaAligner;
    use crate::align::astar::{AstarEngine, DijkstraHeuristic};
    use crate::align::cost_models::affine::Affine;
    use crate::graph::poa::{POAGraph, POANodeIndex};
    use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};

    type G = POAGraph<u32>;
    type N = POANodeIndex<u32>;

    fn rng(state: &mut u64, n: usize) -> usize {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        (*state % n as u64) as usize
    }

    fn random_graph(state: &mut u64) -> (G, Vec<u8>) {
        let alphabet = b"ACGT";
        let template: Vec<u8> = (0..10 + rng(state, 30))
            .map(|_| alphabet[rng(state, 4)])
            .collect();
        let mut aligner = PoastaAligner::new(
            AstarEngine::<G, Affine>::new(Affine::new(4, 6, 2)).with_heuristic(DijkstraHeuristic),
            POAGraph::new(),
        );
        for s in 0..1 + rng(state, 5) {
            let seq: Vec<u8> = template
                .iter()
                .flat_map(|&b| match rng(state, 12) {
                    0 => vec![],
                    1 => vec![b, alphabet[rng(state, 4)]],
                    2 => vec![alphabet[rng(state, 4)]],
                    _ => vec![b],
                })
                .collect();
            if !seq.is_empty() {
                aligner.align_one_named(&format!("s{s}"), &seq).unwrap();
            }
        }
        let query: Vec<u8> = template
            .iter()
            .map(|&b| {
                if rng(state, 10) == 0 {
                    alphabet[rng(state, 4)]
                } else {
                    b
                }
            })
            .collect();
        (aligner.into_graph(), query)
    }

    /// Every path of k real nodes.
    fn all_paths<A: AlignableGraph>(g: &A, k: usize) -> Vec<Vec<A::Node>> {
        let mut out = Vec::new();
        for v in g.all_nodes_iter() {
            if v == g.start_node() || v == g.end_node() {
                continue;
            }
            let mut stack = vec![vec![v]];
            while let Some(path) = stack.pop() {
                if path.len() == k {
                    out.push(path);
                    continue;
                }
                for s in g.successors(*path.last().unwrap()) {
                    if s != g.end_node() {
                        let mut p = path.clone();
                        p.push(s);
                        stack.push(p);
                    }
                }
            }
        }
        out
    }

    fn brute<A: AlignableGraph>(g: &A, q: &[u8], k: usize, r: usize) -> Vec<(u32, usize, usize)> {
        let paths = all_paths(g, k);
        let mut out = Vec::new();
        for j in 0..q.len() / k {
            let seed = &q[j * k..(j + 1) * k];
            for path in &paths {
                let mm = path
                    .iter()
                    .zip(seed)
                    .filter(|&(&v, &b)| !g.is_symbol_equal(v, b))
                    .count();
                if mm <= r {
                    out.push((j as u32, path[0].index(), path[k - 1].index()));
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn seeds(q: &[u8], k: usize, r: usize) -> Seeds {
        let mut seeds = Seeds::default();
        seeds.reset(q, k, r, q.len() / k);
        seeds
    }

    fn found<A: AlignableGraph>(g: &A, q: &[u8], k: usize, r: usize) -> Vec<(u32, usize, usize)> {
        let seeds = seeds(q, k, r);
        let mut trie = SeedTrie::default();
        trie.rebuild(&seeds);
        let mut finder = SeedFinder::default();
        let mut matches = SeedMatches::default();
        finder.find(g, &seeds, &trie, &mut matches);
        let mut out = Vec::new();
        for j in 0..seeds.n {
            out.extend(
                matches
                    .of_seed(j)
                    .iter()
                    .map(|m| (j as u32, m.start.index(), m.end.index())),
            );
            // `starting_at` finds each match by its start node.
            for (id, m) in (matches.seed_start[j]..).zip(matches.of_seed(j)) {
                assert!(matches.starting_at(j, m.start).contains(&id));
            }
        }
        out
    }

    #[test]
    fn finder_matches_brute_force() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        for _ in 0..200 {
            let (g, q) = random_graph(&mut state);
            for (k, r) in [(3, 0), (4, 0), (4, 1), (5, 1)] {
                assert_eq!(found(&g, &q, k, r), brute(&g, &q, k, r), "k={k} r={r}");
            }
        }
    }

    #[test]
    fn rest_matches_agrees_with_brute_force() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        for _ in 0..100 {
            let (g, q) = random_graph(&mut state);
            let k = 4;
            let seeds = seeds(&q, k, 1);
            for v in g.all_nodes_iter() {
                if v == g.end_node() {
                    continue;
                }
                for i in (0..q.len() / k * k).filter(|i| i % k != 0) {
                    let end = (i / k + 1) * k;
                    let want = all_paths(&g, end - i).iter().any(|path| {
                        g.successors(v).any(|x| x == path[0])
                            && path
                                .iter()
                                .zip(&q[i..end])
                                .filter(|&(&w, &b)| !g.is_symbol_equal(w, b))
                                .count()
                                <= 1
                    });
                    assert_eq!(rest_matches(&g, &seeds, v, i, &mut Vec::new()), want);
                }
            }
        }
    }

    /// A graph whose `N` nodes equal every base, like an IUPAC-aware
    /// graph. Byte comparison would miss matches through `N`.
    struct Iupac(G);

    impl GraphBase for Iupac {
        type NodeType = N;
        type NodeIter<'a> = <G as GraphBase>::NodeIter<'a>;
        type Successors<'a> = <G as GraphBase>::Successors<'a>;
        type Predecessors<'a> = <G as GraphBase>::Predecessors<'a>;
        fn all_nodes_iter(&self) -> Self::NodeIter<'_> {
            self.0.all_nodes_iter()
        }
        fn node_count(&self) -> usize {
            self.0.node_count()
        }
        fn is_empty(&self) -> bool {
            self.0.is_empty()
        }
        fn successors(&self, node: N) -> Self::Successors<'_> {
            self.0.successors(node)
        }
        fn predecessors(&self, node: N) -> Self::Predecessors<'_> {
            self.0.predecessors(node)
        }
        fn out_degree(&self, node: N) -> usize {
            self.0.out_degree(node)
        }
        fn in_degree(&self, node: N) -> usize {
            self.0.in_degree(node)
        }
    }

    impl GraphWithNodeOrdering for Iupac {
        fn start_node(&self) -> N {
            self.0.start_node()
        }
        fn end_node(&self) -> N {
            self.0.end_node()
        }
        fn node_rank(&self, node: N) -> usize {
            self.0.node_rank(node)
        }
        fn rank_to_node(&self, rank: usize) -> N {
            self.0.rank_to_node(rank)
        }
    }

    impl AlignableGraph for Iupac {
        type Node = N;
        fn get_node_symbol(&self, p: N) -> u8 {
            self.0.get_node_symbol(p)
        }
        fn is_symbol_equal(&self, p: N, symbol: u8) -> bool {
            let s = self.0.get_node_symbol(p);
            s == b'N' || s == symbol
        }
    }

    #[test]
    fn matches_go_through_nodes_that_equal_every_base() {
        use crate::graph::alignment::AddAlignment;
        let mut g = POAGraph::<u32>::new();
        g.add_alignment("s", b"ACNTTTT", None, &[1; 7]).unwrap();
        let g = Iupac(g);
        // Seed ACGT only matches through the N.
        let q = b"ACGT";
        let got = found(&g, q, 4, 0);
        assert_eq!(got.len(), 1);
        assert_eq!(got, brute(&g, q, 4, 0));
        // Sentinels are never part of a match.
        for (_, s, e) in got {
            for v in [s, e] {
                assert!(v != g.start_node().index() && v != g.end_node().index());
            }
        }
    }
}
