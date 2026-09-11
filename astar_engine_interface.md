# A* Alignment Engine — Interface Design & Implementation Plan

Design for porting the A* aligner from `stable` into the `banded-poa` branch's
engine/trait architecture. This document is the review artifact requested by
`astar_engine_port.md`; no code has been written yet.

---

## 1. Goals & scope

Port the gap-affine (and 2-piece / linear) A* aligner so that it:

1. Implements the existing `AlignmentEngine<&[u8]>` trait and drops into
   `PoastaAligner<E, G>` with **zero changes to caller code** — the same slot
   currently filled by `CanonicalDP` and `BandDoublingEngineScalar`.
2. Supports **linear, affine, and 2-piece** scoring.
3. Supports **flexible heuristic configuration** (Dijkstra / gap-cost),
   parameterised by alignment mode (global vs semi-global).
4. Constructs and reuses a **bubble index** for the gap-cost heuristic and for
   the pruning ("reached bubble exits") optimisation.
5. Later unlocks **semi-global (ends-free)** alignment, which the band-doubling
   engine cannot do.

Non-goals for this port: SIMD kernels (A* is inherently scalar/priority-driven),
and changing the DP/banded engines.

---

## 2. How A* differs from the DP engines (why the interface is separate)

The banded/canonical engines are **kernel-driven**: `AlignmentCostModel` picks a
`DPKernel` that fills dense band columns. A* is **priority-queue-driven**: it
expands individual alignment-graph states in best-first order, greedily extends
matches depth-first (0-cost match edges), and terminates when the end state pops.

So A* needs its own internal traits (`AlignmentGraph`, `AstarVisited`,
`AstarQueue`, `AstarHeuristic`) that have **no DP-kernel analogue**. These stay
*private to the A* engine module*; the only public surface is the
`AlignmentEngine` impl. This mirrors how `stable` was organised, but re-homed
under `src/align/`.

Crucially, the `DPKernel` associated type on `AlignmentCostModel` is simply
**ignored** by A* — the cost structs (`Affine`, `TwoPiece`, `Linear`) already
carry the scalar cost getters A* needs.

---

## 3. Module layout

Bubble analysis lives under `src/graph/` (it is a property of the graph, reusable
beyond A*); the engine lives under `src/align/`.

```
src/
├── graph/
│   ├── mod.rs                   # + `pub mod bubbles;`
│   └── bubbles/                 # NEW — ported ~verbatim from stable
│       ├── mod.rs
│       ├── finder.rs            # SuperbubbleFinder (Gärtner et al.)
│       ├── index.rs             # BubbleIndex: node→bubble map, dist-to-end
│       └── reached.rs           # ReachedBubbleExitsMatch pruning helper
└── align/
    ├── mod.rs                   # + re-export astar engine; reuse AlignState/AlignmentMode
    ├── score.rs                 # NEW/shared — u32 + INF + sat/min3, promoted from engine/dp.rs
    ├── offsets.rs               # NEW — OffsetType (compact query-offset newtype)
    ├── cost_models/
    │   └── mod.rs               # + add provided `gap_cost(state, len)` method
    └── engine/
        ├── mod.rs               # + `pub mod astar;`
        ├── dp.rs                # INF/sat/min3 now re-exported from align::score
        └── astar/               # NEW — the engine
            ├── mod.rs           # AstarEngine<C, H>: AlignmentEngine impl + main loop
            ├── aln_graph.rs     # AlignmentGraphNode, AlignmentGraph trait + affine/2p impls
            ├── visited.rs       # AstarVisited trait + BlockedVisitedStorage ([u32; STATES] cells)
            ├── queue.rs         # LayeredQueue + per-model queue layers
            ├── dfa.rs           # DepthFirstGreedyAlignment (match extension)
            └── heuristic.rs     # AstarHeuristic trait + Dijkstra / MinGapCost
```

Rationale for `graph/bubbles/`: the bubble index is a property of the graph, is
independently testable, and is cached on the `POAGraph` across alignments (the
retained `with_existing_bubbles` reuse path).

---

## 4. Graph-trait bridging (the only change to existing traits)

`AlignableGraph` on `banded-poa` already supplies everything A* needs except a
symbol-equality helper. Add one **provided** method (no impl changes required):

```rust
// src/align/traits.rs
pub trait AlignableGraph: GraphWithNodeOrdering<NodeType = Self::Node> {
    type Node: GraphNodeId;
    fn get_node_symbol(&self, p: Self::Node) -> u8;

    /// NEW: provided. A* compares query bytes against graph symbols in its
    /// hottest loop; give impls a chance to specialise, default is trivial.
    #[inline]
    fn is_symbol_equal(&self, p: Self::Node, symbol: u8) -> bool {
        self.get_node_symbol(p) == symbol
    }
}
```

Name mapping used throughout the port (stable → banded-poa):

| stable                              | banded-poa                          |
| ----------------------------------- | ----------------------------------- |
| `AlignableRefGraph`                 | `AlignableGraph`                    |
| `NodeIndexType`                     | `GraphNodeId`                       |
| `G::NodeIndex`                      | `G::Node`                           |
| `node_count_with_start_and_end()`  | `node_count()` (already inc. start/end) |
| `get_node_ordering() -> Vec<usize>` | `node_rank(n)` / `rank_to_node(r)`  |
| `is_symbol_equal`                   | `is_symbol_equal` (added above)     |
| `graphs::tools::rev_postorder_nodes`| new helper in `graph` (see §9, Task 1) |
| `src/bubbles/`                      | `src/graph/bubbles/`                |

`SuccessorIterator`/`successors`/`predecessors` map 1:1
(`GraphBase::Successors`/`Predecessors`).

---

## 5. Cost-model extension

A* needs the "minimum gap cost from a state over a length" function. Add it as a
**provided method** on the existing trait so all three cost structs get it free:

```rust
// src/align/cost_models/mod.rs
pub trait AlignmentCostModel {
    type Kernel: DPKernel<Costs = Self>;   // unchanged; ignored by A*
    fn mismatch(&self) -> u8;
    fn gap_open(&self) -> u8;
    fn gap_extend(&self) -> u8;
    fn gap_open2(&self) -> u8;
    fn gap_extend2(&self) -> u8;
    fn min_gap_extend(&self) -> u8;

    /// NEW: cost of a gap of `length` given the current state.
    /// `length == 0` ⇒ 0. If already in a gap state, no open cost is charged.
    /// Two-piece models return the cheaper of the two pieces.
    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize;
}
```

Default provided impls can cover linear/affine; `TwoPiece` overrides to take the
min over both pieces. (On stable this lived on the `AlignmentCosts` trait; here it
folds into the one cost trait so there is a single source of cost parameters.)

---

## 6. Core A* types & traits

### 6.1 Offset type (`offsets.rs`)

Query offsets are stored in every visited cell and every queued state, so their
width dominates memory traffic. Port a trimmed `OffsetType`:

```rust
pub trait OffsetType: Copy + Ord + Hash + Debug + Default {
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn increase_one(&self) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
}
impl OffsetType for u32 { /* ... */ }   // u16/u64 also, chosen by seq length
```

The engine selects `O = u32` by default and can fall back to `u64` for very long
queries (dispatch in the `AlignmentEngine::align` entry point based on
`query.len()`). Trimmed vs stable: drop the `num`-crate supertraits
(`Unsigned + SaturatingAdd + …`) — A* only uses the four ops above, so the heavy
bounds are unnecessary.

### 6.2 Score — harmonized with the DP engines (`align/score.rs`)

**Decision: drop stable's `Score` enum / `NonMaxU32` and use the DP engines'
`u32 + INF` representation everywhere.** Promote the trio already living in
`engine/dp.rs` into a shared `align::score` module:

```rust
// src/align/score.rs  (moved from engine/dp.rs; dp.rs re-exports for compatibility)
pub const INF: u32 = 1 << 28;              // "unreachable / not yet reached"
#[inline(always)] pub fn sat(a: u32, b: u32) -> u32 { a.wrapping_add(b).min(INF) }
#[inline(always)] pub fn min3(a: u32, b: u32, c: u32) -> u32 { a.min(b).min(c) }
```

A* consequences:
- Visited cells default to `INF` (== stable's `Unvisited`); `update_score_if_lower`
  is a plain `new < *slot` on `u32`.
- Priority `f = g + h` stays `u32` (`sat(g, h_u32)`), no enum unwrapping.
- Backtrace equalities become exact `u32` arithmetic (`pred == curr - cost`),
  identical in spirit to the DP kernels' `backtrace_op`.

This removes the `nonmax` dependency from the A* path and means **one score
convention across all three engines** — see §11.1.

### 6.3 Alignment-graph state (`aln_graph.rs`)

```rust
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AlignmentGraphNode<N, O> { node: N, offset: O }   // (graph node, query offset)

pub trait AlignmentGraph {
    type CostModel: AlignmentCostModel;
    fn get_costs(&self) -> &Self::CostModel;

    fn initial_states<G, O>(&self, g: &G) -> Vec<AlignmentGraphNode<G::Node, O>>;
    fn is_end<G, O>(&self, g: &G, seq: &[u8],
                    n: &AlignmentGraphNode<G::Node, O>, s: AlignState) -> bool;

    // Expansion of successor states; `f` receives (score_delta, succ_node, succ_state)
    fn expand_all<V, G, O, F>(&self, v: &mut V, g: &G, seq: &[u8], score: Score,
                              n: &AlignmentGraphNode<G::Node, O>, s: AlignState, f: F);
    fn expand_ref_graph_end<V, N, O, F>(&self, v: &mut V, parent: &..., score, f);
    fn expand_query_end<V, N, O, F>(&self, v: &mut V, parent, child, score, f);
    fn expand_mismatch<V, N, O, F>(&self, v: &mut V, parent, child, score, f);
}
```

`AlignmentMode` is threaded into `initial_states`/`is_end`: `Global` seeds the
start sentinel at offset 0 and ends at (end node, `seq.len()`); `EndsFree` seeds
the extra free-indel start states and accepts early end states. **The `EndsFree`
arms are where semi-global support lands** — the band-doubling engine has no
equivalent.

### 6.4 Visited storage (`visited.rs`)

```rust
pub trait AstarVisited<N, O> {
    fn get_score(&self, n: &AlignmentGraphNode<N,O>, s: AlignState) -> u32;   // INF = unvisited
    fn set_score(&mut self, n: &AlignmentGraphNode<N,O>, s: AlignState, sc: u32);
    fn update_score_if_lower(&mut self, n, s, parent, parent_s, sc) -> bool;
    fn mark_reached(&mut self, sc, n, s);
    fn dfa_match(&mut self, sc, parent, child);
    fn prune(&self, sc, n, s) -> bool;
    fn backtrace<G>(&self, g: &G, seq: &[u8], end: &AlignmentGraphNode<N,O>)
        -> Vec<AlignedPair<N>>;                 // ← produces banded-poa's AlignedPair
    fn write_tsv<W: Write>(&self, w: &mut W) -> Result<(), PoastaError>;
}
```

Backing store `BlockedVisitedStorage<N, O, const S: usize, const B: usize = 8>`:
node-rank × query-offset space tiled into **B×B blocks** allocated on demand
(`Vec<FxHashMap<offset_block, [[Cell<S>; B]; B]>>`). This is the "BxB blocks at
once" allocation strategy called out in the plan — keeps a reached region's cells
contiguous and cache-local, and one hash lookup amortises over B² cells.

**Cell layout harmonized with DP columns:** a cell is `[u32; S]` where `S =
DPKernel::STATES` and slot `i` is `AlignState as u8` — the *exact* layout the DP
kernels already use for a band column (`[M, D, I]` for affine, verified in
`kernels/affine.rs`). So `get_score(n, state)` is `cell[state as usize]`, and the
same `state as u8` convention is shared by both engines (§11.2). Affine `S = 3`,
two-piece `S = 5`, linear `S = 1`.

### 6.5 Priority queue (`queue.rs`)

Port `LayeredQueue<L> + QueueLayer`: a bucket queue keyed by priority
`f = g + h`. Per-model layer keeps **separate vecs per align state** (M/I/D…) to
cut branch mispredicts when popping. LIFO within a bucket (stack-like) preserves
the depth-first flavour that pairs with the greedy extender.

### 6.6 Greedy match extension (`dfa.rs`)

Port `DepthFirstGreedyAlignment` ~verbatim. Because match edges cost 0, once a
Match state pops we DFS through matching (node, offset+1) successors for free,
only re-entering the queue at mismatches / graph-end / query-end. This is the
single biggest constant-factor win and must be preserved.

### 6.7 Heuristics (`heuristic.rs`)

```rust
pub trait AstarHeuristic<N, O> {
    fn h(&self, n: &AlignmentGraphNode<N, O>, s: AlignState) -> usize;
}
pub struct Dijkstra;                         // h ≡ 0  (exact, slow, always valid)
pub struct MinimumGapCostAffine<N, O> {      // admissible gap-cost lower bound
    costs: /* GapAffine-equivalent */, bubble_index: Arc<BubbleIndex<N>>, seq_len: usize,
}
```

`MinimumGapCostAffine` uses `BubbleIndex::get_{min,max}_dist_to_end` to bound the
minimum remaining indels and multiplies by gap cost. **Mode-awareness**: under
`EndsFree`, trailing/leading gaps are free, so the heuristic must return 0 for the
free dimension rather than charging gap cost — the heuristic struct is
constructed with the `AlignmentMode` and branches on it (this is the "alter
behaviour depending on alignment mode" requirement).

---

## 7. Public engine surface (`engine/astar/mod.rs`)

The whole machinery hides behind one struct implementing the existing trait:

```rust
pub struct AstarEngine<C, H = MinGapCostHeuristic> {
    costs: C,
    mode: AlignmentMode,
    heuristic_cfg: H,                 // zero-sized selector; see below
    existing_bubbles: Option<Arc<BubbleIndex<...>>>,   // optional reuse
    enable_pruning: bool,
}

impl<C, H> AstarEngine<C, H> {
    pub fn new(costs: C) -> Self;                      // global, gap-cost heuristic, pruning on
    pub fn with_mode(self, mode: AlignmentMode) -> Self;
    pub fn with_heuristic<H2>(self, h: H2) -> AstarEngine<C, H2>;
    pub fn with_existing_bubbles(self, idx: Arc<BubbleIndex<..>>) -> Self;
    pub fn without_pruning(self) -> Self;
}

impl<C, H, G> AlignmentEngine<&[u8]> for AstarEngine<C, H>
where C: AlignmentCostModel, H: HeuristicConfig, G: AlignableGraph
{
    type Graph = G;
    type Success = AlignOutput<G>;    // reuse existing result type → PoastaAligner works unchanged
    type Error = AstarError;          // DECIDED: fallible — empty queue before end ⇒ Err, not panic

    fn align(&self, graph: &G, query: &[u8]) -> Result<AlignOutput<G>, Self::Error> {
        // 1. offset-width dispatch (u32 vs u64) based on query.len()
        // 2. build/reuse BubbleIndex
        // 3. build AlignmentGraph (from costs + mode), visited storage, heuristic
        // 4. run astar_alignment(...) main loop  (ported from stable)
        // 5. backtrace → Vec<AlignedPair>, fill AlignmentStats, wrap in AlignOutput
    }
}

/// Returned when the queue empties before the end state is reached
/// (e.g. a disconnected graph). Implements `std::error::Error` so it
/// satisfies `AlignmentEngine::Error` and flows through `PoastaAlignerError`.
#[derive(Debug)]
pub enum AstarError { NoAlignmentFound }
```

### Heuristic configuration

Replace stable's four hand-written `AlignmentConfig` structs
(`AffineDijkstra`, `AffineMinGapCost`, `Affine2Piece*`, …) with a small
`HeuristicConfig` selector trait carried as the `H` type parameter. It maps
`(cost model, mode)` → a concrete `AstarHeuristic` and says whether a bubble
index is required:

```rust
pub trait HeuristicConfig {
    fn needs_bubble_index(&self) -> bool;
    fn build<C, N, O>(&self, costs: &C, mode: AlignmentMode,
                      bubbles: Option<Arc<BubbleIndex<N>>>, seq_len: usize)
        -> impl AstarHeuristic<N, O>;
}
struct Dijkstra;             // no bubble index
struct MinGapCost;           // needs bubble index; 2-piece uses cheaper piece
```

This collapses the 4× config-struct combinatorial explosion on stable into
`AstarEngine<CostModel, HeuristicSelector>`, and makes "affine + Dijkstra" vs
"2-piece + min-gap-cost" a two-type choice instead of a bespoke struct each.

---

## 8. Performance considerations (per the plan's Notes)

- **B×B block allocation** for visited storage (§6.4): amortised hashing,
  contiguous cache lines, allocate a block only when its region is reached.
- **Compact `OffsetType`** (§6.1): `u32` offsets halve state size vs `usize`,
  doubling states per cache line in both the queue and visited cells.
- **`u32 + INF` scores** (§6.2): a cell slot is a bare `u32` with `INF`
  encoding "unvisited" — no enum tag, no occupancy bitmap, same as the DP band.
- **Per-state queue vectors** (§6.5): pop M/I/D from separate vecs → predictable
  branches in the hottest loop.
- **Greedy DFA extension** (§6.6): 0-cost match runs never touch the queue.
- **`#[inline]` on `is_symbol_equal`, offset ops, block-index math** — all tiny,
  all hot.
- **Bubble pruning** gated by `enable_pruning` so it can be measured on/off.

---

## 9. Step-by-step implementation task list

Ordered so each step compiles and is testable before the next. Bubbles and
plumbing first (pure, unit-testable), engine last.

1. **Graph helper: reverse post-order.** Add `rev_postorder_nodes(&G) ->
   Vec<G::Node>` (and rank inverse) to the graph module, matching what the
   bubble finder consumes. Reuse the existing topological-sort DFS. *Test:*
   ordering on the mock graphs.

2. **Port bubbles into `src/graph/bubbles/`.** `finder.rs`, `index.rs`,
   `reached.rs`, `mod.rs` with the trait renames from §4. Wire `graph/mod.rs`
   (`pub mod bubbles;`). *Test:* re-enable stable's commented-out `BubbleIndex`
   unit tests against `graph::mock` graphs (dist-to-end, node→bubble map).

2a. **Promote `align::score`.** Move `INF`/`sat`/`min3` from `engine/dp.rs` to a
   new `align/score.rs`; re-export from `dp.rs` so the banded/canonical engines are
   untouched. This is the shared score representation both engines build on (§11.1).
   *Test:* existing DP/banded tests still pass (pure move).

3. **`offsets.rs`.** Trimmed `OffsetType` + `u32` impl (default). *Test:* trivial
   round-trip.

4. **`cost_models`: add `gap_cost`.** Provided method + `TwoPiece` override.
   *Test:* affine vs 2-piece gap-cost values (port stable's heuristic doctest
   numbers).

5. **`AlignableGraph::is_symbol_equal`.** Provided method (§4). No test needed.

6. *(folded into Task 2a — no separate A* `Score` type; use `u32 + INF`.)*

7. **`aln_graph.rs`.** `AlignmentGraphNode`, `AlignmentGraph` trait, and the
   **affine** impl (`initial_states`/`is_end`/`expand_*`) for `Global` mode
   only. Leave `EndsFree` arms as `todo!()`.

8. **`visited.rs`.** `AstarVisited` trait + `BlockedVisitedStorage` with
   `[u32; STATES]` cells, INF-initialised, `AlignState as u8` slot order (§11.2).
   *Test:* set/get/update-if-lower on a small grid.

9. **`queue.rs`.** `LayeredQueue`, `QueueLayer`, affine layer. *Test:*
   monotone pop order by priority.

10. **`dfa.rs`.** `DepthFirstGreedyAlignment`. (Compiles against traits from
    7–9; exercised via the engine.)

11. **`heuristic.rs`.** `AstarHeuristic`, `Dijkstra`, `MinimumGapCostAffine`
    (global mode). *Test:* port stable's `test_min_gap_cost`.

12. **`engine/astar/mod.rs`.** `AstarEngine`, `HeuristicConfig`, the main
    `astar_alignment` loop, backtrace → `AlignOutput`, `AlignmentStats`
    (`cells_computed` = states visited, etc.). Offset-width dispatch.
    *Test (the big one):* **equivalence to `CanonicalDP`** — same score on the
    existing test corpus (`tests/`), reusing the oracle pattern already used for
    the banded engine. Run Dijkstra and MinGapCost; both must match the oracle.

13. **2-piece + linear.** Add the 2-piece cell/queue-layer/expansion and the
    linear specialisation. *Test:* equivalence to `CanonicalDP` with the
    respective cost models.

14. **Semi-global (`EndsFree`).** Fill the `todo!()` arms in `aln_graph.rs` and
    make the heuristic mode-aware (§6.7). *Test:* semi-global cases (no DP oracle
    on `banded-poa` for this yet — validate against hand-checked cases or the
    `stable` engine's output).

15. **CLI wiring.** Expose engine choice (`--engine astar`) and heuristic flag in
    `src/cli/poasta.rs`; default remains the band-doubling engine.

16. **Bench & tune.** `poasta_vs_spoa` / existing stats plumbing; sweep block
    size `B`, offset width, pruning on/off. Confirm no regressions on the
    `AlignmentStats` fraction-of-matrix metric.

### Suggested milestones
- **M1 (Tasks 1–6):** bubbles + primitives landed and unit-tested.
- **M2 (Tasks 7–12):** affine global A* matches the canonical DP oracle. ← core.
- **M3 (Tasks 13–14):** all cost models + semi-global.
- **M4 (Tasks 15–16):** CLI + performance tuning.

---

## 10. Resolved decisions

1. **`Error` type → fallible.** Empty queue before the end state returns
   `Err(AstarError::NoAlignmentFound)`; no `panic!`. (§7)
2. **Bubble-index reuse → keep the reuse path.** `AstarEngine` retains
   `with_existing_bubbles(Arc<BubbleIndex>)`; the `POAGraph` may cache and hand in
   an index so it is not rebuilt every alignment. `BubbleIndex` stays `Arc`-shared
   between the visited-storage pruning and the heuristic.
3. **Offset width → `u32` default.** `align` uses `O = u32`; a `u64` fallback for
   >~4 Gbp queries can be added later but is not built now.
4. **No feature-gating.** `score`, `offsets`, and `graph::bubbles` are always
   compiled — they are cheap and the shared `AlignmentEngine` surface is the goal.

---

## 11. Harmonizing the A* and DP/banded engines

The two engines are algorithmically different (best-first priority search vs.
dense band-column DP), so they should *not* share their inner loops. But they
model the **same alignment problem with the same costs**, and today they diverge
in several representations that need not differ. Unifying those shrinks the
surface area, lets tests and debug tooling target both engines, and is the
natural follow-on to the score/state harmonization already folded in above.

Ordered from "already decided / low-risk" to "worth doing but bigger".

### 11.1 Score representation — **decided, do it first**
Promote `INF`/`sat`/`min3` out of `engine/dp.rs` into `align::score` and have all
three engines use `u32` scores with `INF` as the sentinel (replacing stable's
`Score`/`NonMaxU32` enum). One arithmetic convention, one place to change the
overflow bound, and A* backtrace equalities read like the DP `backtrace_op`
comparisons. (Tasks 2a, 6.)

### 11.2 `AlignState ⇄ state-id ⇄ cell slot` — **low-risk, high payoff**
The DP kernels already treat state ids as `AlignState as u8` and lay a band
column out as `[u32; STATES]` with slot order `[M, D, I(, D2, I2)]`
(`kernels/affine.rs`). Have A*'s visited cell use the *identical* `[u32; STATES]`
layout and the same `state as usize` indexing (§6.4). Benefits:
- `AstarVisited::get_score(n, state)` becomes `cell[state as usize]` — no per-state
  `match`, fewer branches in the hot loop.
- `STATES` / `STATE_NAMES` come from **one** source (the kernel/model), so adding a
  scoring model updates both engines' state counts together.
- The debug-TSV writers (`write_tsv`) can share a single column-formatting helper
  keyed by `STATE_NAMES`.

### 11.3 A single "scoring model descriptor" both engines consume
Right now `AlignmentCostModel::Kernel` is a DP-only association. Factor the *model
shape* (state count, state names, `AlignState`→slot map, per-transition costs)
into a small trait — call it `ScoringScheme` — that `DPKernel` builds on and the
A* `AlignmentGraph` also consumes:

```rust
pub trait ScoringScheme {
    const STATES: usize;
    const STATE_NAMES: &'static [&'static str];
    type Costs: AlignmentCostModel;
    // transition-cost helpers used by BOTH engines, sourced from Costs:
    fn open_extend_del(costs: &Self::Costs) -> u32;   // gap_open + gap_extend
    fn extend_del(costs: &Self::Costs) -> u32;
    // …insertion / second-piece equivalents…
}
```

This removes the duplicated "gap_open + gap_extend" arithmetic scattered through
stable's `expand_*` and the DP kernels, and guarantees the two engines charge the
*same* transition costs by construction — which is exactly what the
CanonicalDP-equivalence test asserts at runtime. (Introduce in M3, once the affine
A* is proven, so it is refactor-under-green rather than up-front abstraction.)

### 11.4 Cost model is already the shared source of truth — keep it that way
`AlignmentCostModel` (mismatch/gap_open/gap_extend/gap_open2/gap_extend2/
min_gap_extend + the new `gap_cost`, §5) is the one type both engines read costs
from. No engine should hard-code costs; the only addition is `gap_cost`, needed by
the heuristic. This is already harmonized — the note is to *not* regress it by
adding an A*-private cost struct (stable had `GapAffine` separate from the DP
costs; here they collapse into one).

### 11.5 Result & stats — already unified, tighten the semantics
Both engines emit `AlignOutput<G>` = `{ score, Vec<AlignedPair>, AlignmentStats }`,
so `PoastaAligner`, GFA/FASTA output, and consensus all work unchanged. Two small
alignments to make:
- **`AlignmentStats` for A*:** `cells_computed` = alignment states visited
  (`num_visited`), which keeps `fraction_of_full_matrix` meaningful and directly
  comparable to the DP/banded number. `max_bandwidth` has no A* analogue —
  document it as `0`/not-applicable rather than inventing a value.
- **Backtrace decisions:** A*'s backtrace makes the same four moves the DP
  `BacktraceOp` enum encodes (Diagonal/Delete/Insert/SwitchState). Reuse that enum
  (or at least the `AlignedPair`-emitting tail) so both backtraces share the
  indel-vs-match push logic and the "skip zero-cost indel-close" rule.

### 11.6 `AlignmentMode` — one mode enum, honored by both
`align::AlignmentMode` already exists and both engines take it. Give all engine
constructors the same shape (`new(costs)` = Global; `.with_mode(mode)` builder) so
callers and the CLI switch engines without special-casing. A* implements
`EndsFree`; the banded engine returns a clear "unsupported mode" error instead of
silently doing Global — so the CLI can fall back to A* for semi-global.

### 11.7 Debug/trace tooling
`AstarVisited::write_tsv` and the DP debug writer both dump `(node, offset, state,
score)` rows. With 11.1+11.2 in place they can share one TSV schema and one
`Score`-formatting path, so the existing Python plotting scripts work against
either engine's output unchanged.

### What is *not* worth unifying
The forward computation itself: SIMD dense column folds (`accumulate_pred`/
`finalize`) vs. scalar best-first expansion + greedy DFA are fundamentally
different traversals. Forcing a shared abstraction over both would pessimize the
DP kernels' vectorization for no correctness benefit. Keep the inner loops
separate; share only the *vocabulary* (score, state ids, costs, results, mode).

### Recommended sequencing
- **Now (part of the port):** 11.1, 11.2, 11.4, 11.6 — cheap, and the A* code is
  written against the shared vocabulary from day one.
- **After affine A* is green (M3):** 11.3, 11.5, 11.7 — refactors best done with
  the equivalence test already guarding behaviour.
