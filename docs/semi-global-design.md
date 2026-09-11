# Semi-global alignment in POASTA's band-doubling engine

Design notes from a Claude conversation. Captures (1) how SPOA implements semi-global,
(2) how to port it to POASTA's band-doubling engine while minimizing bands and cells,
(3) where optimality is preserved vs. sacrificed, and (4) the pragmatic mode for the
"small query into a larger graph" use case.

---

## 1. SPOA's reference implementation

SPOA (rvaser/spoa) exposes three alignment types in
`include/spoa/alignment_engine.hpp`:

```cpp
enum class AlignmentType { kSW, kNW, kOV };
```

`kOV` (overlap) is asymmetric: free end-gaps on the **graph** side, gaps still
charged on the **query** side.

- **First row** `H[0][j]`: penalized like NW. Query-start gaps cost.
- **First column** `H[i][0]`: zeroed for *every* node (sisd_alignment_engine.cpp:230-253):

  ```cpp
  for (std::uint32_t i = 1; i < matrix_height; ++i) {
      pimpl_->H[i * matrix_width] = 0;
  }
  ```

  This is the "start anywhere in graph" mechanism.

- **Termination**: for each graph sink (`it->outedges.empty()`) the entire row
  is scanned for the best score; traceback then walks back until
  `i == 0 || j == 0` (the `ov_condition` lambda). So query can end early on a
  sink, but the graph must reach a sink.

- **No sentinels**: SPOA's graph has no `#`/`$` sentinels. Per-node handling is
  the source of correctness. SIMD and SISD paths share the algorithm.

---

## 2. Porting to POASTA's band-doubling engine

POASTA today seeds bands only at the start sentinel (band_doubling.rs:277-387)
and reads the final score from predecessors of the end sentinel
(band_doubling.rs:698-724). `AlignmentMode::EndsFree` is honored by the
canonical DP but **not** by the band-doubling engine.

### Naive port: seed everywhere
Push a free-entry band `[0, k+1)` onto every real node and initialize its
`qi=0` cell to 0. This faithfully reproduces SPOA's first-column zero.

Problem: if you let those seeds merge with the inherited diagonal band, band
widths grow linearly with depth (the seed anchors at `qlo=0`, the inherited
band extends right, and they overlap at `q=0`). Cells blow up.

### Minimum-band port
1. Seed `[0, k+1)` only at nodes whose inherited band wouldn't already include
   `q=0` (i.e., depth > k). Closer nodes get free entry automatically from the
   start sentinel.
2. **Do not** force-merge the free-entry seed with a disjoint inherited
   diagonal band. Keep them as two narrow bands (width ≈ k) instead of one wide
   band (width ≈ d+k). The change is at band_doubling.rs:323 — only merge when
   the inherited band genuinely touches `q=0`.
3. Add a kernel hook `init_free_entry` that zeroes the M state at `qi=0` for
   any band tagged free-entry.
4. Replace the end-predecessor termination scan with a scan over all bands on
   all real nodes (for `graph_free_end = Unbounded`), filtered by the mode's
   query-end and graph-end constraints.

Worst-case cells ≈ `O(|V| · k)`. That's the honest floor: every node that
*could* be an anchor needs a band touching `q=0`.

---

## 3. Cell-count optimizations

### What preserves optimality

| Technique | Optimal? | Why |
|---|---|---|
| Two-phase budget: global pre-pass gives `s_g`, semi-global runs with `k₂ = s_g / min_ge + 1` | ✅ | Semi-global score `s* ≤ s_g`. Ukkonen certificate catches under-estimates and triggers doubling. |
| Don't force-merge disjoint bands | ✅ | Same cells, just different bookkeeping. |
| Static band-construction filter using admissible `h(v,q)` | ✅ | Shrinks band widths without skipping any optimal cell. |
| Termination scan honoring the mode's freeness flags | ✅ | Just expands the set of acceptable end states. |
| Sparse free-entry seeding (k-mer / heuristic anchors) | ❌ | True anchor may live at an unseeded node. |
| Free-entry via init-rule only, no extra seeds at deep nodes | ❌ | Bands drift away from `q=0` past depth `k`; deep anchors unreachable. |

### A*-style live pruning has a real limitation

Because bands are processed in topological order and end-candidates live at the
tail of the topo order, `s_best` doesn't tighten during a single forward sweep
**except when `graph_free_end = Unbounded`** (every real node is then an
end-candidate, so every completed band tightens `s_best`).

Where pruning *does* pay:
1. **Across band-doubling iterations** — each attempt's achieved score becomes
   the next attempt's UB, collapsing the budget. This is the primary cell
   minimizer.
2. **`graph_free_end = Unbounded` only** — live `s_best` tightening works.
3. **Static band-construction filter** — using a precomputed reverse LB
   `h(v, q)` plus the best-predecessor `g` already in hand (predecessors are
   topologically earlier and therefore computed), prune individual `qi` values
   during band construction. Shrinks widths, not just band lists.

---

## 4. Admissible reverse bound `h(v, q)`

For a cell `(v, q)`, the cost to reach a valid end state is at least:

```
remaining_q = m - q
excess_q = max(0, remaining_q - max_hops_to_end[v])   // if !qry_free_end
excess_g = max(0, min_hops_to_end[v] - remaining_q)   // if !graph_free_end
h = (excess_q > 0) · gap_open_ins + excess_q · gap_extend_ins
  + (excess_g > 0) · gap_open_del + excess_g · gap_extend_del
```

`min_hops_to_end[v]` and `max_hops_to_end[v]` are computed in one reverse
topological pass over the DAG. For `graph_free_end = Unbounded`, both are zero
for every real node, so `excess_g = 0` and `h` reduces to a function of `q`
alone — precompute as a length-`m+1` table.

This bound is admissible (matches and non-surplus gaps contribute zero,
underestimating). Pair with an upper bound (Phase-1 `s_g`, or `s_best` if it
has tightened) for safe pruning.

---

## 5. Mode for "small query → larger graph"

In POASTA's `AlignmentMode::EndsFree`:

```rust
AlignmentMode::EndsFree {
    qry_free_begin:   Bound::Excluded(0),  // query starts at position 0
    qry_free_end:     Bound::Excluded(0),  // query ends at position m
    graph_free_begin: Bound::Unbounded,    // match can start anywhere in graph
    graph_free_end:   Bound::Unbounded,    // match can end anywhere in graph
}
```

This is stricter than SPOA's `kOV` (which is graph-start + query-end free, not
graph-end free). It's the right semantics for "find the best substring of the
graph that matches this whole query."

Consequences:
- `min_hops_to_end[v] = max_hops_to_end[v] = 0` for every real node.
- `h(v, q) = if m > q { gap_open_ins + (m - q) · gap_extend_ins } else { 0 }`.
- Every real node is an end candidate → `s_best` tightens live during the
  topological sweep, so A*-style pruning actually works.
- Termination scan picks the min over (node, qi=m) across all live bands;
  graph-end freeness means any node counts, query-end strictness means only
  `q = m` cells qualify.

---

## 6. Implementation checklist

1. Add `init_free_entry` to the kernel trait alongside `init_start`.
2. Add a `dead: bool` and/or `free_entry: bool` flag to `Band`.
3. Extend `Bands::for_global_alignment` → `Bands::for_semi_global`:
   - seed free-entry bands at nodes with depth > k,
   - skip the force-merge of disjoint inherited + free-entry bands.
4. Add `ReverseBounds` (precomputed `min_hops_to_end`, `max_hops_to_end`).
5. Use `ReverseBounds::h` during band construction to trim cell ranges.
6. In the forward DP loop, after computing each band, tighten `s_best` if the
   mode and node qualify as an end candidate (cheap check on
   `graph_free_end == Unbounded` or `is_sink`).
7. Replace the end-predecessor termination scan with a mode-aware scan over
   all live bands.
8. Drive the outer band-doubling loop with `k₂ = s_g / min_ge + 1` after a
   Phase-1 global pre-pass.

The Ukkonen optimality certificate (`required_k = forced + s/min_ge + 1`)
carries over without change; semi-global `s*` is only smaller than the global
`s_g`, so a certificate passing under global passes under semi-global too.

---

## References in the codebase

- `src/align/mod.rs:17-29` — `AlignmentMode` enum.
- `src/align/engine/band_doubling.rs:171-225` — `Band` struct.
- `src/align/engine/band_doubling.rs:242-251` — `BandEdge`.
- `src/align/engine/band_doubling.rs:277-387` — current global-only band
  seeding/propagation.
- `src/align/engine/band_doubling.rs:323` — the merge call site to gate on
  free-entry compatibility.
- `src/align/engine/band_doubling.rs:518-540` — `DynPredIter`.
- `src/align/engine/band_doubling.rs:615` — kernel `init_start` call site;
  add `init_free_entry` alongside.
- `src/align/engine/band_doubling.rs:655-695` — forward DP loop, where the
  live `s_best` update goes.
- `src/align/engine/band_doubling.rs:698-724` — end-predecessor termination
  scan to replace.
- `src/align/engine/dp.rs` — canonical DP, already mode-aware; the oracle for
  testing.
