# CLAUDE.md

Guidance for Claude Code when working in this repository. The full
walkthrough of the code is in [docs/architecture.md](docs/architecture.md);
this file is the short map.

## What POASTA is

A partial order aligner in Rust. Sequences are aligned to a POA graph (a DAG
of symbols) with an A* search that supports linear, affine and two-piece
gap costs, global and ends-free modes, and a bubble-based heuristic and
pruning scheme. The result is folded back into the graph.

## Build and test

```bash
cargo build                                   # debug
RUSTFLAGS="-C target-cpu=native" cargo build --release
cargo test                                    # all unit tests (in-tree)
cargo test <name>                             # one test
cargo test --lib astar::                      # one module
cargo build --features disable_tracing        # trace events compiled out
cargo build --features viz                    # chain heuristic inspection API (poasta-viz)
cargo run --features bench --bin poasta_vs_spoa -- <args>   # needs spoa C++ build
cargo doc --no-deps
vale src docs CLAUDE.md                       # prose lint (comments only in .rs)
```

Binaries after a release build: `target/release/poasta`, `target/release/lasagna`.

## Map

- `src/align/`: aligning a query to a graph. `mod.rs` has `PoastaAligner`,
  `AlignmentMode` and `AlignState`; `traits.rs` has `AlignableGraph`,
  `AlignmentEngine` and `AlignOutput`. `cost_models/` has the three cost models.
- `src/align/astar/`: the engine. `mod.rs` (`AstarEngine`, `AstarCostModel`,
  main loop), `aln_graph.rs` (state machines), `visited.rs` (block storage,
  reached bits and backtrace), `queue.rs` (bucket queue), `dfa.rs` (free match extension),
  `heuristic/` (`mod.rs`: `AstarHeuristic`, `HeuristicConfig`, Dijkstra;
  `min_gap.rs`; `chain/`: the gap-chaining seed heuristic, see
  `docs/chain-heuristic-explained.md`; `chain/inspect.rs` is the `viz`
  inspection API), `pruning.rs`, and the test-only
  `dp_oracle.rs` and `differential_tests.rs`.
- `src/graph/`: `POAGraph` (`poa.rs`), graph traits, `AddAlignment`,
  consensus, superbubbles (`bubbles/`), and I/O (`io/`: FASTA, GFA, GAF, DOT).
- `src/cli/` and `src/bin/`: clap structs and the `poasta`, `lasagna` and
  `poasta_vs_spoa` drivers.
- `src/debug/`: `--debug-output-dir` support, a `tracing` layer that writes
  per-sequence TSVs of the search.

## Invariants

- The graph is read-only during alignment. `AddAlignment` runs afterwards.
- `graph/` never imports from `align/` except `poa.rs` implementing
  `AddAlignment` for `AlignOutput`.
- Every A* configuration, global and semi-global, must match `DpOracle`'s
  score on the test graphs, and every backtrace must re-score to the
  reported score (`test_support::score_from_alignment`). The randomized
  test (`differential_tests.rs`) checks this on random graphs and fails on
  any consistency violation of the heuristics.
- Symbol comparisons go through `AlignableGraph::is_symbol_equal`, never
  byte equality. A seed match missed that way makes the chain heuristic
  overestimate.
- An `AlignableGraph` node's index is its rank. The visited storage is
  indexed by it. `POAGraph` nodes are ranks, `NodeId` names store nodes.
- Query offsets are a generic `O: OffsetType` (`u16`, `u32`, `u64`) on
  `AstarEngine`. The CLI uses the default `Offset = u32`. Visited slots and
  queue vectors are indexed by `AlignState as usize`.
- `AlignmentEngine::align` takes `&mut self` so the engine can keep scratch
  storage between alignments. Every thread owns its engine.
- Queue pop order within a bucket is LIFO, Match before deletions before
  insertions. It decides which equal-score alignment is returned, so
  changing it changes MSA output.

## Style

- Comments are prose. Use short sentences and avoid em dashes and filler. Run vale.
- Deliberate shortcuts carry a comment that states their ceiling.
- Keep the module docs (`//!`) in each `mod.rs` accurate when moving code.
