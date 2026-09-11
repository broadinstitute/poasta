# CLAUDE.md

Guidance for Claude Code when working in this repository. The full
walkthrough of the code lives in [docs/architecture.md](docs/architecture.md);
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
cargo run --features bench --bin poasta_vs_spoa -- <args>   # needs spoa C++ build
cargo doc --no-deps
vale src docs CLAUDE.md README.md             # prose lint (comments only in .rs)
```

Binaries after a release build: `target/release/poasta`, `target/release/lasagna`.

## Map

- `src/align/`: aligning a query to a graph. `mod.rs` has `PoastaAligner`,
  `AlignmentMode` and `AlignState`; `traits.rs` has `AlignableGraph`,
  `AlignmentEngine` and `AlignOutput`; `cost_models/` the three cost models.
- `src/align/astar/`: the engine. `mod.rs` (`AstarEngine`, `AstarCostModel`,
  main loop), `aln_graph.rs` (state machines), `visited.rs` (scores and
  backtrace), `queue.rs` (bucket queue), `dfa.rs` (free match extension),
  `heuristic.rs`, `pruning.rs`, and the test-only `dp_oracle.rs`.
- `src/graph/`: `POAGraph` (`poa.rs`), graph traits, `AddAlignment`,
  consensus, superbubbles (`bubbles/`), and I/O (`io/`: FASTA, GFA, GAF, DOT).
- `src/cli/` and `src/bin/`: clap structs and the `poasta`, `lasagna` and
  `poasta_vs_spoa` drivers.
- `src/debug/`: `--debug-output-dir` support, a `tracing` layer that writes
  per-sequence TSVs of the search.

## Invariants

- The graph is read-only during alignment; `AddAlignment` runs afterwards.
- `graph/` never imports from `align/` except `poa.rs` implementing
  `AddAlignment` for `AlignOutput`.
- Every A* configuration must match `DpOracle`'s score on the test graphs,
  and every backtrace must re-score to the reported score
  (`test_support::score_from_alignment`).
- Query offsets are `u32`. Visited slots and queue vectors are indexed by
  `AlignState as usize`.
- Queue pop order within a bucket is LIFO, Match before deletions before
  insertions. It decides which equal-score alignment is returned, so
  changing it changes MSA output.

## Style

- Comments are prose: no em dashes, short sentences, no filler. Run vale.
- Deliberate shortcuts carry a `ponytail:` comment naming the ceiling.
- Keep the module docs (`//!`) in each `mod.rs` accurate when moving code.
