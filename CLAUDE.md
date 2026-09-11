# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**POASTA** is a fast, gap-affine sequence-to-graph and partial order aligner (POA) written in Rust. It implements an Ukkonen-style band-doubling alignment for DNA/protein sequences against partial order alignment graphs (DAGs).

## Build & Test Commands

```bash
# Debug build
cargo build

# Optimized release build (native CPU)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Run all tests
cargo test

# Run a single test by name
cargo test <test_name>

# Run tests in a specific module
cargo test --lib graph::

# Run benchmarks against spoa (requires spoa-rs)
cargo run --features test --bin poasta_vs_spoa -- <args>
```

The CLI binary is built by default (via the `cli` feature). After release build: `./target/release/poasta align -h`.

## Architecture

### Core Modules

**`src/graph/`** — POA graph data structure
- `poa.rs`: `POAGraph<Ix>` — a `StableDiGraph` (petgraph) where each node is a byte (sequence symbol). Has sentinel `start` (#) and `end` ($) nodes. Nodes carry topological rank; edges carry weight and a set of sequence IDs traversing them.
- `traits.rs`: Trait hierarchy — `GraphBase` → `GraphWithNodeOrdering` → `AlignableGraph`. Alignment engines are generic over `AlignableGraph`.
- `alignment.rs`: `AddAlignment` trait — consumed by `POAGraph` to ingest an `AlignResult` and update the graph.
- `utils.rs`: Topological sort via reverse post-order DFS.

**`src/align/`** — Alignment engines (Gotoh 3-state: Match/Deletion/Insertion)
- `engine/dp.rs`: `CanonicalDP` — O(N×m) reference implementation with full DP matrix. Used as ground-truth oracle in tests.
- `engine/band_doubling.rs`: `BandDoublingEngineScalar` — Ukkonen-style per-node multi-band tracking with band doubling and merging.
- `cost_models/affine.rs`: `AffineAlignmentCosts` — standard affine gap model (match, mismatch, gap_open, gap_extend).
- `traits.rs`: `AlignmentEngine<ToAlign>` trait, `AlignableGraph` requirements, `AlignedPair<N>` type, `AlignmentMode` (Global / EndsFree).

**`src/cli/poasta.rs`** — Clap-derived CLI with subcommands: `align`, `view`, `stats`. Alignment subcommand accepts FASTA input, optional existing graph (`-g`), and outputs `.poasta` binary format or DOT/FASTA/GFA.

**`src/errors.rs`** — `GraphError<Ix>`, `PoastaIOError`, `PoastaError<Ix>`.

**`src/debug/`** — Conditional tracing for intermediate aligner state (only active in non-optimized builds).

### Key Design Patterns

- **Engines are generic**: `AlignmentEngine<G: AlignableGraph>` — cost models and graph types are pluggable.
- **Canonical DP is the oracle**: All optimized engines should produce results equivalent to `CanonicalDP` on the same inputs. Tests validate this.
- **Band doubling propagation**: Bands are tracked per graph node, seeded at successors of the start sentinel, propagated topologically, and merged when they are adjacent or overlap.
- **Graph mutation is post-alignment**: `POAGraph` is modified via `AddAlignment` only after the aligner returns an `AlignResult` — the graph itself is read-only during alignment.

### Active Development Context

The branch (`banded-poa`) implements the band-doubling optimization, and includes migration of the stable A* engine into the new code organization and trait hierarchy.

## Banded Alignment Algorithm Implementation

* First, compute bands for each node. Starting with an initial band on the start node, compute bands in topological order. Bands from all predecessor nodes are considered. We keep track per band what query range originates from which predecessor band, in `BandEdge` objects.
* Allocate a single chunk of memory for storing alignment scores, enough to hold scores for previously computed bands.
* Iterate over all bands, and compute alignment scores based on its predecessor bands.
* Check if we were able to reach the alignment end
