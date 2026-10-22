# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed

- `POAGraph::in_edges` and `POAGraph::out_edges` are public (#37).

### Fixed

- The consensus follows each node's heaviest incoming edge, as in Lee's heaviest
  bundling. It used to maximise the summed edge weight, so a long insertion in a few
  sequences could outweigh the majority (#37). This can change consensus output.

## [1.0.0](https://github.com/broadinstitute/poasta/releases/tag/v1.0.0) - 2026-10-18

A rewrite of POASTA, twice as fast as the previous version.

### Added

- Linear, affine and two-piece affine gap costs.
- Gap-chaining seed heuristic (`--heuristic chain`).
- Semi-global alignment (`-m semi-global`, seed heuristic only).
- Consensus output (`--include-consensus`, `--consensus-only`).

## [0.1.0](https://github.com/broadinstitute/poasta/releases/tag/v0.1.0) - 2025-07-08

The first release of POASTA. POASTA is still very much in active development,
and APIs and CLIs are subject to change.

This version includes

- `poasta` CLI for computing partial order alignments, supporting gap-affine and two-piece cost models (thanks @ekg).
- `lasagne` CLI for computing sequence-to-DAG alignments, outputing alignments in GAF format. **(ALPHA)**
