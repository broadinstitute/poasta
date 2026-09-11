# Plan port of A* alignment engine to the new crate structure

The banded-poa branch of Poasta has many improved features,
such as better logging using the tracing crate, better trait
organization, and additional alignment engines.

For this refactor, we want to port the A* alignment algorithm
implemented on the `stable` branch to the new `banded-poa`
branch.

This enables more flexibility across settings: for example,
the banded-poa engine does not work for semi-global
alignment, while the A* engine could support that.

## Preliminaries

- Use modern, idiomatic Rust
- Format with `rustfmt`

## Review A* engine and plan interface design

The A* algorithm is quite different from a more typical
(banded) DP algorithm. Explore the current implementation
on `stable`: review the core A* loop, extension, and how
heuristics can be implemented. Review the bubble index
and how it's computed and stored. Heuristics might need
data from this index.

Design an updated trait and struct interface to implement
the A* engine in the current branch. It should support:

- Linear, affine, and 2-piece scoring
- Flexible heuristic configuration, that can alter behaviour
  depending on alignment mode (global vs semi-global).
- Bubble index construction
- Fit in with the current alignment engine traits

Write the interface to a file for me to review.

## Notes

During the design phase, consider performance. Aim to minimize
allocations, branches, and cache misses in hot loops.

Example: the A* engine allocates BxB blocks of the matrix at once.

