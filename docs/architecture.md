# POASTA architecture

POASTA aligns a sequence to a partial order alignment (POA) graph with an A*
search over the alignment graph, then folds the result back into the POA
graph. This document explains how the crate is laid out, how data flows
through it, and how the A* engine works, in the order the code runs. It is
written for someone who knows sequence alignment and Rust but has not read
the code.

## Crate layout

```
src/
├── lib.rs                  crate root: align, graph, errors, cli, debug
├── errors.rs               GraphError, PoastaIOError, PoastaError
│
├── align/                  aligning a query to a graph
│   ├── mod.rs              AlignmentMode, AlignState, PoastaAligner, RunStats
│   ├── traits.rs           AlignableGraph, AlignmentEngine, AlignOutput, AlignedPair
│   ├── score.rs            u32 scores with INF as "unreached"
│   ├── offsets.rs          OffsetType: the query offset width (u32)
│   ├── cost_models/        AlignmentCostModel and its three models
│   │   ├── linear.rs       mismatch, gap_extend
│   │   ├── affine.rs       mismatch, gap_open, gap_extend
│   │   └── two_piece.rs    two affine pieces, cheaper one wins
│   └── astar/              the A* engine
│       ├── mod.rs          AstarEngine, AstarCostModel, HeuristicConfig, main loop
│       ├── aln_graph.rs    alignment-graph state machines (1, 3 or 5 states)
│       ├── visited.rs      BlockedVisitedStorage, AstarData, backtrace
│       ├── queue.rs        LayeredQueue bucket queue, StateQueueLayer
│       ├── dfa.rs          DepthFirstGreedyAlignment: free match extension
│       ├── heuristic.rs    Dijkstra, MinimumGapCost
│       ├── pruning.rs      ReachedBubbleExitsMatch: bubble-exit pruning
│       └── dp_oracle.rs    test-only full-matrix DP used as the score oracle
│
├── graph/                  the POA graph and everything that works on it
│   ├── mod.rs
│   ├── traits.rs           GraphBase, GraphWithNodeOrdering, GraphNodeId
│   ├── poa.rs              POAGraph on a petgraph StableDiGraph
│   ├── alignment.rs        AddAlignment: fold an alignment into a graph
│   ├── consensus.rs        heaviest-bundle consensus
│   ├── utils.rs            reverse post-order helper
│   ├── bubbles/            superbubble detection and index
│   │   ├── finder.rs       SuperbubbleFinder (Gärtner et al. 2018)
│   │   └── index.rs        BubbleIndex: node to bubble map, distances to end
│   ├── io/                 FASTA/FASTQ, FASTA MSA, GFA, GAF, DOT
│   └── mock.rs             test-only graph builders
│
├── cli/                    clap argument structs (feature "cli")
│   ├── mod.rs              AlignmentSpan, LogLevel, CostModelKind
│   ├── poasta.rs           poasta align / view / stats
│   ├── lasagna.rs          lasagna
│   └── poasta_vs_spoa.rs   benchmark harness args (feature "bench")
│
├── debug/                  --debug-output-dir support (feature "cli")
│   ├── mod.rs              write_graph_dot, sanitise_filename
│   └── subscriber.rs       AlignStateLayer: A* trace events to TSV
│
└── bin/                    thin drivers: poasta.rs, lasagna.rs, poasta_vs_spoa.rs
```

Dependencies point one way: `graph` knows nothing about aligning, `align`
builds on `graph`, and `cli`, `debug` and the binaries build on both. The
one exception is deliberate: `graph/poa.rs` implements `AddAlignment` for
the engine's `AlignOutput`, so it imports `AlignOutput` and `AlignedPair`
from `align::traits`. Nothing else in `graph/` may import from `align/`.

## Data flow

`poasta align` builds or extends a graph from a FASTA file:

1. `bin/poasta.rs` parses arguments, installs the `tracing` subscriber,
   reads the sequences and an optional seed graph (a FASTA MSA), and picks
   a cost model and heuristic.
2. It constructs an `AstarEngine<POAGraph<u32>, C, H>` and hands it, with
   the graph, to `PoastaAligner`.
3. `PoastaAligner::align_one` aligns one sequence with `AlignmentEngine::align`
   and then calls `AddAlignment::add_alignment` on the graph. The first
   sequence skips alignment and seeds the graph.
4. The graph is written as a FASTA MSA or GFA, optionally with a consensus.

`lasagna` aligns reads to a fixed GFA graph. It builds the graph and its
`BubbleIndex` once, shares one engine across worker threads with
`with_existing_bubbles`, and emits one GAF record per read. The graph is
never modified.

The engine only ever reads the graph. Mutation happens afterwards, through
`AddAlignment`, and only in the POA-construction flow.

## The graph

`POAGraph<Ix>` wraps a petgraph `StableDiGraph`. Every node carries one
symbol byte; every edge carries a weight and the IDs of the sequences that
traverse it. Two sentinel nodes, start (`#`) and end (`$`), bracket the real
nodes so every path has a single source and sink. `post_process` assigns a
topological rank to each node and the shortest and longest depth from the
start, both of which the DOT writer shows.

Three traits describe what the engine needs from a graph:

- `GraphBase` (`graph/traits.rs`): node iteration, successors, predecessors,
  degrees.
- `GraphWithNodeOrdering`: start and end node, `node_rank` and
  `rank_to_node`. Ranks are a topological order, so a predecessor always has
  a lower rank than its successor. The visited storage, the bubble index and
  consensus all rely on this.
- `AlignableGraph` (`align/traits.rs`): the node symbol, and
  `is_symbol_equal`, which the greedy match extension calls in its hottest
  loop.

`consensus.rs` implements Lee's heaviest bundling: a forward pass in rank
order that maximises accumulated edge weight, then a walk back from the end
sentinel. `io/` reads FASTA and FASTQ (gzip aware), reads and writes FASTA
MSAs and GFA 1.0, writes GAF for lasagna and DOT for debugging.

## The A* engine

Alignment is a shortest-path problem over the alignment graph: the product
of the POA graph with the query positions, with one copy per alignment
state. A node in that graph is an `AlignmentGraphNode` (graph node plus
query offset) together with an `AlignState`. Edge costs come from the cost
model. The end state is the end sentinel at query offset `m`.

### Cost models and `AstarCostModel`

`AlignmentCostModel` (`cost_models/mod.rs`) exposes mismatch, gap open and
gap extend costs for up to two gap pieces, plus `gap_cost(state, length)`,
which the heuristic and the pruning check use. Single-piece models report
zero for piece two.

`AstarCostModel` (`astar/mod.rs`) attaches the engine machinery shaped for
each model: how many visited slots a cell needs (`STATES`), which alignment
graph expands its states, which visited storage holds its scores and which
queue layer buckets them. One `AstarEngine` serves all three models through
this trait.

| model | states | alignment graph | slots per cell |
|---|---|---|---|
| linear | Match | `LinearAlignmentGraph` | 1 |
| affine | Match, Deletion, Insertion | `GapAlignmentGraph<Affine>` | 3 |
| two-piece | plus Deletion2, Insertion2 | `GapAlignmentGraph<TwoPieceAffine>` | 5 |

Slots are indexed by `AlignState as usize`, in the order Match, Deletion,
Insertion, Deletion2, Insertion2.

### Alignment graphs

`AlignmentGraph` (`aln_graph.rs`) provides `initial_states` and `is_end`
from the alignment mode, and requires four expansion methods. The main loop
never expands a `Match` state directly. It runs the match DFA from it and
calls the expansion for the reason the DFA stopped:

- `expand_mismatch`: the next graph node does not match the next query base.
  Queue the mismatch and open a gap in each direction.
- `expand_ref_graph_end`: the graph ran out before the query. Only an
  insertion is possible.
- `expand_query_end`: the query ran out before the graph. Only a deletion is
  possible.
- `expand_all`: a gap state popped from the queue. Close it for free into
  `Match` and extend it by one.

The linear graph has a single state. With no gap-open cost there is nothing
to remember between steps, so a gap is an edge that costs `gap_extend` and
moves along one axis, and `expand_all` is unreachable.

`GapAlignmentGraph` handles one or two gap pieces. Each piece opens directly
from `Match` under its own open plus extend cost and only ever extends
within itself. The two pieces are independent tracks with no transitions
between them, which is what the `min(M + go2 + ge2, D2 + ge2)` recurrence
means.

Ends-free alignment (semi-global) is shared by all graphs: initial states
are every graph node within the free leading-graph budget crossed with every
query offset within the free leading-query budget, all at score zero, and a
`Match` state is an end once its offset is within the free trailing-query
budget and its node can reach the end within the free trailing-graph budget.

### Visited storage

`BlockedVisitedStorage<N, O, S, B>` (`visited.rs`) stores one `[u32; S]` per
`(node rank, query offset)` cell, in `B × B` blocks (`B = 8`) allocated on
demand in a hash map per block row. Blocking keeps the cells a search
touches next close in memory. Unreached cells read as `INF`
(`score.rs`), and `update_score_if_lower` is a plain comparison against it.

`AstarData<N, O, C, S>` wraps that storage with the cost model, the bubble
index and, per bubble exit, the sorted query offsets at which the exit was
reached in `Match`. It implements `AstarVisited`, the interface the main
loop, the DFA and the heuristic see.

### The queue

`LayeredQueue` (`queue.rs`) is a bucket queue: a deque of layers indexed by
priority `score + h`. Push and pop are constant time, and the front layer is
always the minimum. `StateQueueLayer<N, O, S>` is one bucket, with a vector
per alignment state so the main loop's state dispatch stays predictable.
Within a bucket, pops are LIFO and ordered Match, Deletion, Deletion2,
Insertion, Insertion2, which keeps the search depth-first among equals and
pairs with the greedy extension below.

### Main loop and greedy match extension

`run_astar_search` (`astar/mod.rs`) is the loop. It seeds the queue with the
initial states, then repeatedly pops the best state, skips it if a better
score has since been recorded, stops if it is an end state, drops it if
pruning says it cannot improve on a reached bubble exit, and otherwise marks
it reached and expands it.

Match edges cost nothing, so a popped `Match` state is not expanded one step
but extended as far as it will go. `DepthFirstGreedyAlignment` (`dfa.rs`)
walks matching successors depth first, recording each as visited at the
same score, and returns only when it hits a mismatch, the graph end, the
query end or an ends-free end state. This is the single largest
constant-factor win in the engine.

### Heuristics

`AstarHeuristic::h` returns a lower bound on the cost still to pay.
`Dijkstra` returns zero. `MinimumGapCost` uses the bubble index's shortest
and longest distance from the node to the graph end: if even the shortest
remaining graph path is longer than the remaining query, deletions are
unavoidable, and if even the longest is shorter, insertions are. The bound
is the cost model's `gap_cost` for that length, without a reopen charge
when the state is already in a gap of that kind. Under ends-free modes a
gap that fits inside a free region costs nothing, so the bound is zero
there, which keeps it admissible.

`HeuristicConfig` (`DijkstraHeuristic`, `MinGapCostHeuristic`) is the engine's
`H` type parameter and builds the heuristic per alignment.

### Bubble index and pruning

`SuperbubbleFinder` (`bubbles/finder.rs`) finds superbubbles in the DAG.
`BubbleIndex` (`bubbles/index.rs`) records, per node, which bubbles it lies
in and its distance to each exit, plus the shortest and longest distance to
the graph end. It is built once per alignment, or once per graph for
lasagna.

Every path through a bubble reaches the bubble's exit. Once the exit has
been reached in `Match` at some offset and score, a state inside the bubble
can only matter if it can reach the exit at an offset and score no reached
entry already dominates. `ReachedBubbleExitsMatch` (`pruning.rs`) makes
that check using the minimum gap cost between the state and the reached
offsets. States that fail it are dropped from the queue and from DFA
extension. Pruning changes the number of states visited but never the
score; the tests assert both.

### Backtrace

Scores alone determine the alignment. `AstarData::backtrace` walks back
from the end state, at each step finding the predecessor whose recorded
score is consistent with the edge into the current state. For the gap
models the walk carries the state along and treats an equal score in a gap
slot at the same cell as the free gap-closing edge. For the linear model
each step returns which kind of move led into the current cell. The result
is a `Vec<AlignedPair>`, oldest first.

### Statistics

`AlignmentStats::cells_computed` is the number of states expanded or
reached by extension. `fraction_of_full_matrix` divides that by
`(m + 1) × real nodes × states`, the size of the full dynamic programming
matrix, so a value of one means as much work as a full DP.

## Alignment modes

`AlignmentMode::Global` aligns the whole query to a whole path from start
to end. `AlignmentMode::EndsFree` carries four `Bound<usize>` budgets: free
query prefix and suffix, and free graph prefix and suffix, each measured in
bases or hops. The CLI's `--alignment-span semi-global` maps to unbounded
free graph ends and no free query ends.

## Tracing and debug output

The engine logs through the `tracing` crate. At info level the binaries
report configuration and run statistics. Debug level adds per-sequence
scores. Trace level carries the search itself:

| target | fired when | fields |
|---|---|---|
| `poasta::astar::set_visited` | a state is popped and expanded, or accepted as the end | score, node, offset, state |
| `poasta::astar::extend` | the DFA reaches a state by a free match | score, node, offset, extended_offset, state |
| `poasta::astar::queue_item` | a state is pushed onto the queue | score, node, offset, state, priority |

`--debug-output-dir` installs `AlignStateLayer` (`debug/subscriber.rs`),
which opens `<seq>_visited.tsv` and `<seq>_queued.tsv` whenever a span with
a `seq_id` field starts (the `align_seq` span in `PoastaAligner` and in
lasagna's workers) and appends every event fired inside it. The binaries
also write `graph_for_<seq>.dot`, the graph the sequence was aligned
against. The visited row count equals `cells_computed`.

The level filter is attached to the log formatter only, so trace events
reach the file layer without being printed. Trace calls cost a cached
branch when no subscriber wants them; the `disable_tracing` cargo feature
compiles them out for benchmark builds, in which case the debug directory
receives only DOT files.

`python/poasta-tools` predates this instrumentation and reads a different
file format. It needs updating before it can plot these files.

## Command-line tools

- `poasta align`: build or extend an MSA. Cost model, costs, heuristic,
  alignment span, output format, consensus and debug output are flags.
  `view` and `stats` are declared but not implemented on this branch.
- `lasagna`: align reads to a GFA graph and write GAF, with `-j` worker
  threads.
- `poasta_vs_spoa`: compare scores against SPOA on a FASTA file. Needs
  `--features bench`, which pulls in the `spoa_rs` git dependency and its
  C++ build.

## Testing

Everything is unit tested in-tree; `tests/` holds only fixture files.

- `DpOracle` (`astar/dp_oracle.rs`) is a plain full-matrix Gotoh DP for all
  three cost models. Every A* configuration must reproduce its score on the
  test graphs, with either heuristic.
- `test_support::score_from_alignment` replays an alignment through the
  cost model, so a backtrace is checked against the score the engine
  reported, independently of how the scores were computed.
- Ends-free cases are hand-computed, since the oracle is global only.
- Pruning is checked to leave scores unchanged and never visit more states.
- `graph/mock.rs` builds small topologies for the bubble and storage tests
  without symbols or alignment.

## Invariants worth knowing

- The graph is read-only during alignment.
- Query offsets are `u32` (`Offset` in `astar/mod.rs`); queries above 4 Gbp
  would need a wider `OffsetType`.
- `INF` is `1 << 28`, so adding any realistic cost to it cannot overflow.
- Visited slots and queue vectors are indexed by `AlignState as usize`.
- Within a queue bucket, order is LIFO and Match before deletions before
  insertions. Changing it changes which of several equal-score alignments
  is returned, and therefore the MSA layout.
