<h1 align="center">🍝 POASTA</h1>
<h2 align="center">Fast, gap-affine sequence-to-graph and partial order aligner</h2>

POASTA is a fast and optimal partial order aligner that support gap-affine alignment penalties. Inspired by 
a recent [algorithm for pairwise alignment](https://github.com/smarco/WFA2-lib), it can exploit exact matches
between the query and the graph, greatly speeding up the alignment process.

## Installation

### Pre-built binaries

TODO

### Conda

TODO

### Building POASTA from source

#### Rust compiler

POASTA is written in Rust, and to build and install it, you'll need a recent version of the Rust compiler. The 
minimum supported Rust version is 1.70.

1. Download and install `rustup`: https://rustup.rs/
2. Run `rustup update`

#### Building POASTA

1. Clone the repository. 
    
   ```bash
   git clone https://github.com/broadinstitute/poasta
   ```
 
2. Move into the directory. 

   ```bash
   cd poasta
   ```
 
3. Build using `cargo`. We enable a flag to ensure the compiler uses all features of your machine's CPU. 
   To maximize portability of the binary, however, remove the `RUSTFLAGS="..."` part.
    
   ```bash
   RUSTFLAGS="-C target-cpu=native" cargo build --release
   ```
 
4. The built `poasta` executable will be available in the directory `target/release/`


## Usage

### Creating an alignment from scratch

To create a multiple sequence alignment from scratch, simply give it a FASTA. The FASTA file can be compressed 
with gzip (filename should have a `.gz` extension).

```bash
poasta align -o graph.poasta sequences.fasta
```

This will output the graph to a binary file called `graph.poasta`. POASTA can reuse this file to later align
additional sequences to it.

### Re-using an earlier alignment

To align additional sequences to an earlier created partial order graph, specify the existing graph using the 
`-g` option.

```bash
poasta align -g graph.poasta -o graph_updated.poasta new_sequences.fasta
```

This will import the graph stored in `graph.poasta`, then align the additional sequences in `new_sequences.fasta` to
this graph, and outputs the updated graph to `graph_updated.poasta`.

### Importing an existing multiple sequence alignment in FASTA format

POASTA can import an existing multiple sequence alignment stored in columnar FASTA format (e.g., those 
created by other tools like `mafft` or `spoa`), create the equivalent partial order graph from the existing alignment,
and then align new sequences to it. To achieve this, specify the FASTA MSA with extension .fa, .fna, or .fasta with
the `-g` option (file is also allowed to be compressed with gzip if it has a `.gz` suffix).

```bash
poasta align -g msa.fasta -o graph_updated.poasta new_sequences.fasta
```

### Other output formats

The default output format is POASTA's binary file format storing the graph data structure. With the `-O` option you
can specify a different output format.

Supported formats:

* DOT (GraphViz): Specify with `-O dot`
* TODO: ~~FASTA MSA~~
* TODO: ~~Graph GFA~~

For example, to visualize the graph with GraphViz:

```bash
poasta align -Odot sequences.fasta | dot -Tpng -o graph.png
```

Note that we did not specify an output file for `poasta align` (we did not use the `-o` option). If no output filename
is given, standard output will be used, so the output can be directly piped to `dot` to create the visualization.






