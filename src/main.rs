extern crate petgraph;
extern crate num;

use petgraph::prelude::*;
use petgraph::dot::Dot;

mod alignment;
mod graph;
mod wavefront;

use alignment::AlignedPair;
use graph::POAGraph;
use wavefront::aligner::WavefrontPOAligner;
use wavefront::compute::gap_affine::WFComputeGapAffine;


fn main() -> Result<(), String> {
    let mut poa_graph = POAGraph::new();

    let seq1 = b"AATGGTTGTCACGTCAGTAA";
    let weights1 = vec![1; seq1.len()];

    let seq2 = b"ATTGTAAAGTCTCGTCGGTTT";
    let weights2 = vec![1; seq2.len()];

    // EMBOSS_001         1 AATGGT--TGTCACGTCAGT     18
    //                       ||.||  .|||.||||.||
    // EMBOSS_001         1 -ATTGTAAAGTCTCGTCGGT     19

    let alignment = vec![
        AlignedPair{ rpos: Some(NodeIndex::from(0)), qpos: None},
        AlignedPair{ rpos: Some(NodeIndex::from(1)), qpos: Some(0)},
        AlignedPair{ rpos: Some(NodeIndex::from(2)), qpos: Some(1)},
        AlignedPair{ rpos: Some(NodeIndex::from(3)), qpos: Some(2)},
        AlignedPair{ rpos: Some(NodeIndex::from(4)), qpos: Some(3)},
        AlignedPair{ rpos: Some(NodeIndex::from(5)), qpos: Some(4)},
        AlignedPair{ rpos: None, qpos: Some(5)},
        AlignedPair{ rpos: None, qpos: Some(6)},
        AlignedPair{ rpos: Some(NodeIndex::from(6)), qpos: Some(7)},
        AlignedPair{ rpos: Some(NodeIndex::from(7)), qpos: Some(8)},
        AlignedPair{ rpos: Some(NodeIndex::from(8)), qpos: Some(9)},
        AlignedPair{ rpos: Some(NodeIndex::from(9)), qpos: Some(10)},
        AlignedPair{ rpos: Some(NodeIndex::from(10)), qpos: Some(11)},
        AlignedPair{ rpos: Some(NodeIndex::from(11)), qpos: Some(12)},
        AlignedPair{ rpos: Some(NodeIndex::from(12)), qpos: Some(13)},
        AlignedPair{ rpos: Some(NodeIndex::from(13)), qpos: Some(14)},
        AlignedPair{ rpos: Some(NodeIndex::from(14)), qpos: Some(15)},
        AlignedPair{ rpos: Some(NodeIndex::from(15)), qpos: Some(16)},
        AlignedPair{ rpos: Some(NodeIndex::from(16)), qpos: Some(17)},
        AlignedPair{ rpos: Some(NodeIndex::from(17)), qpos: Some(18)},
        AlignedPair{ rpos: Some(NodeIndex::from(18)), qpos: None},
        AlignedPair{ rpos: Some(NodeIndex::from(19)), qpos: None},
    ];

    poa_graph.add_alignment_with_weights(seq1, None, &weights1)?;
    poa_graph.add_alignment_with_weights(seq2, Some(&alignment), &weights2)?;

    let transformed = poa_graph.graph.map(
        |ix, data|
            format!("{:?} ({:?})", char::from(data.symbol), poa_graph.get_node_rank(ix)),
        |_, data|
            format!("{}, {:?}", data.weight, data.sequence_ids)
    );

    let dot = Dot::new(&transformed);
    println!("{:?}", dot);

    //        let seq1 = b"AATGGTTGTCACGT------CAGT";
    let seq3 = b"TTGTCAACGTCAGTAAAAA";

    let mut aligner: WavefrontPOAligner<u32, WFComputeGapAffine<u32>> = WavefrontPOAligner::new(&poa_graph, "output");
    let alignment = aligner.align(seq3);

    let (ref_aln, qry_aln): (String, String) = alignment
        .into_iter()
        .map(|pair|
            (pair.rpos.map(|nix| char::from(poa_graph.graph[nix].symbol)).unwrap_or('-'),
             pair.qpos.map(|qpos| char::from(seq3[qpos])).unwrap_or('-'))
        )
        .unzip();

    eprintln!("  Graph: {:?}", ref_aln);
    eprintln!("  Query: {:?}", qry_aln);

    Ok(())
}
