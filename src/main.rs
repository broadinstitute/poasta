extern crate petgraph;
extern crate num;

mod graph;
mod wavefront;

use std::convert::identity;

use graph::{Alphabet, ASCIIAlphabet, POAGraph, AlignedPair};

use petgraph::prelude::*;
use crate::wavefront::aligner::WavefrontPOAligner;
use crate::wavefront::compute::gap_affine::WFComputeGapAffine;


fn seq_to_alphabet_code<T, C, A: Alphabet<T, C>>(alphabet: &A, seq: &[T]) -> Result<Vec<C>, String>
where
    C: Into<char>
{
    let (coded_sequence_opt, failed): (Vec<_>, Vec<_>) = seq.iter()
        .map(|e| alphabet.encode(e))
        .partition(Option::is_some);

    if !failed.is_empty() {
        let failed_chars: Vec<(usize, char)> = failed.into_iter().filter_map(identity).map(|t| t.into()).enumerate().collect();
        return Err(format!("Found symbols not in the alphabet: {:?}", failed_chars));
    }

    Ok(coded_sequence_opt.into_iter().filter_map(identity).collect())
}


fn main() -> Result<(), String> {
    let alphabet = ASCIIAlphabet::new(b"ACGT");
    let mut poa_graph = POAGraph::new();

    let seq1 = b"AATGGTTGTCACGTCAGT";
    let seq1_coded = seq_to_alphabet_code(&alphabet, seq1)?;
    let weights1 = vec![1; seq1.len()];

    let seq2 = b"ATTGTAAAGTCTCGTCGGT";
    let seq2_coded = seq_to_alphabet_code(&alphabet, seq2)?;
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
        AlignedPair{ rpos: Some(NodeIndex::from(17)), qpos: Some(18)}
    ];

    poa_graph.add_alignment_with_weights(&seq1_coded, None, &weights1)?;
    poa_graph.add_alignment_with_weights(&seq2_coded, Some(&alignment), &weights2)?;

    let seq3 = b"AATGGTTGTCACGTCAGT";
    let seq3_coded = seq_to_alphabet_code(&alphabet, seq3)?;

    println!("Ranked nodes: {:?}", poa_graph.rank_to_node);
    println!("Node labels: {:?}", poa_graph.graph.node_indices().map(|v| poa_graph.graph[v].code).collect::<Vec<u8>>());
    println!("Seq 3 coded: {:?}", seq3_coded);

    let mut aligner: WavefrontPOAligner<u32, WFComputeGapAffine<u32>, u8> = WavefrontPOAligner::new(&poa_graph);
    aligner.align(&seq3_coded);

    // let transformed = poa_graph.graph.map(
    //     |_, data| char::from(alphabet.decode(data.code).unwrap()),
    //     |_, data| format!("{}, {:?}", data.weight, data.sequence_ids));
    //
    // let dot = Dot::new(&transformed);
    //
    // println!("{:?}", dot);

    Ok(())
}
