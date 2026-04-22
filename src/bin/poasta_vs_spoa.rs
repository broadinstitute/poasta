use std::error::Error;
use std::fs::{create_dir_all, File};
use std::io::{BufReader, BufWriter};
use std::path::Path;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use noodles::fasta::record::Sequence;
use noodles::fasta::{self, Record};

use poasta::align::cost_models::affine::Affine;
use poasta::align::cost_models::linear::Linear;
use poasta::align::cost_models::two_piece::TwoPieceAffine;
use poasta::align::engine::band_doubling::BandDoublingEngineScalar;
use poasta::align::engine::AlignResult;
use poasta::align::traits::AlignmentEngine;
use poasta::cli::poasta::CostModelKind;
use poasta::cli::poasta_vs_spoa::PoastaVsSpoaArgs;
use poasta::graph::io::fasta::load_graph_from_fasta_msa;
use poasta::graph::poa::POAGraph;

fn main() -> Result<(), Box<dyn Error + 'static>> {
    let args = PoastaVsSpoaArgs::parse();

    let error_output_dir = args.output_dir.join("errors");
    create_dir_all(&error_output_dir)?;

    let is_gzipped = args
        .sequences
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    let reader_inner: Box<dyn std::io::BufRead> = if is_gzipped {
        Box::new(
            File::open(&args.sequences)
                .map(MultiGzDecoder::new)
                .map(BufReader::new)?,
        )
    } else {
        Box::new(File::open(&args.sequences).map(BufReader::new)?)
    };

    let mut reader = fasta::io::Reader::new(reader_inner);
    let mut records = reader.records();

    // Build the SPOA alignment engine matching the requested cost model
    let match_score: i8 = 0;
    let mismatch_score = -(args.cost_mismatch as i8);
    let gap_open_score = -(args.cost_gap_open as i8) - (args.cost_gap_extend as i8);
    let gap_extend_score = -(args.cost_gap_extend as i8);

    let mut spoa_engine = match args.cost_model {
        CostModelKind::Linear => spoa_rs::AlignmentEngine::new_linear(
            spoa_rs::AlignmentType::kNW,
            match_score,
            mismatch_score,
            gap_extend_score,
        ),
        CostModelKind::Affine => spoa_rs::AlignmentEngine::new_affine(
            spoa_rs::AlignmentType::kNW,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        ),
        CostModelKind::TwoPiece => spoa_rs::AlignmentEngine::new_convex(
            spoa_rs::AlignmentType::kNW,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
            -(args.cost_gap_open2 as i8) - (args.cost_gap_extend2 as i8),
            -(args.cost_gap_extend2 as i8),
        ),
    };

    // Use the first N sequences to build a SPOA graph
    eprintln!("Creating graph using SPOA...");
    let mut spoa_graph = spoa_rs::Graph::new();
    let mut seq_names = vec![];
    for _ in 0..args.num_graph {
        if let Some(r) = records.next() {
            let r = r?;
            let sequence = str::from_utf8(r.sequence().as_ref())?;

            let (_, aln) = spoa_engine.align(sequence, &spoa_graph);
            spoa_graph.add_alignment(aln, sequence);

            seq_names.push(r.definition().clone());
        } else {
            eprintln!("Not enough sequences in FASTA file to build the graph");
            return Ok(());
        }
    }

    // Write SPOA MSA to FASTA, used as seed graph for POASTA
    eprintln!("Writing graph...");
    let graph_fname = args.output_dir.join("graph.fasta");
    let mut writer = File::create(&graph_fname)
        .map(BufWriter::new)
        .map(fasta::io::Writer::new)?;

    let msa_aln = spoa_graph.generate_msa();
    assert_eq!(msa_aln.len(), args.num_graph);

    for (i, aln_seq) in msa_aln.into_iter().enumerate() {
        let seq = Sequence::from_iter(aln_seq.as_bytes().iter().copied());
        let record = Record::new(seq_names[i].clone(), seq);
        writer.write_record(&record)?;
    }
    drop(writer);

    // Import SPOA MSA into POASTA
    eprintln!("Reading POASTA graph...");
    let poasta_graph = load_graph_from_fasta_msa::<u32, _>(
        File::open(graph_fname).map(BufReader::new)?,
    )?;

    // Align remaining sequences with both engines and compare scores
    let (total_tested, total_incorrect) = match args.cost_model {
        CostModelKind::Affine => {
            let costs = Affine::new(
                args.cost_match,
                args.cost_mismatch,
                args.cost_gap_open,
                args.cost_gap_extend,
            );
            let engine = BandDoublingEngineScalar::<Affine, u32>::new(costs);
            run_comparison(&engine, &poasta_graph, &mut records, &mut spoa_engine, &spoa_graph, &error_output_dir)?
        }
        CostModelKind::Linear => {
            let costs = Linear::new(args.cost_match, args.cost_mismatch, args.cost_gap_extend);
            let engine = BandDoublingEngineScalar::<Linear, u32>::new(costs);
            run_comparison(&engine, &poasta_graph, &mut records, &mut spoa_engine, &spoa_graph, &error_output_dir)?
        }
        CostModelKind::TwoPiece => {
            let costs = TwoPieceAffine::new(
                args.cost_match,
                args.cost_mismatch,
                args.cost_gap_open,
                args.cost_gap_extend,
                args.cost_gap_open2,
                args.cost_gap_extend2,
            );
            let engine = BandDoublingEngineScalar::<TwoPieceAffine, u32>::new(costs);
            run_comparison(&engine, &poasta_graph, &mut records, &mut spoa_engine, &spoa_graph, &error_output_dir)?
        }
    };

    eprintln!("Total sequences tested: {total_tested}");
    eprintln!("Total incorrect:        {total_incorrect}");

    Ok(())
}

fn run_comparison<E>(
    engine: &E,
    poasta_graph: &POAGraph<u32>,
    records: &mut impl Iterator<Item = std::io::Result<fasta::Record>>,
    spoa_engine: &mut spoa_rs::AlignmentEngine,
    spoa_graph: &spoa_rs::Graph,
    error_output_dir: &Path,
) -> Result<(usize, usize), Box<dyn Error>>
where
    E: for<'s> AlignmentEngine<&'s [u8], Graph = POAGraph<u32>, Success = AlignResult<POAGraph<u32>>>,
    for<'s> <E as AlignmentEngine<&'s [u8]>>::Error: std::fmt::Display,
{
    let mut total_tested = 0;
    let mut total_incorrect = 0;

    while let Some(r) = records.next() {
        let r = r?;
        let seq_bytes: &[u8] = r.sequence().as_ref();

        eprintln!("Aligning {}..", str::from_utf8(r.name())?);

        let poasta_aln = engine
            .align(poasta_graph, seq_bytes)
            .map_err(|e| format!("poasta alignment failed: {e}"))?;

        let seq_str = str::from_utf8(seq_bytes)?;
        let (spoa_score, _) = spoa_engine.align(seq_str, spoa_graph);
        let spoa_penalty = (-spoa_score) as u32;

        if spoa_penalty != poasta_aln.score {
            eprintln!(
                "POASTA score {} != SPOA score: {}",
                poasta_aln.score,
                spoa_penalty,
            );

            let seq_name = String::from_utf8(r.name().to_owned())?;
            let mut output_file = File::create(error_output_dir.join(seq_name + ".fasta"))
                .map(BufWriter::new)
                .map(fasta::io::Writer::new)?;

            output_file.write_record(&r)?;
            total_incorrect += 1;
        }

        total_tested += 1;
    }

    Ok((total_tested, total_incorrect))
}
