use std::error::Error;
use std::fs::{create_dir_all, File};
use std::io::{BufReader, BufWriter, Write};

use clap::Parser;
use flate2::read::MultiGzDecoder;
use noodles::fasta::record::Sequence;
use noodles::fasta::{self, Record};

use poasta::align::astar::heuristic::{self, MinGapCost};
use poasta::align::cost_models::affine::Affine;
use poasta::align::PoastaAligner;
use poasta::cli::poasta_vs_spoa::PoastaVsSpoaArgs;
use poasta::graph::poa::POASeqGraph;

fn main() -> Result<(), Box<dyn Error + 'static>> {
    let args = PoastaVsSpoaArgs::parse();

    // Make directory to write sequences with faulty alignments
    let error_output_dir = args.output_dir.join("errors");
    create_dir_all(&error_output_dir)?;

    let is_gzipped = args
        .sequences
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    // Check if we have a gzipped file
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

    // Use the first N sequences to construct a graph
    let mut graph = spoa_rs::Graph::new();
    let mismatch = args.cost_mismatch.unwrap_or(4) as i8;
    let gap_open = args.cost_gap_open.unwrap_or(6) as i8;
    let gap_extend = args.cost_gap_extend.unwrap_or(2) as i8;
    let mut engine = spoa_rs::AlignmentEngine::new_affine(
        spoa_rs::AlignmentType::kNW,
        0,
        -mismatch,
        -gap_open - gap_extend,
        -gap_extend,
    );

    eprintln!("Creating graph using SPOA...");
    let mut seq_names = vec![];
    for _ in 0..args.num_graph {
        if let Some(r) = records.next() {
            let r = r?;
            let sequence = str::from_utf8(r.sequence().as_ref())?;

            let (_, aln) = engine.align(sequence, &graph);
            graph.add_alignment(aln, sequence);

            seq_names.push(r.definition().clone());
        } else {
            eprintln!("Not enough sequences in FASTA graph");

            // TODO: Better error return
            return Ok(());
        }
    }

    // Write graph MSA to FASTA, to be used as input for POASTA
    eprintln!("Writing graph...");
    let graph_fname = args.output_dir.join("graph.fasta");
    let mut writer = File::create(&graph_fname)
        .map(BufWriter::new)
        .map(fasta::io::Writer::new)?;

    let msa_aln = graph.generate_msa();

    assert_eq!(msa_aln.len(), args.num_graph);

    for (i, aln_seq) in msa_aln.into_iter().enumerate() {
        let seq = Sequence::from_iter(aln_seq.as_bytes().iter().copied());
        let record = Record::new(seq_names[i].clone(), seq);

        writer.write_record(&record)?;
    }

    drop(writer);

    // Import MSA as POASTA graph
    eprintln!("Reading POASTA graph...");
    let mut reader = File::open(graph_fname).map(BufReader::new)?;

    let poasta_graph = POASeqGraph::<u32>::try_from_fasta_msa(&mut reader)?;

    let scoring = Affine::<i32, u32>::new(mismatch as u8, gap_open as u8, gap_extend as u8);

    let heuristic = MinGapCost::new(scoring);
    let aligner = PoastaAligner::new(heuristic);

    // Align remaining sequences and compare alignment scores
    let mut total_tested = 0;
    let mut total_incorrect = 0;
    while let Some(r) = records.next() {
        let r = r?;

        eprintln!("Aligning {}..", str::from_utf8(r.name())?);

        let poasta_aln = aligner.align(
            &poasta_graph,
            r.sequence().as_ref(),
            poasta::align::AlignmentMode::Global,
        )?;

        let seq = str::from_utf8(r.sequence().as_ref())?;
        let (spoa_score, _) = engine.align(seq, &graph);

        if -spoa_score != poasta_aln.score.as_usize() as i32 {
            eprintln!(
                "POASTA score {} != SPOA score: {}",
                poasta_aln.score.as_usize(),
                -spoa_score
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

    eprintln!("Total sequences tested: {total_tested}");
    eprintln!("Total incorrect:        {total_incorrect}");

    Ok(())
}
