//! Unified FASTA/FASTQ sequence reading with transparent gzip support.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use noodles::{fasta, fastq};

use crate::errors::PoastaIOError;

/// Sequence file format, inferred from the filename extension.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

impl SeqFormat {
    /// Detect the format from a path's extension. Strips a trailing `.gz` first.
    pub fn from_path(path: &Path) -> Result<Self, PoastaIOError> {
        let name = path
            .file_name()
            .ok_or(PoastaIOError::InvalidFormat)?
            .to_string_lossy()
            .to_lowercase();

        let stripped = name.strip_suffix(".gz").unwrap_or(&name);

        if stripped.ends_with(".fa")
            || stripped.ends_with(".fasta")
            || stripped.ends_with(".fna")
            || stripped.ends_with(".msa")
        {
            Ok(SeqFormat::Fasta)
        } else if stripped.ends_with(".fq") || stripped.ends_with(".fastq") {
            Ok(SeqFormat::Fastq)
        } else {
            Err(PoastaIOError::InvalidFormat)
        }
    }
}

/// Open `path` (gzipped or plain) and return a `BufRead` handle.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead>, PoastaIOError> {
    let file = File::open(path).map_err(|source| PoastaIOError::FileReadError { source })?;
    let is_gzipped = path
        .file_name()
        .map(|v| v.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false);

    if is_gzipped {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Iterate sequences from a FASTA or FASTQ file (plain or gzipped).
///
/// Yields `(name, sequence)` pairs. Format is auto-detected from the filename.
pub fn open_sequences(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Result<(String, Vec<u8>), PoastaIOError>>>, PoastaIOError> {
    let format = SeqFormat::from_path(path)?;
    let reader = open_reader(path)?;

    tracing::debug!(path = %path.display(), ?format, "opening sequence file");

    match format {
        SeqFormat::Fasta => {
            let mut fasta_reader = fasta::io::Reader::new(reader);
            let iter = std::iter::from_fn(move || {
                let mut def_buf = String::new();
                match fasta_reader.read_definition(&mut def_buf) {
                    Ok(0) => return None,
                    Ok(_) => {}
                    Err(e) => return Some(Err(PoastaIOError::FileReadError { source: e })),
                }
                let name = match def_buf
                    .trim_start_matches('>')
                    .split_ascii_whitespace()
                    .next()
                {
                    Some(n) => n.to_owned(),
                    None => return Some(Err(PoastaIOError::InvalidFormat)),
                };

                let mut seq_buf: Vec<u8> = Vec::new();
                match fasta_reader.read_sequence(&mut seq_buf) {
                    Ok(_) => Some(Ok((name, seq_buf))),
                    Err(e) => Some(Err(PoastaIOError::FileReadError { source: e })),
                }
            });
            Ok(Box::new(iter))
        }
        SeqFormat::Fastq => {
            let mut fastq_reader = fastq::io::Reader::new(reader);
            let iter = std::iter::from_fn(move || {
                let mut record = fastq::Record::default();
                match fastq_reader.read_record(&mut record) {
                    Ok(0) => None,
                    Ok(_) => match std::str::from_utf8(record.name()) {
                        Ok(n) => Some(Ok((n.to_owned(), record.sequence().to_vec()))),
                        Err(e) => Some(Err(PoastaIOError::InvalidUtf8 { source: e })),
                    },
                    Err(e) => Some(Err(PoastaIOError::FileReadError { source: e })),
                }
            });
            Ok(Box::new(iter))
        }
    }
}
