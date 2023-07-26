#!/usr/bin/env python3
"""
Extract sequences for a specific gene name from RefSeq downloaded
genomes.

Use [genome_updater](https://github.com/pirovc/genome_updater) to download
RefSeq genomes, and use this script to extract specific gene sequences from them.
"""

import argparse
import gzip
import re
from pathlib import Path
from tqdm import tqdm

import skbio

gene_name_re = re.compile(r"\[gene=(.*?)\]")


def extract_genes(gene_name, refseq_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Searching for gene", gene_name)

    num_files_to_search = len(list(refseq_dir.glob("*_cds_from_genomic.fna.gz")))
    with (gzip.open(output_dir / f"{gene_name}_sequences.fna.gz", "wt") as o,
          open(output_dir / f"{gene_name}_accession_map.tsv", "w") as acc_o):
        for fname in tqdm(refseq_dir.glob("*_cds_from_genomic.fna.gz"), total=num_files_to_search):
            parts = fname.name.split('_', maxsplit=2)
            accession = "_".join(parts[:2])
            asm_name = parts[2].replace("_cds_from_genomic.fna.gz", "")

            with gzip.open(fname, "rt") as f:
                for r in skbio.io.read(f, "fasta"):
                    if (match := gene_name_re.search(r.metadata['description'])):
                        entry_gene = match.group(1)

                        if entry_gene.strip() == gene_name:
                            r.metadata['description'] += f" [accession={accession}] [asm_name={asm_name}]"
                            skbio.io.write(r, "fasta", into=o)
                            print(r.metadata['id'], accession, sep='\t', file=acc_o)


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'refseq', type=Path,
        help="Path to directory with RefSeq genomes (as downloaded by genome_updater.sh)."
    )

    parser.add_argument(
        'genes', nargs='+', metavar='GENE',
        help='The gene name(s) for which to extract the DNA sequence(s)'
    )

    parser.add_argument(
        '-o', '--output-dir', type=Path,
        help="Output directory. A subdirectory per gene will be created."
    )

    args = parser.parse_args()

    for gene in args.genes:
        extract_genes(gene, args.refseq, args.output_dir / gene)


if __name__ == '__main__':
    main()



