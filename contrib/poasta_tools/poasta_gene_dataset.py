#!/usr/bin/env python3
"""
Extract sequences for a specific gene name from RefSeq downloaded
genomes.

Use [genome_updater](https://github.com/pirovc/genome_updater) to download
RefSeq genomes, and use this script to extract specific gene sequences from them.
"""

import argparse
import gzip
import shlex
from xml.etree.ElementTree import XMLPullParser
from pathlib import Path
import subprocess
import re
import tqdm

import pandas
import skbio


def get_aliases(gene):
    proc = subprocess.Popen(
        "esearch -db gene -query {} | efetch -format docsum".format(shlex.quote(f"{gene}[Gene Name]")),
        shell=True, stdout=subprocess.PIPE,
    )

    all_aliases = set()
    xmlparser = XMLPullParser(events=('end'))
    for line in proc.stdout:
        xmlparser.feed(line)

        for event, elem in xmlparser.read_events():
            if event == 'end' and elem.tag == 'DocumentSummary':
                other = elem.find('OtherAliases').text
                all_aliases.update((v.strip() for v in other.split(',')))

                elem.clear()

    return all_aliases


gene_name_re = re.compile(r"\[gene=(.*?)\]")


def extract_genes(gene_name, refseq_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    aliases = get_aliases(gene_name)
    aliases.add(gene_name)
    print("Searching for gene", gene_name, f" (all aliases: {aliases})")

    with gzip.open(output_dir / f"{gene_name}_sequences.fna.gz", "wt") as o:
        for fname in tqdm(refseq_dir.glob("*_cds_from_genomic.fna.gz")):
            with gzip.open(fname, "rt") as f:
                for r in skbio.io.read(f, "fasta", validate=False):
                    if (match := gene_name_re.search(r.metadata['description'])):
                        print("Found gene name:", match)
                        entry_gene = match.group(1)

                        if entry_gene.strip() in aliases:
                            print("writing..")
                            skbio.io.write(r, "fasta", into=o)


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

    asm_summary_cols = [
        'asm_accession',
        'bioproject',
        'biosample',
        'wgs_master',
        'refseq_category',
        'taxid',
        'species_taxid',
        'organism_name',
        'infraspecific_name',
        'isolate',
        'version_status',
        'assembly_level',
        'release_type',
        'genome_rep',
        'seq_rel_date',
        'asm_name',
        'asm_submitter',
        'gbrs_paired_asm',
        'paired_asm_comp',
        'ftp_path',
        'excluded_from_refseq',
        'relation_to_type_material',
        'asm_not_live_date',
        'assembly_type',
        'group',
        'genome_size',
        'genome_size_ungapped',
        'gc_percent',
        'replicon_count',
        'scaffold_count',
        'contig_count',
        'annotation_provider',
        'annotation_name',
        'annotation_date',
        'total_gene_count',
        'protein_coding_gene_count',
        'non_coding_gene_count',
        'pubmed_id',
    ]

    df = pandas.read_csv(args.refseq / "assembly_summary.txt", sep='\t',
                         names=asm_summary_cols)

    for gene in args.genes:
        extract_genes(gene, args.refseq, args.output_dir / gene)


if __name__ == '__main__':
    main()



