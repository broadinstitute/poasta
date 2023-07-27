#!/usr/bin/env python3
"""
Cluster gene sequences using Mash sketches, and create smalller subsets of the
gene sequences for benchmarking.
"""

import argparse
import gzip
from pathlib import Path

import numpy
import pandas
import skbio
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster


def phylip_to_matrix(fname):
    with open(fname) as f:
        fiter = iter(f)
        num_entries = int(next(fiter))
        reference_names = []

        matrix = numpy.zeros((num_entries, num_entries))
        for i, line in enumerate(fiter):
            first_tab = line.find('\t')
            ref_name = line[:first_tab].strip()
            reference_names.append(ref_name)

            if first_tab >= 0:
                row = numpy.fromstring(line[first_tab+1:], sep='\t')
                matrix[i, 0:len(row)] = row

        index_upper = numpy.triu_indices(num_entries)
        matrix[index_upper] = matrix.T[index_upper]

    dmatrix = pandas.DataFrame(matrix, index=reference_names,
                               columns=reference_names)

    print(dmatrix)
    return dmatrix


def create_dataset(asm_summary, acc_map, output_dir, dmatrix, sequences_fasta, max_dist, graph_ix, align_ix):
    output_dir.mkdir(exist_ok=True, parents=True)

    graph_ids = set(dmatrix.index[graph_ix])
    align_ids = set(dmatrix.index[align_ix])

    accessions_graph = []
    accessions_align = []
    lengths_graph = []
    lengths_align = []

    with (gzip.open(sequences_fasta, "rt") as ifasta,
          gzip.open(output_dir / "graph_seq.fna.gz", "wt") as graph_out,
          gzip.open(output_dir / "align_seq.fna.gz", "wt") as align_out):

        for r in skbio.io.read(ifasta, "fasta"):
            accession = acc_map.loc[r.metadata['id'], 'accession']

            if r.metadata['id'] in graph_ids:
                skbio.io.write(r, "fasta", into=graph_out)
                accessions_graph.append(accession)
                lengths_graph.append(len(r))
            elif r.metadata['id'] in align_ids:
                skbio.io.write(r, "fasta", into=align_out)
                accessions_align.append(accession)
                lengths_align.append(len(r))

    with open(output_dir / "meta.toml", "w") as meta_out:
        print(f"clustering_max_dist = {max_dist:g}", file=meta_out)
        print(file=meta_out)

        print("[graph_set]", file=meta_out)
        print('fname = "graph_seq.fna.gz"', file=meta_out)
        print(f"num_seqs = {len(graph_ix)}", file=meta_out)
        print(f"avg_seq_len = {sum(lengths_graph)/len(graph_ix):.1f}", file=meta_out)
        avg_dist = dmatrix.iloc[graph_ix, graph_ix].mean(axis=None)
        print(f"avg_pairwise_dist = {avg_dist:g}", file=meta_out)
        print(file=meta_out)
        print("[graph_set.species]", file=meta_out)

        species_graph = asm_summary.loc[accessions_graph, 'organism_name'].value_counts().sort_values(ascending=False)
        print(species_graph)
        for ix in species_graph.index:
            count = species_graph.loc[ix]
            print(f'"{ix}" = {count}', file=meta_out)

        print(file=meta_out)
        print("[align_set]", file=meta_out)
        print('fname = "align_seq.fna.gz"', file=meta_out)
        print(f"num_seqs = {len(align_ix)}", file=meta_out)
        print(f"avg_seq_len = {sum(lengths_align)/len(align_ix):.1f}", file=meta_out)

        avg_dist = dmatrix.iloc[align_ix, align_ix].mean(axis=None)
        avg_min_dist = dmatrix.iloc[align_ix, graph_ix].min(axis=0).mean(axis=None)
        print(f"avg_pairwise_dist = {avg_dist:g}", file=meta_out)
        print(f"avg_min_dist_to_graph = {avg_min_dist:g}", file=meta_out)
        print(f"avg_min_dist_to_graph = {avg_min_dist:g}")
        print(file=meta_out)

        print("[align_set.species]", file=meta_out)

        species_align = asm_summary.loc[accessions_align, 'organism_name'].value_counts().sort_values(ascending=False)
        for ix in species_align.index:
            count = species_align.loc[ix]
            print(f'"{ix}" = {count}', file=meta_out)


def load_assembly_summary_meta(fname):
    col_names = [
        'assembly_acc',
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

    return pandas.read_csv(fname, sep='\t', names=col_names, index_col=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip())

    parser.add_argument(
        'assembly_summary_list', type=Path,
        help="List to assembly summary metadata list as obtained from NCBI"
    )

    parser.add_argument(
        'sequences_fasta', type=Path,
        help="Path to FASTA with gene sequences"
    )
    parser.add_argument(
        'sequences_acc_map', type=Path,
        help="Path to a TSV mapping gene sequence IDs to their assembly accessions"
    )

    parser.add_argument(
        'mash_traingle', type=Path,
        help="Distance matrix for the gene sequences as created by `mash traingle`"
    )

    parser.add_argument(
        '-o', '--output-dir', type=Path,
        help="Output directory"
    )

    args = parser.parse_args()

    print("Loading assembly summary list...")
    asm_summary = load_assembly_summary_meta(args.assembly_summary_list)
    print(asm_summary)
    acc_map = pandas.read_csv(args.sequences_acc_map, sep='\t', index_col=0, names=['seq_id', 'accession'])

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading distance matrix...")
    dmatrix = phylip_to_matrix(args.mash_traingle)

    print("Creating linkage matrix...")
    Z = linkage(squareform(dmatrix.values), method='single')

    print("Deduplicating sequences...")
    duplicates = fcluster(Z, 0, criterion='distance')
    duplicate_clusters = numpy.unique(duplicates)

    to_keep = []
    for cluster in duplicate_clusters:
        cluster_select = duplicates == cluster
        cluster_size = numpy.count_nonzero(cluster_select)

        first = dmatrix.index[cluster_select][0]
        if cluster_size > 1:
            print(f"Found duplicates cluster (size: {cluster_size}), selecting", first)
            to_keep.append(first)
        else:
            to_keep.append(first)

    dmatrix_reduced = dmatrix.loc[to_keep, to_keep]
    print(dmatrix_reduced)
    print("Creating new linkage matrix...")
    Z = linkage(squareform(dmatrix_reduced.values), method='single')

    for max_dist in [0.01, 0.1, 0.25]:
        clusters_assigned = fcluster(Z, max_dist, criterion='distance')
        all_clusters, counts = numpy.unique(clusters_assigned, return_counts=True)
        print(f"Clustering with max dist = {max_dist}, num_clusters =", len(all_clusters))
        print("Cluster size histogram:")
        cluster_sizes, freq_hist = numpy.unique(counts, return_counts=True)
        max_freq = numpy.max(freq_hist)
        for size, freq in zip(cluster_sizes, freq_hist):
            print(f"{size}: {freq}", "#" * int(round((freq * 100) / max_freq)))

        dist_output = args.output_dir / f"maxdist_{max_dist:g}"
        dist_output.mkdir(parents=True, exist_ok=True)
        non_cluster_ids = []

        dataset_num = 1
        for cluster in all_clusters:
            cluster_select = clusters_assigned == cluster
            cluster_ix = numpy.nonzero(cluster_select)[0]
            cluster_size = numpy.count_nonzero(cluster_select)

            if cluster_size >= 1000:
                graph_sizes = [100, 250, 500]
            elif cluster_size >= 500:
                graph_sizes = [100, 250]
            elif cluster_size >= 100:
                graph_sizes = [5, 50]
            else:
                non_cluster_ids.extend(cluster_ix)
                continue

            for num_seq_graph in graph_sizes:
                graph_seq_select = numpy.random.choice(numpy.arange(len(cluster_ix)), size=num_seq_graph, replace=False)
                graph_seq_ix = cluster_ix[graph_seq_select]
                align_seq_select = numpy.ones((len(cluster_ix),), dtype=bool)
                align_seq_select[graph_seq_select] = False
                align_seq_ix = cluster_ix[align_seq_select]

                print("Creating dataset", dataset_num, "maxdist", max_dist, "graphsize", num_seq_graph, "align_seq",
                      len(align_seq_ix))
                create_dataset(asm_summary, acc_map, dist_output / f"dataset{dataset_num}_graphsize{num_seq_graph}",
                               dmatrix_reduced, args.sequences_fasta, max_dist, graph_seq_ix, align_seq_ix)

                dataset_num += 1

        with (gzip.open(dist_output / "non_cluster_seq.fna.gz", "wt") as out,
              gzip.open(args.sequences_fasta, "rt") as ifile):
            non_cluster_seqs = set(dmatrix_reduced.index[non_cluster_ids])

            for r in skbio.io.read(ifile, "fasta"):
                if r.metadata['id'] in non_cluster_seqs:
                    skbio.io.write(r, "fasta", into=out)


if __name__ == '__main__':
    main()
