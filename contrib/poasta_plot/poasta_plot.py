#!/usr/bin/env python3
"""
Create a plot showing the computed cells by POASTA.

The plot will show both the graph (read from the accompanying DOT file),
and the DP matrix.
"""

import re
import sys
import json
import argparse
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
from collections import defaultdict

import numpy
import pandas
import matplotlib.pyplot as plt
import seaborn
import networkx as nx


def load_spoa_matrix(fname):
    """
    This function loads a dynamic programming matrix as computed by SPOA
    from a TSV file
    """

    xlabels = []
    ylabels = []
    row_data = []
    col_data = []
    score_data = []

    with open(fname) as f:
        for i, line in enumerate(f):
            line = line.strip()

            if not line:
                continue

            parts = line.split('\t')
            if i == 0:
                xlabels = [f"{pos}\n{c}" for pos, c in enumerate(parts)]
            else:
                ylabels.append(f"{parts[0]} ({parts[1]})")
                for col, score in enumerate(parts[2:]):
                    row_data.append(i - 1)
                    col_data.append(col)
                    score_data.append(int(score) if score != "2147482624" else float('NaN'))

    return xlabels, ylabels, pandas.DataFrame({"rank": row_data, "offset": col_data, "score": score_data})


poasta_node_label = re.compile(r"'(\w)' \((\d+)\)")
spoa_node_label = re.compile(r"(\d+) - (\w)")


def load_graph(fname):
    g = nx.nx_agraph.read_dot(fname)

    g.graph['rankdir'] = 'TB'
    g.graph['graph']['rankdir'] = 'TB'

    for n, ndata in g.nodes(data=True):
        if (match := poasta_node_label.search(ndata['label'])):
            ndata['rank'] = int(match.group(2))
            ndata['symbol'] = match.group(1)
        elif (match := spoa_node_label.search(ndata['label'])):
            ndata['rank'] = int(match.group(1)) + 1
            ndata['symbol'] = match.group(2)
        else:
            print("Could not parse node label:", ndata['label'], file=sys.stderr)

    edges_to_delete = []
    for u, v, data in g.edges(data=True):
        if data.get('style') == "dotted":
            edges_to_delete.append((u, v))

    g.remove_edges_from(edges_to_delete)

    return g


def poa_graph_layout(g):
    """
    Use graphviz to layout the x position of a node, but use its rank for the y position
    """

    pos = nx.nx_agraph.graphviz_layout(g, prog='dot')

    new_pos = {
        n: (pos[n][0], ndata['rank'] + 0.5) for n, ndata in g.nodes(data=True)
    }

    return new_pos


def poa_matrix_discontinuieties(g):
    hlines = set()
    for n in nx.topological_sort(g):
        if g.out_degree(n) > 1:
            for neighbor in g.successors(n):
                hlines.add(g.nodes[neighbor]['rank'])
        elif g.in_degree(n) > 1 and g.out_degree(n) > 0:
            hlines.add(g.nodes[n]['rank'])

    return list(sorted(hlines))


@dataclass
class PoastaIteration:
    score: Optional[int] = 0
    extend_paths: Optional[list[tuple]] = field(default_factory=list)
    wavefronts: Optional[dict[str, dict[str, list[int | bool]]]] = field(default_factory=dict)

    def to_dataframe(self, kind):
        nodes = []
        offsets = []
        is_endpoints = []

        if kind == "deletion":
            nodes.append(0)
            offsets.append(None)
            is_endpoints.append(False)

        for rank, (offset, is_endpoint) in self.wavefronts[kind].items():
            nodes.append(rank)
            offsets.append(offset)
            is_endpoints.append(is_endpoint)

        return pandas.DataFrame({
            "score": pandas.Series([self.score] * len(nodes), dtype=int),
            "kind": pandas.Series([kind] * len(nodes)),
            "rank": pandas.Series(nodes, dtype=int),
            "offset": pandas.Series(offsets, dtype=float),
            "is_endpoint": pandas.Series(is_endpoints, dtype=bool),
        })

    def to_dataframe_after_extend(self, kind):
        if kind == "match":
            # Apply offsets from path extension
            data = self.wavefronts.get(kind, {}).copy()
            for path in self.extend_paths:
                for node, offset in path[:-1]:
                    data[int(node)] = (int(offset), False)

                data[int(path[-1][0])] = (int(path[-1][1]), True)

        else:
            data = self.wavefronts.get(kind, {})

        nodes = []
        offsets = []
        is_endpoints = []

        if kind == "deletion":
            nodes.append(0)
            offsets.append(None)
            is_endpoints.append(False)

        for rank, (offset, is_endpoint) in data.items():
            nodes.append(rank)
            offsets.append(offset)
            is_endpoints.append(is_endpoint)

        return pandas.DataFrame({
            "score": pandas.Series([self.score] * len(nodes), dtype=int),
            "kind": pandas.Series([kind] * len(nodes)),
            "rank": pandas.Series(nodes, dtype=int),
            "offset": pandas.Series(offsets, dtype=float),
            "is_endpoint": pandas.Series(is_endpoints, dtype=bool),
        })


@dataclass
class PoastaAlignment:
    seq_name: Optional[str] = None
    seq_length: Optional[int] = 0
    max_rank: Optional[int] = 0
    iterations: list[PoastaIteration] = field(default_factory=list)


def load_poasta_log(logfile):
    alignment = PoastaAlignment()
    curr_iteration = PoastaIteration(score=0)
    curr_iteration.wavefronts["match"] = {0: (0, True)}

    with open(logfile) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            entry = json.loads(line)
            kind, data = next(iter(entry.items()))

            if kind == "NewSequence":
                alignment.seq_name = data['seq_name']
                alignment.seq_length = data['seq_length']
                alignment.max_rank = data['max_rank']
            elif kind == "ExtendPath":
                curr_iteration.extend_paths.append([data['start'], *data['path']])
            elif kind == "CurrWavefront":
                if curr_iteration.score != data['score']:
                    alignment.iterations.append(curr_iteration)
                    curr_iteration = PoastaIteration(score=int(data['score']))

                curr_iteration.wavefronts[data['wf_type']] = defaultdict(dict)

                for rank, offset, is_endpoint in data['node_offsets']:
                    curr_iteration.wavefronts[data['wf_type']][int(rank)] = (int(offset), bool(is_endpoint))

        alignment.iterations.append(curr_iteration)

    return alignment


def main():
    parser = argparse.ArgumentParser(description="Plot POASTA aligner computation state")

    parser.add_argument('graph', type=Path,
                        help="The graph used for alignment in DOT format")
    parser.add_argument('poasta_log', type=Path,
                        help="The POASTA debug log for a specific sequence")
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help="Output directory")

    args = parser.parse_args()

    if args.output is None:
        print("ERROR: no output directory specified!", file=sys.stderr)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)

    print("Loading graph...", file=sys.stderr)
    g = load_graph(args.graph)
    graph_layout = poa_graph_layout(g)

    print("Loading poasta data...", file=sys.stderr)
    alignment = load_poasta_log(args.poasta_log)

    print("Creating plots...", file=sys.stderr)
    iter_dfs_after_ext = defaultdict(list)

    max_score = alignment.iterations[-1].score
    xlabels = list(range(alignment.seq_length + 1))
    yticks = numpy.arange(alignment.max_rank) + 0.5
    ylabels = ["0 (S)", *(f"{ndata['rank']} ({ndata['symbol']})" for n, ndata in g.nodes(data=True))]
    hlines = poa_matrix_discontinuieties(g)

    xlabels_score = [str(i.score) for i in alignment.iterations]

    for aln_iter in alignment.iterations:
        fig, axes = plt.subplots(3, 5, figsize=(20, 20), width_ratios=[2, 8.75, 0.25, 8.75, 0.25],
                                 constrained_layout=True)

        for row, kind in enumerate(["insertion", "deletion", "match"]):
            nx.draw(g, pos=graph_layout, ax=axes[row, 0], node_size=75,
                    labels={n: ndata['symbol'] for n, ndata in g.nodes(data=True)},
                    font_size=6)

            df_after_ext = aln_iter.to_dataframe_after_extend(kind)
            iter_dfs_after_ext[kind].append(df_after_ext)
            df = pandas.concat(iter_dfs_after_ext[kind], ignore_index=True).set_index(['rank', 'offset'])

            print("Processing score", aln_iter.score, f"({kind})", file=sys.stderr)
            if len(df) > 0:
                # Only keep cells with lowest score
                deduplicated = df[~df.index.duplicated(keep='first')].reset_index()
                seaborn.heatmap(deduplicated.pivot(index="rank", columns="offset", values="score"),
                                cmap="plasma", vmin=0, vmax=max_score,
                                xticklabels=xlabels, yticklabels=ylabels,
                                annot=True, fmt="g", annot_kws={"fontsize": "x-small"},
                                ax=axes[row, 1], cbar_ax=axes[row, 2])

            if kind == "match":
                for path in aln_iter.extend_paths:
                    path_x = numpy.array([p[1] + 0.6 for p in path])
                    path_y = numpy.array([p[0] + 0.6 for p in path])

                    axes[row, 1].plot(path_x, path_y, color='white', marker='o', alpha=0.5)

            for y in hlines:
                axes[row, 1].axhline(y, color='black')

            axes[row, 0].set_axis_on()
            axes[row, 0].yaxis.set_visible(True)
            axes[row, 0].invert_yaxis()
            axes[row, 0].tick_params(axis="y", left=True, labelleft=True)
            axes[row, 0].set_ylim(0, len(ylabels))
            axes[row, 0].set_yticks(numpy.arange(len(ylabels)) + 0.5)
            axes[row, 0].set_yticklabels(ylabels)
            axes[row, 0].set_ylabel("Rank (node)")
            axes[row, 1].set_title(f"{kind.upper()} matrix")
            axes[row, 1].sharey(axes[row, 0])
            axes[row, 1].set_yticks(yticks)
            axes[row, 1].set_yticklabels(ylabels)
            axes[row, 1].tick_params(left=False, labelleft=False)
            axes[row, 1].set_ylabel("")
            axes[row, 1].set_ylim(0, len(ylabels))
            axes[row, 1].invert_yaxis()
            axes[row, 2].set_ylabel("Score")
            axes[row, 2].yaxis.set_ticks_position('left')
            axes[row, 2].yaxis.set_label_position('left')

            if len(df) > 0:
                seaborn.heatmap(df.reset_index().pivot(index="rank", columns="score", values="offset"),
                                cmap="plasma", vmin=0, vmax=alignment.seq_length + 1,
                                xticklabels=xlabels_score,
                                annot=True, fmt="g", annot_kws={"fontsize": "x-small"},
                                ax=axes[row, 3], cbar_ax=axes[row, 4])

            axes[row, 3].set_title(f"{kind.upper()} wavefront")
            axes[row, 3].set_ylim(0, len(ylabels))
            axes[row, 3].sharey(axes[row, 0])
            axes[row, 3].set_ylabel("")
            axes[row, 3].tick_params(left=True, labelleft=False)
            axes[row, 4].set_ylabel("Offset")

            # Highlight alignment end points
            scatter_x = []
            scatter_y = []
            for record in df_after_ext[df_after_ext['is_endpoint']].itertuples():
                scatter_x.append(record.offset + 0.5)
                scatter_y.append(record.rank + 0.5)

            axes[row, 1].scatter(scatter_x, scatter_y, s=150, facecolors='none', edgecolors='white', alpha=0.5)

            scatter_x = []
            scatter_y = []
            df = df.reset_index()
            score_to_x = {s: x for x, s in enumerate(df['score'].unique())}
            for record in df[df['is_endpoint']].itertuples():
                scatter_x.append(score_to_x[record.score] + 0.5)
                scatter_y.append(record.rank + 0.5)

            axes[row, 3].scatter(scatter_x, scatter_y, s=100, facecolors='none', edgecolors='white', alpha=0.5)

        fig.savefig(args.output / "score{}.after_extend.png".format(aln_iter.score), dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
