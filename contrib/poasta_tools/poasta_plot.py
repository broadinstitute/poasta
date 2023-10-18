#!/usr/bin/env python3
"""
Create a plot showing the computed cells by POASTA.

The plot will show both the graph (read from the accompanying DOT file),
and the DP matrix.
"""

import re
import sys
import argparse
from pathlib import Path
from typing import NamedTuple, Any
from collections import defaultdict

import numpy
import numpy as np
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


poasta_node_label = re.compile(r"'(\w|#)' \((\d+)\)")


def load_graph(fname):
    g = nx.nx_agraph.read_dot(fname)
    tsorted = list(nx.dfs_postorder_nodes(g))
    tsorted.reverse()

    g.graph['rankdir'] = 'TB'
    g.graph['graph']['rankdir'] = 'TB'

    node_ix_to_rank = {}

    for rank, n in enumerate(tsorted):
        ndata = g.nodes[n]
        ndata['rank'] = rank

        if match := poasta_node_label.search(ndata['label']):
            node_ix = int(match.group(2))
            ndata['node_ix'] = node_ix
            node_ix_to_rank[node_ix] = rank

            ndata['symbol'] = match.group(1)
        else:
            print("Could not parse node label:", ndata['label'], file=sys.stderr)

    edges_to_delete = []
    for u, v, data in g.edges(data=True):
        if data.get('style') == "dotted":
            edges_to_delete.append((u, v))

    g.remove_edges_from(edges_to_delete)

    return g, node_ix_to_rank


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


def load_state_graph_data(fname, score):
    kv = {}
    with open(fname) as f:
        first = next(f)[2:]
        fields = first.split(" - ")

        for field in fields:
            k, v = field.split(':', maxsplit=1)
            k = k.strip()
            v = v.strip()

            if v.isnumeric():
                v = int(v)

            kv[k] = v

    df = pandas.read_csv(fname, sep='\t', comment='#')
    df['queued'] = df['score'] > score

    return df.set_index('matrix').sort_index(), kv


score_regexp = re.compile(r"score(\d+)")


class GraphLayoutData(NamedTuple):
    node_pos: dict[Any, tuple[float | int, float | int]]
    ylabels: list[str]
    yticks: list[float | int] | np.ndarray
    hlines: list[float | int]


def main():
    parser = argparse.ArgumentParser(description="Plot POASTA aligner computation state")

    parser.add_argument('graph', type=Path,
                        help="The graph used for alignment in DOT format")
    parser.add_argument('state_graph_tsvs', type=Path, metavar='STATE_TSV', nargs='+',
                        help="List of state graph TSVs")
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help="Output directory")
    parser.add_argument('-w', '--fig-width', default=None, required=False, type=int)

    args = parser.parse_args()

    if args.output is None:
        print("ERROR: no output directory specified!", file=sys.stderr)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)

    print("Loading graph...", file=sys.stderr)
    g, node_ix_to_rank = load_graph(args.graph)

    print("Creating graph layout...")
    graph_layout = poa_graph_layout(g)

    nodes_sorted = [n for n in sorted(g.nodes(data=True), key=lambda d: d[1]['rank'])]
    ylabels = [f"{ndata['rank']} ({ndata['symbol']})" for n, ndata in nodes_sorted]
    yticks = numpy.arange(len(ylabels)) + 0.5
    hlines = poa_matrix_discontinuieties(g)

    graph_layout_data = GraphLayoutData(graph_layout, ylabels, yticks, hlines)

    print("Loading poasta data...", file=sys.stderr)
    state_graph_scores = {
        fname.name: int(score_regexp.search(fname.name).group(1))
        for fname in args.state_graph_tsvs
    }

    sorted_state_graphs: list[Path] = list(sorted(args.state_graph_tsvs, key=lambda v: state_graph_scores[v.name]))

    if not sorted_state_graphs:
        print("No state graph data given!", file=sys.stderr)
        return 1

    max_score = max(state_graph_scores.values())

    for state_graph_path in sorted_state_graphs:
        score = state_graph_scores[state_graph_path.name]
        state_df, kv = load_state_graph_data(state_graph_path, score)
        state_df['rank'] = state_df['node_id'].map(node_ix_to_rank).astype(int)

        print("SCORE:", score)
        print(state_graph_path)
        print(state_df)
        print(kv)

        create_plot_for_score(g, graph_layout_data, kv['seq'], state_df, score, max_score,
                              args.output, args.fig_width)
        print()


def create_plot_for_score(g, graph_layout_data, sequence, state_df, score, max_score, output_dir, fig_width=None):

    num_nodes = len(g)
    max_offset = len(sequence) + 1
    seq_chars = ["-"] + list(sequence)
    xlabels = [f"{seq_chars[i]}\n{i}" for i in range(max_offset)]

    fig_width = fig_width if fig_width is not None else int(round(len(xlabels) / 4))
    fig_height = int(round((fig_width * num_nodes) / max_offset))

    for row, kind in enumerate(["match", "deletion", "insertion"]):
        score_df = state_df.loc[[kind]].copy()
        queued_df = (score_df[score_df['score'] > score]
                     .set_index(['rank', 'offset']))
        score_df = (score_df[score_df['score'] <= score]
                    .set_index(['rank', 'offset']))

        print("Processing score", score, f"({kind})", file=sys.stderr)
        if len(score_df) > 0:
            # Only keep cells with lowest score
            # duplicated = score_df.index.duplicated(keep=False)
            # if numpy.count_nonzero(duplicated) > 0:
            #     print("DUPLICATED:", file=sys.stderr)
            #     print(score_df[duplicated], file=sys.stderr)

            # deduplicated = df[~df.index.duplicated(keep='first')].reset_index()
            deduplicated = score_df.reset_index()

            # Make heatmap DF
            pivot = deduplicated.pivot(index="rank", columns="offset", values="score")
            pivot_queued = queued_df.reset_index().pivot(index="rank", columns="offset", values="score")

            # Ensure rows/cols exists even when no data exists for those
            rows = pivot.index.union(list(range(num_nodes + 1)), sort=True)
            cols = pivot.columns.union(list(range(max_offset)), sort=True)
            pivot = pivot.reindex(index=rows, columns=cols).sort_index()
            pivot.sort_index(axis=1, inplace=True)
            pivot_queued = pivot_queued.reindex(index=rows, columns=cols).sort_index()
            pivot_queued.sort_index(axis=1, inplace=True)

            if kind == "match":
                print(score_df)
                print(pivot.index)
                print(pivot.columns)
                print(pivot)

            fig, axes = make_dp_matrix_plot(g, graph_layout_data, pivot, pivot_queued, max_score,
                                            xlabels, figsize=(fig_width, fig_height))

            fig.savefig(output_dir / f"score{score}.{kind}.png", dpi=150)
            plt.close(fig)


def make_dp_matrix_plot(g, graph_layout_data: GraphLayoutData, matrix, matrix_queued, max_score, xlabels, figsize):
    fig, axes = plt.subplots(1, 3, figsize=figsize, width_ratios=[1, 8.75, 0.25],
                             constrained_layout=True)

    nx.draw(g, pos=graph_layout_data.node_pos, ax=axes[0], node_size=75,
            labels={n: ndata['symbol'] for n, ndata in g.nodes(data=True)},
            font_size=6)

    seaborn.heatmap(matrix,
                    cmap="plasma", vmin=0, vmax=max_score,
                    xticklabels=xlabels, yticklabels=graph_layout_data.ylabels,
                    annot=True, fmt="g", annot_kws={"fontsize": "x-small"},
                    ax=axes[1], cbar_ax=axes[2])

    seaborn.heatmap(matrix_queued,
                    cmap="plasma", vmin=0, vmax=max_score,
                    xticklabels=xlabels, yticklabels=graph_layout_data.ylabels,
                    annot=True, fmt="g", annot_kws={"fontsize": "x-small"}, alpha=0.4,
                    ax=axes[1], cbar_ax=axes[2])

    for y in graph_layout_data.hlines:
        axes[1].axhline(y, color='black')

    axes[0].set_axis_on()
    axes[0].yaxis.set_visible(True)
    axes[0].invert_yaxis()
    axes[0].tick_params(axis="y", left=True, labelleft=True)
    axes[0].set_ylim(0, len(graph_layout_data.ylabels))
    axes[0].set_yticks(graph_layout_data.yticks)
    axes[0].set_yticklabels(graph_layout_data.ylabels)
    axes[0].set_ylabel("Rank (node)")

    axes[1].set_xticklabels(xlabels, rotation=0)
    axes[1].sharey(axes[0])
    axes[1].set_yticks(graph_layout_data.yticks)
    axes[1].set_yticklabels(graph_layout_data.ylabels)
    axes[1].grid(True, which='both')
    axes[1].set_axisbelow(True)
    axes[1].tick_params(left=False, labelleft=False)
    axes[1].set_ylabel("")
    axes[1].set_ylim(0, len(graph_layout_data.ylabels))
    axes[1].invert_yaxis()

    axes[2].set_ylabel("Score")

    return fig, axes


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
