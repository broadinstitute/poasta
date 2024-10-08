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

import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.transforms import IdentityTransform
import seaborn
import networkx as nx


def load_visited_states(fname):
    visited_df = pandas.read_csv(fname, sep='\t')
    visited_df['node_pos'] = visited_df['offset'] - visited_df['diag']

    for col in ['node_pos', 'offset', 'diag', 'score', 'node']:
        visited_df[col] = pandas.to_numeric(visited_df[col], downcast='signed')

    visited_df = visited_df.set_index(['state', 'node']).sort_index()

    return visited_df


label_regexp = re.compile(r"<font.*>([\w\#\$]+)</font>")


def load_graph(fname):
    g = nx.nx_agraph.read_dot(fname)

    # Parse any key-value pairs in the comments
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if not line.startswith('#'):
                break

            parts = line[1:].strip().split(':', maxsplit=1)
            if len(parts) == 2:
                k, v = parts
                g.graph[k.strip()] = v.strip()

    tsorted = list(nx.dfs_postorder_nodes(g))
    tsorted.reverse()

    g.graph['rankdir'] = 'TB'
    g.graph['graph']['rankdir'] = 'TB'

    for rank, n in enumerate(tsorted):
        # Remove GraphViz HTML formatting from node label
        seq = "".join(g.group(1) for g in label_regexp.finditer(g.nodes[n]['label']))
        g.nodes[n]['label'] = seq
        g.nodes[n]['rank'] = rank
        print(g.nodes[n])

    # Remove any aligned interval edges
    edges_to_delete = []
    for u, v, data in g.edges(data=True):
        if data.get('style') == "dotted":
            edges_to_delete.append((u, v))

    g.remove_edges_from(edges_to_delete)

    return g, tsorted


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


def main():
    parser = argparse.ArgumentParser(description="Plot POASTA aligner computation state")

    parser.add_argument('graph', type=Path,
                        help="The graph used for alignment in DOT format")
    parser.add_argument('visited', type=Path, metavar='VISITED_TSV',
                        help="TSV with visited states")
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help="Output directory")
    parser.add_argument('-w', '--fig-width', default=None, required=False, type=int)

    args = parser.parse_args()

    if args.output is None:
        print("ERROR: no output directory specified!", file=sys.stderr)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)

    print("Loading graph...", file=sys.stderr)
    g, tsorted = load_graph(args.graph)

    if not g.graph.get('seq_to_align'):
        print("No sequence to align metadata in graph file...")
        return 1

    seq = g.graph['seq_to_align']
    matrix_seq = f"-{seq} "
    seq_len = len(seq)
    matrix_width = seq_len + 2
    total_node_seq_len = sum(len(g.nodes[n]['label']) for n in tsorted)
    xlabels = [f"{matrix_seq[i]}\n{i}" for i in range(matrix_width)]

    fig_width = args.fig_width if args.fig_width is not None else int(round(matrix_width / 4))
    fig_height = int(round((fig_width * total_node_seq_len) / matrix_width))
    node_height_ratios = [len(g.nodes[n]['label']) for n in tsorted]

    print("Loading visited states...")
    visited = load_visited_states(args.visited)
    max_score = visited['score'].max()
    prefix = args.visited.stem.replace("_visited", "")

    align_states = ["Match", "Deletion", "Insertion"]
    for state in align_states:
        subdf = visited.loc[[state]]

        fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
        gs = fig.add_gridspec(ncols=3, nrows=len(node_height_ratios),
                              height_ratios=node_height_ratios,
                              width_ratios=[1, 8, 0.25])

        graph_axes = [
            fig.add_subplot(gs[i, 0]) for i in range(len(node_height_ratios))
        ]
        node_axes = [
            fig.add_subplot(gs[i, 1]) for i in range(len(node_height_ratios))
        ]
        cbar_ax = fig.add_subplot(gs[:, 2])

        for i in range(0, len(node_height_ratios)-1):
            graph_axes[i].sharex(graph_axes[-1])
            node_axes[i].sharex(node_axes[-1])

        for node_ix in subdf.index.unique(level=1):
            node_df = subdf.loc[(state, [node_ix]), :].sort_values('score')
            node_seq = g.nodes[str(node_ix)]['label']
            node_seq_len = len(node_seq)
            node_rank = g.nodes[str(node_ix)]['rank']
            yticks = list(range(node_seq_len))
            yticklabels = [f"{node_seq[i]} - {i}" for i in yticks]

            indexed = node_df.reset_index().set_index(['node_pos', 'offset'])
            duplicated = indexed.index.duplicated(False)
            print("duplicated:")
            print(indexed[duplicated].sort_index())
            deduplicated = indexed[~indexed.index.duplicated()]

            # pivot = node_df.pivot(index="node_pos", columns="offset", values="score")
            pivot = deduplicated.reset_index().pivot(index='node_pos', columns='offset', values='score')

            # Ensure rows/cols exists even when no data exists for those
            rows = pivot.index.union(list(range(node_seq_len)), sort=True)
            cols = pivot.columns.union(list(range(matrix_width)), sort=True)
            pivot = pivot.reindex(index=rows, columns=cols).sort_index()
            pivot.sort_index(axis=1, inplace=True)

            node_axes[node_rank].grid('both')
            seaborn.heatmap(pivot,
                            cmap="turbo", vmin=0, vmax=max_score,
                            xticklabels=xlabels, yticklabels=yticklabels,
                            annot=True, fmt="g", annot_kws={"fontsize": "x-small"},
                            ax=node_axes[node_rank], cbar_ax=cbar_ax)

            node_axes[node_rank].set_yticklabels(yticklabels, rotation=0)
            node_axes[node_rank].set_ylabel('')

        node_axes[-1].set(xticks=list(range(matrix_width)), xticklabels=xlabels)
        for i, ax in enumerate(node_axes):
            if i < len(node_height_ratios) - 1:
                ax.xaxis.set_tick_params(labelbottom=False)
                ax.set_xlabel('')

        draw_graph(fig, graph_axes, node_axes, g)
        fig.savefig(args.output / f'{prefix}.{state}.pdf', bbox_inches='tight')


def draw_graph(fig, graph_axes, node_axes, g):
    # First compute layout with `dot` to determine node x positions
    pos = nx.nx_agraph.graphviz_layout(g, prog='dot')
    min_x = min(p[0] for p in pos.values())
    max_x = max(p[0] for p in pos.values())
    node_width = max(0.5, 0.2 * (max_x - min_x))
    node_width_half = node_width / 2

    fig_pos = {}
    for n, (x, y) in pos.items():
        rank = g.nodes[n]['rank']
        graph_axes[rank].set_xlim(min_x - node_width_half, max_x + node_width_half)
        graph_axes[rank].set_xticks([])
        graph_axes[rank].set_yticks([])
        graph_axes[rank].spines[['left', 'right', 'top', 'bottom']].set_visible(False)
        graph_axes[rank].patch.set_alpha(0.0)

        trans = graph_axes[rank].get_xaxis_transform()
        node_x = x - node_width_half

        patch = Rectangle((node_x, 0), node_width, 1, facecolor='grey', transform=trans)
        graph_axes[rank].add_patch(patch)
        graph_axes[rank].invert_yaxis()

        fig_pos[n] = [
            trans.transform((x, 0)),
            trans.transform((x, 1))
        ]

    for (u, v) in g.edges():
        print((u, v))
        ranku = g.nodes[u]['rank']
        rankv = g.nodes[v]['rank']

        axu = graph_axes[ranku]
        transu = axu.get_xaxis_transform()
        axv = graph_axes[rankv]
        transv = axv.get_xaxis_transform()

        ux = pos[u][0]
        vx = pos[v][0]
        axu.annotate(
            "",
            xytext=(ux, 0),
            textcoords=transu,
            xy=(vx, 1),
            xycoords=transv,
            arrowprops={"facecolor": "black", "width": 0.5, "headwidth": 8},
            in_layout=False
        )


    # ax.set_axis_on()
    # ax.yaxis.set_visible(True)
    # ax.invert_yaxis()
    # ax.tick_params(axis="y", left=True, labelleft=True)
    # ax.set_ylim(0, len(graph_layout_data.ylabels))
    # ax.set_yticks(graph_layout_data.yticks)
    # ax.set_yticklabels(graph_layout_data.ylabels)
    # ax.set_ylabel("Rank (node)")


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
