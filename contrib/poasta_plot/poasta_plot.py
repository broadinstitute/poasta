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


poasta_node_label = re.compile(r"'(\w|#)' \((\d+)\)")


def load_graph(fname):
    g = nx.nx_agraph.read_dot(fname)
    tsorted = list(nx.topological_sort(g))

    g.graph['rankdir'] = 'TB'
    g.graph['graph']['rankdir'] = 'TB'

    node_ix_to_rank = {}

    for rank, n in enumerate(tsorted):
        print(rank, n)
        ndata = g.nodes[n]
        ndata['rank'] = rank

        if (match := poasta_node_label.search(ndata['label'])):
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


def state_tree_to_df(tree):
    root = next(nx.topological_sort(tree))

    print("root", tree.nodes[root])

    data = [
        {"graph_ix": tree.nodes[root]['graph_node_ix'],
         "offset": tree.nodes[root]['offset'],
         "matrix": "Match",
         "score": 0}
    ]
    scores = {root: 0}

    for prev, curr in nx.bfs_edges(tree, root):
        prev_offset = tree.nodes[prev]['offset']
        prev_state = tree.nodes[prev]['state']
        prev_score = scores[prev]

        curr_offset = tree.nodes[curr]['offset']
        curr_state = tree.nodes[curr]['state']

        match prev_state:
            case "Start" | "Match" | "Mismatch":
                score_delta = {
                    "Match": 0,
                    "Mismatch": 4,
                    "Deletion": 8,
                    "Insertion": 8
                }
            case "Deletion":
                score_delta = {
                    "Match": 0,
                    "Deletion": 2,
                }
            case "Insertion":
                score_delta = {
                    "Match": 0,
                    "Insertion": 2,
                }

        curr_score = prev_score + score_delta[curr_state]
        scores[curr] = curr_score

        if curr_state == "Match":
            edge_data = tree.edges[prev, curr]
            o = prev_offset + 1

            if edge_data['type'] == "extend_matches" and edge_data['ext_matching_nodes'] != "_networkx_list_start":
                for n in edge_data['ext_matching_nodes']:
                    data.append({
                        "graph_ix": n,
                        "offset": o,
                        "score": curr_score,
                        "matrix": "Match"
                    })

                    o += 1

        data.append({
            "graph_ix": tree.nodes[curr]['graph_node_ix'],
            "offset": curr_offset,
            "matrix": {
                "Start": "Match",
                "Mismatch": "Match",
            }.get(curr_state, curr_state),
            "score": curr_score
        })

    return pandas.DataFrame(data).set_index('matrix')



def main():
    parser = argparse.ArgumentParser(description="Plot POASTA aligner computation state")

    parser.add_argument('graph', type=Path,
                        help="The graph used for alignment in DOT format")
    parser.add_argument('poasta_aln_tree', type=Path,
                        help="The POASTA debug log for a specific sequence")
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help="Output directory")

    args = parser.parse_args()

    if args.output is None:
        print("ERROR: no output directory specified!", file=sys.stderr)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)

    print("Loading graph...", file=sys.stderr)
    g, node_ix_to_rank = load_graph(args.graph)
    graph_layout = poa_graph_layout(g)

    print("Loading poasta data...", file=sys.stderr)
    tree = nx.read_gml(args.poasta_aln_tree)
    tree_score_df = state_tree_to_df(tree)
    print(tree_score_df['graph_ix'].map(node_ix_to_rank))
    for gix in tree_score_df['graph_ix'].unique():
        if gix not in node_ix_to_rank:
            print(gix)

    tree_score_df['rank'] = tree_score_df['graph_ix'].map(node_ix_to_rank).astype(int)

    print("Creating plots...", file=sys.stderr)
    max_score = int(tree_score_df['score'].max())
    max_rank = int(tree_score_df['rank'].max()) + 1
    max_offset = int(tree_score_df['offset'].max()) + 1
    xlabels = list(range(max_offset))
    yticks = numpy.arange(max_rank) + 0.5

    nodes_sorted = [n for n in sorted(g.nodes(data=True), key=lambda d: d[1]['rank'])]
    ylabels = [f"{ndata['rank']} ({ndata['symbol']})" for n, ndata in nodes_sorted]
    hlines = poa_matrix_discontinuieties(g)

    xlabels_score = [f"{score:g}" for score in tree_score_df['score'].unique()]

    for score in tree_score_df['score'].unique():
        fig, axes = plt.subplots(1, 3, figsize=(8, 6), width_ratios=[1, 8.75, 0.25],
                                 constrained_layout=True)

        for row, kind in enumerate(["Match"]):
            nx.draw(g, pos=graph_layout, ax=axes[0], node_size=75,
                    labels={n: ndata['symbol'] for n, ndata in g.nodes(data=True)},
                    font_size=6)

            score_df = tree_score_df.loc[kind].copy()
            score_df = (score_df[score_df['score'] <= score]
                        .set_index(['rank', 'offset']))

            print("Processing score", score, f"({kind})", file=sys.stderr)
            if len(score_df) > 0:
                # Only keep cells with lowest score
                duplicated = score_df.index.duplicated(keep=False)
                if numpy.count_nonzero(duplicated) > 0:
                    print("DUPLICATED:", file=sys.stderr)
                    print(score_df[duplicated], file=sys.stderr)

                # deduplicated = df[~df.index.duplicated(keep='first')].reset_index()
                deduplicated = score_df.reset_index()

                # Make heatmap DF
                pivot = deduplicated.pivot(index="rank", columns="offset", values="score")

                # Ensure rows/cols exists even when no data exists for those
                rows = pivot.index.union(list(range(max_rank + 1)), sort=True)
                cols = pivot.columns.union(list(range(max_offset)), sort=True)

                pivot = pivot.reindex(index=rows, columns=cols).sort_index()
                if kind == "Match":
                    print(score_df)
                    print(pivot.index)
                    print(pivot.columns)
                    print(pivot)

                seaborn.heatmap(pivot,
                                cmap="plasma", vmin=0, vmax=max_score,
                                xticklabels=xlabels, yticklabels=ylabels,
                                annot=True, fmt="g", annot_kws={"fontsize": "x-small"},
                                ax=axes[1], cbar_ax=axes[2])

            for y in hlines:
                axes[1].axhline(y, color='black')

            axes[0].set_axis_on()
            axes[0].yaxis.set_visible(True)
            axes[0].invert_yaxis()
            axes[0].tick_params(axis="y", left=True, labelleft=True)
            axes[0].set_ylim(0, len(ylabels))
            axes[0].set_yticks(yticks)
            axes[0].set_yticklabels(ylabels)
            axes[0].set_ylabel("Rank (node)")
            axes[1].set_title(f"{kind.upper()} matrix")
            axes[1].sharey(axes[0])
            axes[1].set_yticks(yticks)
            axes[1].set_yticklabels(ylabels)
            axes[1].tick_params(left=False, labelleft=False)
            axes[1].set_ylabel("")
            axes[1].set_ylim(0, len(ylabels))
            axes[1].invert_yaxis()
            axes[2].set_ylabel("Score")

        fig.savefig(args.output / "score{}.png".format(score), dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
