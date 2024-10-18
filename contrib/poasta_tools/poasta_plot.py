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
from typing import Any

import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn
import networkx as nx


def load_visited_states(fname):
    visited_df = pandas.read_csv(fname, sep='\t')
    visited_df['node_pos'] = visited_df['offset'] - visited_df['diag']

    cols_to_convert = [
        ('node_pos', np.uint32),
        ('offset', np.uint32),
        ('diag', np.int32),
        ('score', np.uint32),
        ('node', np.uint32)
    ]

    for col, dtype in cols_to_convert:
        visited_df[col] = pandas.to_numeric(visited_df[col]).astype(dtype)

    visited_df = visited_df.sort_values(['state', 'node', 'node_pos', 'offset', 'score'])

    return visited_df


label_regexp = re.compile(r"<font.*>([\w\#\$]+)</font>")


def load_graph(fname: str | Path):
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

    # Remove any aligned interval edges
    edges_to_delete = []
    for u, v, data in g.edges(data=True):
        if data.get('style') == "dotted":
            edges_to_delete.append((u, v))

    g.remove_edges_from(edges_to_delete)

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

    return g, tsorted


def compute_node_y(g, tsorted) -> dict[int, int]:
    node_y = {int(tsorted[0]): 0}
    for n in tsorted[1:]:
        max_y_pred = max((pred for pred in map(int, g.predecessors(n))), key=node_y.get)
        node_y[int(n)] = node_y[max_y_pred] + len(g.nodes[str(max_y_pred)]['label'])

    return node_y


def filter_visited(node_y: dict[str, int], offset_max: int, visited_df: pandas.DataFrame,
                   x_range: str | None, y_range: str | None) -> pandas.DataFrame:
    max_y = max(node_y.values()) + 1

    to_return = visited_df
    x_start = 0
    x_end = offset_max
    print("orig shape:", to_return.shape)
    if x_range:
        x_start, x_end = map(int, x_range.split(':', maxsplit=1))
        if x_start < 0:
            x_start = offset_max + x_start
        if x_end < 0:
            x_end = offset_max + x_end

        to_return = to_return[(to_return['offset'] >= x_start) & (to_return['offset'] < x_end)]

    if y_range:
        y_start, y_end = map(int, y_range.split(':', maxsplit=1))
        if y_start < 0:
            y_start = max_y + y_start
        if y_end < 0:
            y_end = max_y + y_end

        to_return = to_return[(to_return['y'] >= y_start) & (to_return['y'] < y_end)]

    print("filtered shape:", to_return.shape)

    return to_return, x_start, x_end


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

    parser.add_argument('-x', '--x-range', default="", required=False, type=str,
                        help="Range of x values to plot. Accepts python slicing notation.")
    parser.add_argument('-y', '--y-range', default="", required=False, type=str,
                        help="Range of y values to plot. Accepts python slicing notation.")

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

    node_y = compute_node_y(g, tsorted)

    seq = g.graph['seq_to_align']
    matrix_seq = f"-{seq} "
    seq_len = len(seq)
    offset_max = seq_len + 2
    xlabels = [f"{matrix_seq[i]}\n{i}" for i in range(offset_max)]

    print("Loading visited states...")
    visited = load_visited_states(args.visited)
    visited['node_y'] = visited['node'].map(node_y)
    visited['y'] = visited['node_y'] + visited['node_pos']
    visited = visited.set_index(['state', 'node'])

    invalid_visited = visited[visited['offset'] > offset_max]
    if invalid_visited.shape[0] > 0:
        print("WARNING: some visited states are out of range!")
        print(invalid_visited)

    max_score = visited['score'].max()
    prefix = args.visited.stem.replace("_visited", "")

    visited, x_start, x_end = filter_visited(node_y, offset_max, visited, args.x_range, args.y_range)
    visited_nodes = set(visited.index.unique(level=1))
    min_node_rank = min(g.nodes[str(n)]['rank'] for n in visited_nodes)
    max_node_rank = max(g.nodes[str(n)]['rank'] for n in visited_nodes)
    valid_nodes = tsorted[min_node_rank:max_node_rank+1]
    print("Valid nodes", valid_nodes)

    node_height_ratios = [len(g.nodes[n]['label']) for n in valid_nodes]
    valid_nodes_len = sum(node_height_ratios)

    matrix_width = x_end - x_start
    xlabels = xlabels[x_start:x_end]
    fig_width = args.fig_width if args.fig_width is not None else int(round(matrix_width / 4))
    fig_height = int(round((fig_width * valid_nodes_len * 1.25) / matrix_width))

    align_states = ["Match", "Deletion", "Insertion"]
    for state in align_states:
        if state not in visited.index.unique(level=0):
            continue

        subdf = visited.loc[[state]]

        fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
        fig.get_layout_engine().set(hspace=0.00)
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
            node_df = subdf.loc[(state, [node_ix]), :]

            node_seq = g.nodes[str(node_ix)]['label']
            node_seq_len = len(node_seq)
            yticks = list(range(node_seq_len))
            yticklabels = [f"{node_seq[i]} - {i}" for i in yticks]

            indexed = node_df.reset_index().set_index(['node_pos', 'offset'])
            duplicated = indexed.index.duplicated(False)
            print("duplicated:")
            print(indexed[duplicated])
            deduplicated = indexed[~indexed.index.duplicated()]
            print("num to plot:", len(deduplicated))

            # pivot = node_df.pivot(index="node_pos", columns="offset", values="score")
            pivot = deduplicated.reset_index().pivot(index='node_pos', columns='offset', values='score')

            # Ensure rows/cols exists even when no data exists for those
            rows = pivot.index.union(list(range(node_seq_len)), sort=True)
            cols = pivot.columns.union(list(range(x_start, x_end)), sort=True)
            pivot = pivot.reindex(index=rows, columns=cols).sort_index()
            pivot.sort_index(axis=1, inplace=True)

            node_rank = g.nodes[str(node_ix)]['rank'] - min_node_rank
            print("min_node_rank", min_node_rank, node_rank)
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

        draw_graph(fig, graph_axes, node_axes, g, min_node_rank, max_node_rank)
        fig.savefig(args.output / f'{prefix}.{state}.pdf', bbox_inches='tight')


def draw_graph(fig, graph_axes, node_axes, g, min_node_rank, max_node_rank):
    # First compute layout with `dot` to determine node x positions
    pos = nx.nx_agraph.graphviz_layout(g, prog='dot')
    min_x = min(p[0] for p in pos.values())
    max_x = max(p[0] for p in pos.values())
    node_width = max(0.5, 0.2 * (max_x - min_x))
    node_width_half = node_width / 2

    fig_pos = {}
    for n, (x, y) in pos.items():
        rank = g.nodes[n]['rank']
        if rank < min_node_rank or rank > max_node_rank:
            continue

        rank = rank - min_node_rank
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
        ranku = g.nodes[u]['rank']
        rankv = g.nodes[v]['rank']

        if ranku < min_node_rank or ranku > max_node_rank:
            continue
        if rankv < min_node_rank or rankv > max_node_rank:
            continue

        ranku -= min_node_rank
        rankv -= min_node_rank

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


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
