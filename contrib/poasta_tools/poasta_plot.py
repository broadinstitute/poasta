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
from typing import Optional, NamedTuple

import matplotlib.colors
import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
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

    g.graph['rankdir'] = 'TB'
    g.graph['graph']['rankdir'] = 'TB'

    for n, ndata in g.nodes(data=True):
        if match := poasta_node_label.search(ndata['label']):
            rank = int(match.group(2))
            ndata['rank'] = rank
            ndata['symbol'] = match.group(1)
        else:
            print("Could not parse node label:", ndata['label'], file=sys.stderr)

    edges_to_delete = []
    for u, v, data in g.edges(data=True):
        if data.get('style') == "dotted":
            edges_to_delete.append((u, v))

    g.remove_edges_from(edges_to_delete)

    tsorted = [n for n, ndata in sorted(g.nodes(data=True), key=lambda v: v[1]['rank'])]

    return g, tsorted


def poa_graph_layout(g):
    """
    Use graphviz to layout the x position of a node, but use its rank for the y position
    """

    pos = nx.nx_agraph.graphviz_layout(g, prog='dot')

    new_pos = {
        n: (pos[n][0], ndata['rank'] + 0.5) for n, ndata in g.nodes(data=True)
    }

    return new_pos


def poa_matrix_discontinuieties(g, tsorted):
    hlines = set()
    for n in tsorted:
        if g.out_degree(n) > 1:
            for neighbor in g.successors(n):
                hlines.add(g.nodes[neighbor]['rank'])
        elif g.in_degree(n) > 1 and g.out_degree(n) > 0:
            hlines.add(g.nodes[n]['rank'])

    return list(sorted(hlines))


class Backtrace(NamedTuple):
    prev_diag: Optional[int]
    prev_state: str


class Interval:
    def __init__(self, start: int, stop: int, extendable: bool, backtrace: Backtrace):
        self.start = start
        self.stop = stop
        self.extendable = extendable
        self.backtrace = backtrace

    @classmethod
    def from_json(cls, json_obj):
        start = json_obj['interval'][0]['start']
        stop = json_obj['interval'][0]['end']
        extendable = json_obj['interval'][1]['extendable']
        backtrace = Backtrace(json_obj['interval'][1]['bt']['prev_diag'],
                              json_obj['interval'][1]['bt']['prev_state'])

        return cls(start, stop, extendable, backtrace)

    def __str__(self):
        return f"Interval({self.start}, {self.stop}, extendable={self.extendable}, bt={self.backtrace})"


class IntervalList:
    def __init__(self, intervals: list[Interval]):
        self.intervals = intervals

    @classmethod
    def from_json(cls, json_obj):
        intervals = list(map(Interval.from_json, json_obj['entries']))

        return cls(intervals)


class Layer:
    def __init__(self, k_min: int, k_max: int, diagonals: list[IntervalList]):
        self.k_min = k_min
        self.k_max = k_max
        self.diagonals: list[IntervalList] = diagonals

    @classmethod
    def from_json(cls, json_obj):
        k_min = json_obj['k_min']
        k_max = json_obj['k_max']
        diagonals = list(map(IntervalList.from_json, json_obj['diagonals']))

        return cls(k_min, k_max, diagonals)

    def to_endpoint_df(self):
        data = []
        for k, intervals in zip(range(self.k_min, self.k_max+1), self.diagonals):
            for interval in intervals.intervals:
                if interval.extendable:
                    offset = interval.stop - 1
                    row = offset - k
                    data.append({
                        'diag': int(k),
                        'rank': int(row),
                        'offset': int(offset)
                    })

        return pandas.DataFrame(data, dtype=int)

    def to_arrows(self):
        for k, intervals in zip(range(self.k_min, self.k_max+1), self.diagonals):
            for interval in intervals.intervals:
                start_offset = interval.start
                start_row = start_offset - k

                end_offset = interval.stop
                delta = (end_offset - start_offset)
                print(k, interval, (start_offset, start_row), (end_offset, start_row+delta), file=sys.stderr)

                yield (numpy.array([start_offset, start_row], dtype=float),
                       numpy.array([end_offset, start_row+delta], dtype=float),
                       interval)


class LayerSet:
    def __init__(self):
        self.state = {}

    def items(self):
        return self.state.items()


class PoastaState:
    def __init__(self):
        self.meta = {}
        self.layers: list[tuple[int, LayerSet]] = []

    def last(self) -> LayerSet:
        return self.layers[-1][1]

    def max_score(self) -> int:
        return self.layers[-1][0]

    def new(self, score: int):
        self.layers.append((score, LayerSet()))


def parse_poasta_log(fname) -> PoastaState:
    poasta_state = PoastaState()
    poasta_state.new(0)

    with open(fname) as f:
        state_vars = {}
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith('#'):
                var_name, value = line[2:].split(':', maxsplit=1)
                var_name = var_name.strip()
                value = value.strip()

                state_vars[var_name] = value

                if var_name == 'score':
                    poasta_state.new(int(value))
            else:
                json_obj = json.loads(line)
                layer = Layer.from_json(json_obj)
                poasta_state.last().state[state_vars['state']] = layer

        for k, v in state_vars.items():
            if k not in {"score", "state"}:
                poasta_state.meta[k] = v

    return poasta_state


def main():
    parser = argparse.ArgumentParser(description="Plot POASTA aligner computation state")

    parser.add_argument('graph', type=Path,
                        help="The graph used for alignment in DOT format")
    parser.add_argument('poasta_log', type=Path,
                        help="The POASTA debug log for a specific sequence")
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help="Output directory")
    parser.add_argument('-w', '--fig-width', type=int, default=None, required=False, help="Figure width in inches")
    parser.add_argument('-t', '--truth', type=Path, required=False, default=None, help="Path to SPOA truth matrix")

    args = parser.parse_args()

    if args.output is None:
        print("ERROR: no output directory specified!", file=sys.stderr)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)

    print("Loading graph...", file=sys.stderr)
    g, tsorted = load_graph(args.graph)
    graph_layout = poa_graph_layout(g)

    max_score = None
    true_matrix = None
    xlabels = None
    if args.truth:
        print("Loading SPOA truth...", file=sys.stderr)
        xlabels, ylabels, true_matrix = load_spoa_matrix(args.truth)
        max_score = int(true_matrix['score'].max())
        true_matrix = true_matrix.pivot(index="rank", columns="offset", values="score")

    print("Loading poasta data...", file=sys.stderr)
    poasta_state = parse_poasta_log(args.poasta_log)

    print("Creating plots...", file=sys.stderr)
    if max_score is None:
        max_score = poasta_state.max_score()

    max_rank = int(poasta_state.meta['num_nodes'])
    max_offset = int(poasta_state.meta['length'])

    fig_width = args.fig_width if args.fig_width else int(round(len(poasta_state.meta.get('seq', "-" * 60)) / 4))
    fig_height = int(round((fig_width * max_rank) / max_offset))

    if seq := poasta_state.meta.get('seq'):
        xlabels = [f"{i}\n{c}\n" for i, c in zip(range(max_offset+1), f"-{seq}")]
    elif xlabels is None:
        xlabels = list(range(max_offset + 1))

    yticks = numpy.arange(max_rank) + 0.5

    ylabels = [f"{g.nodes[n]['rank']} ({g.nodes[n]['symbol']})" for n in tsorted]
    hlines = poa_matrix_discontinuieties(g, tsorted)

    end_point_dfs = []
    scores = []
    for score, layer_set in poasta_state.layers:
        for aln_state, layer in layer_set.items():
            df = layer.to_endpoint_df()
            df['state'] = aln_state
            end_point_dfs.append(df)
            scores.append(score)

    all_endpoints_df = (pandas.concat(end_point_dfs, keys=scores, names=['score'])
                        .reset_index(level=1, drop=True)
                        .reset_index()
                        .set_index('state')
                        .astype(int))

    if true_matrix is not None:
        true_fig, true_axes = make_dp_matrix_plot(g, graph_layout, hlines, true_matrix, max_score,
                                                  xlabels, ylabels, yticks, (fig_width, fig_height))

        true_axes[1].set_title("True matrix")
        true_fig.savefig(args.output / f"true_matrix.png", dpi=300)

    for score, layer_set in poasta_state.layers:
        for aln_state, layer in layer_set.items():
            if aln_state == 'MisMatch':
                end_points_before = (all_endpoints_df.loc['Extended']
                                     .query(f'score < {score}')
                                     .set_index(['rank', 'offset']))

                end_points_this_score = (all_endpoints_df.loc[aln_state]
                                         .query(f'score == {score}')
                                         .set_index(['rank', 'offset']))

                score_df = pandas.concat([end_points_before, end_points_this_score])
            else:
                score_df = (all_endpoints_df.loc[aln_state]
                            .query(f'score <= {score}')
                            .set_index(['rank', 'offset']))

            print("Processing score", score, f"({aln_state})", file=sys.stderr)
            if len(score_df) == 0:
                continue

            # Only keep cells with the lowest score
            duplicated = score_df.index.duplicated(keep=False)
            if numpy.count_nonzero(duplicated) > 0:
                print("DUPLICATED:", file=sys.stderr)
                print(score_df[duplicated], file=sys.stderr)

            deduplicated = score_df[~score_df.index.duplicated(keep='first')].reset_index()
            # deduplicated = score_df.reset_index()

            # Make heatmap DF
            pivot = deduplicated.pivot(index="rank", columns="offset", values="score")

            # Ensure rows/cols exists even when no data exists for those
            row_idx = pandas.Index(numpy.arange(max_rank + 1).astype(int))
            col_idx = pandas.Index(numpy.arange(max_offset + 1).astype(int))
            rows = pivot.index.union(row_idx, sort=True)
            cols = pivot.columns.union(col_idx, sort=True)
            pivot = pivot.reindex(index=rows, columns=cols).sort_index().sort_index(axis='columns')

            fig, axes = make_dp_matrix_plot(g, graph_layout, hlines, pivot, max_score,
                                            xlabels, ylabels, yticks, (fig_width, fig_height))

            axes[1].set_title(f"{aln_state.upper()} matrix")

            for start, end, interval in layer.to_arrows():
                marker = 'o' if interval.extendable else 'D'
                color = {
                    'Insertion': 'purple',
                    'Deletion': 'red',
                    'MisMatch': 'green',
                    'Extended': 'darkgreen',
                    'Start': 'black',
                }[interval.backtrace.prev_state]

                # color = matplotlib.colormaps['plasma'](1 - arrow_norm(score))
                axes[1].plot([start[0], end[0]], [start[1] - 0.25, end[1]], f'-{marker}', markevery=[1],
                             color=color)

            fig.savefig(args.output / f"score{score}.{aln_state}.png", dpi=300)
            plt.close(fig)


def make_dp_matrix_plot(g, graph_layout, hlines, matrix, max_score, xlabels, ylabels, yticks, figsize):
    fig, axes = plt.subplots(1, 3, figsize=figsize, width_ratios=[1, 8.75, 0.25],
                             constrained_layout=True)

    nx.draw(g, pos=graph_layout, ax=axes[0], node_size=75,
            labels={n: ndata['symbol'] for n, ndata in g.nodes(data=True)},
            font_size=6)

    seaborn.heatmap(matrix,
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

    axes[1].set_xticklabels(xlabels, rotation=0)
    axes[1].sharey(axes[0])
    axes[1].set_yticks(yticks)
    axes[1].set_yticklabels(ylabels)
    axes[1].grid(True, which='both')
    axes[1].set_axisbelow(True)
    axes[1].tick_params(left=False, labelleft=False)
    axes[1].set_ylabel("")
    axes[1].set_ylim(0, len(ylabels))
    axes[1].invert_yaxis()

    axes[2].set_ylabel("Score")

    return fig, axes


if __name__ == '__main__':
    r = main()
    sys.exit(int(r) if r is not None else 0)
