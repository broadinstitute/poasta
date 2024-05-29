#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path
from typing import BinaryIO
from bisect import bisect_left

import pygraphviz


def parse_seq_meta(seq_str: str):
    seq_name, start_node = seq_str.rsplit(b':', maxsplit=1)

    return seq_name.decode('utf-8'), start_node.decode('ascii')


def parse_poasta_graphviz(file: BinaryIO):
    seq_names = file.readline()
    seq_names = seq_names.replace(b"# seq:\t", b"").split(b'\t')

    seq_meta = {}
    for i, seq in enumerate(seq_names):
        seq_name, start_node = parse_seq_meta(seq)
        seq_meta[seq_name] = (i, start_node)

    dot_graph = file.read().decode('ascii')

    return pygraphviz.AGraph(dot_graph), seq_meta


def contains(list_to_search: list, value) -> bool:
    i = bisect_left(list_to_search, value)

    if i != len(list_to_search) and list_to_search[i] == value:
        return True

    return False


def main():
    parser = argparse.ArgumentParser(description="Create visualizations for a subgraph of the POA graph.")

    parser.add_argument(
        'poasta_dot', type=Path,
        help="Path to a POASTA graph in DOT format."
    )

    parser.add_argument(
        'region',
        help="Subgraph to extract and visualize. Format: seq_name:start-stop."
    )

    parser.add_argument(
        '-p', '--pos-offset', type=int, default=1,
        help="Base position of sequence in the graph"
    )

    parser.add_argument(
        '-H', '--highlight', type=str, action="append", default=None,
        help="Highlight the path of a specific sequence. Format: seq_name:color. "
             "Option can be used multiple times."
    )

    args = parser.parse_args()

    print("Reading graph...", file=sys.stderr)
    with args.poasta_dot.open('rb') as ifile:
        g, seq_meta = parse_poasta_graphviz(ifile)

    try:
        seq_name, region = args.region.rsplit(':', maxsplit=1)
        start, stop = region.split('-', maxsplit=1)
        start = max(0, int(start) - 1)
        stop = int(stop)

        if seq_name not in seq_meta:
            print("Invalid sequence name '{}'".format(seq_name), file=sys.stderr)
            sys.exit(1)

        print("Extracting subgraph belonging to", args.region, file=sys.stderr)
    except (IndexError, ValueError):
        print("Invalid region '{}'".format(args.region), file=sys.stderr)
        print("Region should be in the format of seq_name:start-stop", file=sys.stderr)
        sys.exit(1)

    highlight_seq_ids = {}
    if args.highlight:
        for ref_highlight in args.highlight:
            ref_name, color = ref_highlight.rsplit(':', maxsplit=1)
            if ref_name not in seq_meta:
                print(f"Invalid sequence '{ref_name}', ignoring highlight.", file=sys.stderr)
                continue

            highlight_seq_ids[seq_meta[ref_name][0]] = color

    aligned_nodes_subgraphs = {}
    for subgraph in g.subgraphs():
        for n in subgraph.nodes_iter():
            aligned_nodes_subgraphs[n] = list(subgraph.nodes_iter())

    seq_id, start_node = seq_meta[seq_name]

    curr_pos = max(0, args.pos_offset - 1)
    curr_node = start_node
    seq_path = []
    while curr_node:
        if curr_pos >= start and curr_pos < stop:
            seq_path.append(curr_node)
        elif curr_pos >= stop:
            break

        next_node = None
        for out_edge in g.out_edges_iter(curr_node):
            seq_ids = [int(v[1:]) for v in out_edge.attr['class'].split()]
            if contains(seq_ids, seq_id):
                next_node = out_edge[1]

        curr_node = next_node
        curr_pos += 1

    new_subgraph = pygraphviz.AGraph(directed=True, rankdir='LR', nodesep=0.05, pad=0.025, fontname="Arial")
    new_subgraph.node_attr['shape'] = "square"
    new_subgraph.node_attr['style'] = "filled"
    new_subgraph.node_attr['fillcolor'] = "#e3e3e3"
    new_subgraph.node_attr['penwidth'] = 0
    new_subgraph.node_attr['margin'] = 0.05

    for i, n in enumerate(seq_path):
        if n in aligned_nodes_subgraphs:
            for aln_node in aligned_nodes_subgraphs[n]:
                new_subgraph.add_node(aln_node, **g.get_node(aln_node).attr)

            new_subgraph.add_subgraph(aligned_nodes_subgraphs[n], rank="same")
        else:
            new_subgraph.add_node(n, **g.get_node(n).attr)

        new_subgraph.get_node(n).attr['xlabel'] = f"<<font color='black'>{start+i+1}</font>>"

    all_region_nodes = set(new_subgraph.nodes_iter())
    for e in g.edges_iter(all_region_nodes):
        if e[0] not in all_region_nodes or e[1] not in all_region_nodes:
            continue

        seq_ids = set(int(v[1:]) for v in e.attr['class'].split())
        edge_color_list = []
        for seq_id, color in highlight_seq_ids.items():
            if seq_id in seq_ids:
                edge_color_list.append(color)

        edge_attr = dict(**e.attr)
        if edge_color_list:
            if len(seq_ids) - len(edge_color_list) > 0:
                edge_color_list = ["black", *edge_color_list]

            edge_attr['color'] = ":".join(edge_color_list)

        new_subgraph.add_edge(e[0], e[1], **edge_attr)

    print(new_subgraph.to_string())


if __name__ == '__main__':
    main()
