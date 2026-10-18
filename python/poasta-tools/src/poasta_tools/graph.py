"""Load and lay out the POA graph from a POASTA debug DOT file."""

from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple

import networkx as nx


class NodeInfo(NamedTuple):
    rank: int
    symbol: str


def load_graph(dot_path: Path) -> tuple[nx.DiGraph, dict[str, NodeInfo], str | None]:
    """Parse a POASTA debug DOT file.

    Returns
    -------
    G
        Directed graph with node attribute ``rank`` (int) and ``symbol`` (str).
    node_info
        Mapping from networkx node id → NodeInfo(rank, symbol).
    seq_str
        Query sequence string extracted from the leading ``// sequence:`` comment,
        or None if absent.
    """
    text = dot_path.read_text()

    # Extract query sequence from first comment line
    seq_str: str | None = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("//"):
            m = re.match(r"//\s*sequence:\s*(\S+)", line)
            if m:
                seq_str = m.group(1)
            break

    G: nx.DiGraph = nx.nx_agraph.read_dot(dot_path)

    node_info: dict[str, NodeInfo] = {}
    for node, attrs in G.nodes(data=True):
        rank_raw = attrs.get("xlabel", attrs.get("rank", None))
        symbol_raw = attrs.get("label", "?")
        rank = int(rank_raw) if rank_raw is not None else -1
        symbol = symbol_raw.strip('"')
        node_info[node] = NodeInfo(rank=rank, symbol=symbol)

    return G, node_info, seq_str


def body_nodes(node_info: dict[str, NodeInfo]) -> list[str]:
    """Return non-sentinel nodes sorted by rank (excludes # and $ sentinels)."""
    return sorted(
        (n for n, info in node_info.items() if info.symbol not in ("#", "$")),
        key=lambda n: node_info[n].rank,
    )


def graphviz_layout(G: nx.DiGraph) -> dict[str, tuple[float, float]]:
    """Return {node: (x, y)} using Graphviz 'dot' layout."""
    return nx.nx_agraph.graphviz_layout(G, prog="dot", args="-Grankdir=TB")
