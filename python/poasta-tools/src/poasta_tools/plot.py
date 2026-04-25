"""Generate DP matrix plots from POASTA debug output."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import seaborn as sns

from .data import (
    band_color_matrix,
    pivot_state,
    state_columns,
)
from .graph import NodeInfo, body_nodes, graphviz_layout

import polars as pl

ColorBy = Literal["band", "value"]

# Width ratio: graph panel vs heatmap panel
_GRAPH_PANEL_WEIGHT = 1
_HEAT_PANEL_WEIGHT = 6


def _discontinuity_ranks(
    G: nx.DiGraph, node_info: dict[str, NodeInfo]
) -> list[int]:
    """Ranks where topological order breaks — successors of branch points or
    merge nodes themselves — so a horizontal line can separate non-adjacent rows."""
    hlines: set[int] = set()
    for n in nx.topological_sort(G):
        if n not in node_info:
            continue
        if G.out_degree(n) > 1:
            for neighbor in G.successors(n):
                if neighbor in node_info:
                    hlines.add(node_info[neighbor].rank)
        elif G.in_degree(n) > 1 and G.out_degree(n) > 0:
            hlines.add(node_info[n].rank)
    return sorted(hlines)


def _rank_label(node_info: dict[str, NodeInfo], node: str) -> str:
    info = node_info[node]
    return f"{info.rank}:{info.symbol}"


def _draw_graph_panel(
    ax: plt.Axes,
    G: nx.DiGraph,
    node_info: dict[str, NodeInfo],
    body: list[str],
    pos: dict[str, tuple[float, float]],
) -> None:
    """Draw the mini graph side panel.

    X positions come from Graphviz layout; Y positions are node ranks so they
    align with the heatmap rows.
    """
    rank_of = {n: node_info[n].rank for n in body}
    body_set = set(body)

    # Map graphviz x to panel x, keep rank as y
    gviz_xs = {n: pos[n][0] for n in body if n in pos}
    if not gviz_xs:
        ax.set_visible(False)
        return

    x_min, x_max = min(gviz_xs.values()), max(gviz_xs.values())
    x_range = (x_max - x_min) or 1.0

    def norm_x(n: str) -> float:
        return (gviz_xs.get(n, x_min) - x_min) / x_range

    # Draw edges
    for src, dst in G.edges():
        if src not in body_set or dst not in body_set:
            continue
        xs = [norm_x(src), norm_x(dst)]
        ys = [rank_of[src], rank_of[dst]]
        ax.plot(xs, ys, color="gray", linewidth=0.6, zorder=1)

    # Draw nodes
    xs_nodes = [norm_x(n) for n in body]
    ys_nodes = [rank_of[n] for n in body]
    ax.scatter(xs_nodes, ys_nodes, s=40, color="steelblue", zorder=2)
    for n, xn, yn in zip(body, xs_nodes, ys_nodes):
        ax.text(xn, yn, node_info[n].symbol, fontsize=6, ha="center", va="center",
                color="white", zorder=3)

    all_ranks = [rank_of[n] for n in body]
    ax.set_ylim(min(all_ranks) - 0.5, max(all_ranks) + 0.5)
    ax.invert_yaxis()
    ax.set_xlim(-0.1, 1.1)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")


def make_figure(
    cells_df: pl.DataFrame,
    bands_df: pl.DataFrame,
    G: nx.DiGraph,
    node_info: dict[str, NodeInfo],
    seq_str: str | None,
    state: str,
    color_by: ColorBy,
    k_max: int | None = None,
) -> plt.Figure:
    """Build one figure for a single DP state matrix.

    Parameters
    ----------
    k_max
        When set, only cells/bands with k <= k_max are shown (animation mode).
    """
    if k_max is not None:
        cells_df = cells_df.filter(pl.col("k") <= k_max)

    body = body_nodes(node_info)
    if not body:
        raise ValueError("Graph has no non-sentinel nodes")

    all_ranks = sorted({node_info[n].rank for n in body})
    rank_to_nodes: dict[int, list[str]] = {}
    for n in body:
        rank_to_nodes.setdefault(node_info[n].rank, []).append(n)

    q_vals = sorted(cells_df["q"].unique().to_list())
    if not q_vals:
        raise ValueError(f"No cells found for state '{state}' (k_max={k_max})")

    # ── Build color matrix ──────────────────────────────────────────────────
    if color_by == "band":
        bands_filtered = bands_df if k_max is None else bands_df.filter(pl.col("k") <= k_max)
        cmatrix = band_color_matrix(cells_df, bands_filtered, all_ranks, q_vals)
        unique_k = sorted({int(v) for v in cmatrix[~np.isnan(cmatrix)]})
        palette = sns.color_palette("tab20", max(len(unique_k), 1))
        k_to_color = {k: palette[i % len(palette)] for i, k in enumerate(unique_k)}
        rgba = np.zeros((*cmatrix.shape, 4))
        for ri in range(cmatrix.shape[0]):
            for ci in range(cmatrix.shape[1]):
                v = cmatrix[ri, ci]
                if np.isnan(v):
                    rgba[ri, ci] = (1, 1, 1, 0)
                else:
                    r, g, b = k_to_color[int(v)]
                    rgba[ri, ci] = (r, g, b, 1)
        cmap_img = rgba
        legend_handles = [
            mpatches.Patch(color=k_to_color[k], label=f"k={k}") for k in unique_k
        ]
    else:
        vmatrix = pivot_state(cells_df, state, all_ranks, q_vals)
        valid = vmatrix[~np.isnan(vmatrix)]
        vmin, vmax = (float(valid.min()), float(valid.max())) if valid.size else (0, 1)
        normed = (vmatrix - vmin) / max(vmax - vmin, 1)
        cmap = plt.colormaps["turbo"]
        rgba = cmap(normed)
        rgba[np.isnan(vmatrix)] = (1, 1, 1, 0)
        cmap_img = rgba
        legend_handles = []

    # ── Layout graph ────────────────────────────────────────────────────────
    try:
        pos = graphviz_layout(G)
    except Exception:
        pos = {}

    # ── Figure layout ───────────────────────────────────────────────────────
    n_rows, n_cols = len(all_ranks), len(q_vals)
    cell_px = 8
    fig_w = max(6.0, (n_cols * cell_px + 180) / 72)
    fig_h = max(4.0, (n_rows * cell_px + 120) / 72)

    fig, (ax_graph, ax_heat) = plt.subplots(
        1, 2,
        figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [_GRAPH_PANEL_WEIGHT, _HEAT_PANEL_WEIGHT]},
    )
    fig.suptitle(f"DP matrix — {state}" + (f"  (k ≤ {k_max})" if k_max is not None else ""),
                 fontsize=9)

    # ── Graph side panel ────────────────────────────────────────────────────
    _draw_graph_panel(ax_graph, G, node_info, body, pos)
    ax_graph.set_yticks(all_ranks)
    ax_graph.set_yticklabels([f"{r}" for r in all_ranks], fontsize=6)
    ax_graph.yaxis.set_visible(True)
    ax_graph.set_ylim(min(all_ranks) - 0.5, max(all_ranks) + 0.5)
    ax_graph.invert_yaxis()

    # ── Heatmap ─────────────────────────────────────────────────────────────
    ax_heat.imshow(
        cmap_img,
        aspect="auto",
        origin="upper",
        extent=[-0.5, n_cols - 0.5, max(all_ranks) + 0.5, min(all_ranks) - 0.5],
        interpolation="nearest",
    )

    # Y-axis: rank labels (rank:symbol for first node at that rank)
    y_labels = []
    for r in all_ranks:
        nodes_at_rank = rank_to_nodes.get(r, [])
        sym = node_info[nodes_at_rank[0]].symbol if nodes_at_rank else "?"
        y_labels.append(f"{r}:{sym}")
    ax_heat.set_yticks(all_ranks)
    ax_heat.set_yticklabels(y_labels, fontsize=6)

    # X-axis: query positions with optional sequence characters
    x_ticks = list(range(n_cols))
    def _q_label(q: int) -> str:
        if q == 0:
            return "-"
        if seq_str and 0 < q <= len(seq_str):
            return seq_str[q - 1]
        return str(q)

    x_labels = [_q_label(q) for q in q_vals]

    step = max(1, n_cols // 40)
    shown = x_ticks[::step]
    ax_heat.set_xticks(shown)
    ax_heat.set_xticklabels([x_labels[i] for i in shown], fontsize=6, rotation=90)
    ax_heat.set_xlabel("Query position", fontsize=7)
    ax_heat.set_ylabel("")

    # Overlay cell scores
    score_matrix = pivot_state(cells_df, state, all_ranks, q_vals)
    for ri, r in enumerate(all_ranks):
        for ci in range(n_cols):
            v = score_matrix[ri, ci]
            if np.isnan(v):
                continue
            ax_heat.text(
                ci, r, f"{int(v)}",
                fontsize=4, ha="center", va="center",
                color="black", zorder=5,
            )

    # Draw graph discontinuity lines
    for r in _discontinuity_ranks(G, node_info):
        ax_heat.axhline(r - 0.5, color="black", linewidth=0.6, zorder=4)

    # Draw grid lines at band boundaries
    ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_heat.set_yticks(np.arange(min(all_ranks) - 0.5, max(all_ranks) + 1, 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.3, alpha=0.4)
    ax_heat.tick_params(which="minor", bottom=False, left=False)

    if legend_handles:
        ax_heat.legend(
            handles=legend_handles,
            fontsize=5,
            loc="upper right",
            ncol=max(1, len(legend_handles) // 10),
            framealpha=0.7,
        )

    # Share y-axis range
    ax_graph.set_ylim(ax_heat.get_ylim())

    plt.tight_layout()
    return fig


def save_static(fig: plt.Figure, output_dir: Path, seq_name: str, state: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"{seq_name}.{state}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def save_animation_frame(
    fig: plt.Figure, output_dir: Path, seq_name: str, state: str, k: int
) -> Path:
    anim_dir = output_dir / "animation"
    anim_dir.mkdir(parents=True, exist_ok=True)
    out = anim_dir / f"{seq_name}.k{k:04d}.{state}.png"
    fig.savefig(out, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return out
