"""Load band and cell TSV files produced by POASTA's --debug-output-dir."""

from pathlib import Path

import numpy as np
import polars as pl

INF_THRESHOLD = 1 << 28  # matches Rust const INF = 1 << 28


def load_bands(path: Path) -> pl.DataFrame:
    """Return DataFrame with columns: k, node_rank, qlo, qhi."""
    return pl.read_csv(path, separator="\t")


def load_cells(path: Path) -> pl.DataFrame:
    """Return DataFrame with columns: k, node_rank, q, <state>...

    State columns (e.g. Match, Deletion, Insertion) are discovered from the header.
    Values >= INF_THRESHOLD are replaced with NaN so they are masked in plots.
    """
    df = pl.read_csv(path, separator="\t")
    state_cols = [c for c in df.columns if c not in ("k", "node_rank", "q")]
    # Cast state cols to float and mask INF values
    df = df.with_columns([
        pl.when(pl.col(c) >= INF_THRESHOLD)
        .then(None)
        .otherwise(pl.col(c))
        .cast(pl.Float64)
        .alias(c)
        for c in state_cols
    ])
    return df


def state_columns(cells_df: pl.DataFrame) -> list[str]:
    return [c for c in cells_df.columns if c not in ("k", "node_rank", "q")]


def cell_band_map(bands_df: pl.DataFrame) -> dict[tuple[int, int], int]:
    """Map (node_rank, q) → smallest k whose band covers q for that node.

    A band covers q when qlo <= q < qhi.
    """
    mapping: dict[tuple[int, int], int] = {}
    for row in bands_df.sort("k").iter_rows(named=True):
        k, rank, qlo, qhi = row["k"], row["node_rank"], row["qlo"], row["qhi"]
        for q in range(qlo, qhi):
            key = (rank, q)
            if key not in mapping:
                mapping[key] = k
    return mapping


def pivot_state(
    cells_df: pl.DataFrame,
    state: str,
    all_ranks: list[int],
    all_q: list[int],
) -> np.ndarray:
    """Return a 2-D array [rank_idx, q_idx] for the given state column.

    Rows correspond to `all_ranks` (ascending), columns to `all_q` (ascending).
    Missing cells are NaN.
    """
    rank_idx = {r: i for i, r in enumerate(all_ranks)}
    q_idx = {q: i for i, q in enumerate(all_q)}
    matrix = np.full((len(all_ranks), len(all_q)), np.nan)

    for row in cells_df.select(["node_rank", "q", state]).iter_rows():
        rank, q, val = row
        if val is not None and rank in rank_idx and q in q_idx:
            matrix[rank_idx[rank], q_idx[q]] = val

    return matrix


def band_color_matrix(
    cells_df: pl.DataFrame,
    bands_df: pl.DataFrame,
    all_ranks: list[int],
    all_q: list[int],
) -> np.ndarray:
    """Return a 2-D float array [rank_idx, q_idx] where each cell value is the k
    at which it was first computed. NaN for cells not in any band."""
    bmap = cell_band_map(bands_df)
    rank_idx = {r: i for i, r in enumerate(all_ranks)}
    q_idx = {q: i for i, q in enumerate(all_q)}

    # Only include (rank, q) pairs that actually appear in cells_df
    present = set(
        cells_df.select(["node_rank", "q"]).iter_rows()
    )

    matrix = np.full((len(all_ranks), len(all_q)), np.nan)
    for (rank, q), k in bmap.items():
        if (rank, q) in present and rank in rank_idx and q in q_idx:
            matrix[rank_idx[rank], q_idx[q]] = float(k)
    return matrix
