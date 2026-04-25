"""poasta-plot — visualise POASTA band-doubling debug output."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from .data import load_bands, load_cells, state_columns
from .graph import load_graph
from .plot import ColorBy, make_figure, save_animation_frame, save_static


def _discover_sequences(debug_dir: Path) -> list[str]:
    """Return stems of all *_bands.tsv found in debug_dir."""
    seqs = [p.stem.removesuffix("_bands")
            for p in sorted(debug_dir.glob("*_bands.tsv"))]
    if not seqs:
        raise click.UsageError(
            f"No *_bands.tsv files found in {debug_dir}. "
            "Run poasta with --debug-output-dir to generate them."
        )
    return seqs


def _process_one(
    debug_dir: Path,
    sequence: str | None,
    output_dir: Path,
    animate: bool,
    color_by: ColorBy,
) -> None:
    seqs = [sequence] if sequence else _discover_sequences(debug_dir)
    for seq in seqs:
        click.echo(f"\n  Sequence: {seq}")
        _process_sequence(debug_dir, seq, output_dir, animate, color_by)


def _process_sequence(
    debug_dir: Path,
    seq: str,
    output_dir: Path,
    animate: bool,
    color_by: ColorBy,
) -> None:
    dot_path = debug_dir / f"graph_for_{seq}.dot"
    bands_path = debug_dir / f"{seq}_bands.tsv"
    cells_path = debug_dir / f"{seq}_cells.tsv"

    for p in (dot_path, bands_path, cells_path):
        if not p.exists():
            raise click.UsageError(f"Expected file not found: {p}")

    click.echo(f"  Loading graph from {dot_path.name} …")
    G, node_info, seq_str = load_graph(dot_path)

    click.echo(f"  Loading bands from {bands_path.name} …")
    try:
        bands_df = load_bands(bands_path)
    except Exception as exc:
        click.echo(f"  WARNING: skipping {seq}: could not read {bands_path.name}: {exc}", err=True)
        return

    click.echo(f"  Loading cells from {cells_path.name} …")
    try:
        cells_df = load_cells(cells_path)
    except Exception as exc:
        click.echo(f"  WARNING: skipping {seq}: could not read {cells_path.name}: {exc}", err=True)
        return

    if bands_df.height == 0 or cells_df.height == 0:
        click.echo(f"  WARNING: skipping {seq}: empty bands or cells TSV", err=True)
        return

    states = state_columns(cells_df)
    click.echo(f"  States: {states}")

    if animate:
        k_values = sorted(bands_df["k"].unique().to_list())
        click.echo(f"  Generating animation frames for k in {k_values} …")
        for k in k_values:
            for state in states:
                fig = make_figure(cells_df, bands_df, G, node_info, seq_str,
                                  state, color_by, k_max=k)
                out = save_animation_frame(fig, output_dir, seq, state, k)
                click.echo(f"    {out}")
        anim_dir = output_dir / "animation"
        click.echo(
            f"\n  Frames written to {anim_dir}/\n"
            f"  Combine with ffmpeg, e.g.:\n"
            f"    ffmpeg -framerate 4 -i '{anim_dir}/{seq}_Match_k%04d.png' output.mp4"
        )
    else:
        for state in states:
            fig = make_figure(cells_df, bands_df, G, node_info, seq_str,
                              state, color_by)
            out = save_static(fig, output_dir, seq, state)
            click.echo(f"  Saved {out}")


@click.command()
@click.argument("debug_dirs", nargs=-1, required=True,
                type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("--sequence", "-s", default=None, metavar="NAME",
              help="Sequence name (default: auto-detect from *_bands.tsv).")
@click.option("--output", "-o", required=True,
              type=click.Path(file_okay=False, path_type=Path),
              help="Output directory for plots.")
@click.option("--animate", is_flag=True, default=False,
              help="Generate one PNG per k-iteration per state.")
@click.option("--color-by", "color_by",
              type=click.Choice(["band", "value"], case_sensitive=False),
              default="band", show_default=True,
              help="Color cells by band (categorical) or by DP value (continuous).")
def main(
    debug_dirs: tuple[Path, ...],
    sequence: str | None,
    output: Path,
    animate: bool,
    color_by: str,
) -> None:
    """Visualise POASTA band-doubling debug output.

    DEBUG_DIRS are one or more directories written by 'poasta align
    --debug-output-dir'. Each directory is processed independently.
    When multiple directories are given, plots are written to
    OUTPUT/<dir-name>/ subdirectories.
    """
    color_by_typed: ColorBy = color_by  # type: ignore[assignment]
    multi = len(debug_dirs) > 1

    for debug_dir in debug_dirs:
        out_dir = output / debug_dir.name if multi else output
        click.echo(f"\nProcessing {debug_dir} → {out_dir}")
        try:
            _process_one(debug_dir, sequence, out_dir, animate, color_by_typed)
        except click.UsageError:
            raise
        except Exception as exc:
            click.echo(f"  ERROR: {exc}", err=True)
            sys.exit(1)
