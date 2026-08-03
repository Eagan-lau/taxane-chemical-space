#!/usr/bin/env python3
"""Merge the frozen V5 topology and bridge-hypothesis figures.

This is a presentation-only workflow. It reads frozen source tables, renders
one integrated Figure 4, and writes provenance and numerical audit files. It
does not rebuild reaction rules, enumerate chemical space, or discover new
bridges.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import shutil
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image, ImageChops
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


DEFAULT_SOURCE = Path("source_data")
DEFAULT_OUTPUT = Path("figure4_output")
FIGURE_STEM = "Figure_4_Bridge_Hypotheses_to_Global_Topology_V8"

COLORS = {
    "ink": "#20262e",
    "gray": "#7d858d",
    "light_gray": "#e4e7eb",
    "blue": "#2f78bd",
    "blue_light": "#a8c9e5",
    "teal": "#2a8c7b",
    "orange": "#d87924",
    "orange_light": "#f6dfc8",
    "gold": "#d7a321",
    "red": "#b94b52",
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.4,
            "axes.titlesize": 9.0,
            "axes.labelsize": 7.7,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.4,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def wrap(value: object, width: int) -> str:
    return "\n".join(textwrap.wrap(str(value), width=width))


def short_name(value: object) -> str:
    text = str(value)
    return text.split(" (", 1)[0]


def formula_label(value: object) -> str:
    tokens = []
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", str(value)):
        suffix = rf"_{{{count}}}" if count and count != "1" else ""
        tokens.append(rf"\mathrm{{{element}}}{suffix}")
    return "$" + "".join(tokens) + "$" if tokens else str(value)


def quiet_grid(ax, axis: str = "y") -> None:
    ax.grid(axis=axis, color="#dfe3e7", linewidth=0.55, alpha=0.85)
    ax.set_axisbelow(True)


def panel_label(ax, label: str, x: float = -0.13, y: float = 1.06) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def draw_molecule(ax, smiles: str, *, intermediate: bool = False) -> None:
    ax.set_axis_off()
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        ax.text(
            0.5,
            0.5,
            "structure unavailable",
            ha="center",
            va="center",
            color=COLORS["gray"],
        )
        return

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol, canonOrient=True, clearConfs=True)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    drawer = rdMolDraw2D.MolDraw2DCairo(1000, 700)
    options = drawer.drawOptions()
    options.useBWAtomPalette()
    options.padding = 0.045
    options.bondLineWidth = 1.8
    options.minFontSize = 18
    options.maxFontSize = 30
    options.fixedBondLength = 29
    options.addStereoAnnotation = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    image = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")
    background = Image.new("RGB", image.size, "white")
    bbox = ImageChops.difference(image, background).getbbox()
    if bbox:
        margin = 14
        bbox = (
            max(0, bbox[0] - margin),
            max(0, bbox[1] - margin),
            min(image.width, bbox[2] + margin),
            min(image.height, bbox[3] + margin),
        )
        image = image.crop(bbox)
    ax.imshow(image, interpolation="lanczos", aspect="equal")
    ax.set_xlim(0, image.width)
    ax.set_ylim(image.height, 0)

    if intermediate:
        ax.add_patch(
            Rectangle(
                (0.01, 0.01),
                0.98,
                0.98,
                transform=ax.transAxes,
                facecolor="none",
                edgecolor=COLORS["orange"],
                linewidth=1.15,
                zorder=10,
            )
        )


def display_reaction_type(value: object) -> str:
    text = str(value)
    labels = {
        "deacetylation_or_acetyl_ester_hydrolysis": (
            "deacetylation or\nacetyl-ester hydrolysis"
        ),
        "acetylation_or_deacetylation_like_acyl_transfer": (
            "acetylation-like or\n"
            "deacetylation-like acyl transfer"
        ),
    }
    return labels.get(text, wrap(text.replace("_", " "), 24))


def compact_reaction_type(value: object) -> str:
    text = str(value)
    labels = {
        "deacetylation_or_acetyl_ester_hydrolysis": (
            "deacetylation / acetyl-ester hydrolysis"
        ),
        "acetylation_or_deacetylation_like_acyl_transfer": (
            "acetylation-like / deacetylation-like acyl transfer"
        ),
    }
    return labels.get(text, text.replace("_", " "))


def add_route_panel(fig, spec, row: pd.Series, label: str) -> None:
    subgrid = spec.subgridspec(
        3,
        3,
        height_ratios=[0.13, 0.65, 0.22],
        hspace=0.01,
        wspace=0.08,
    )
    title_ax = fig.add_subplot(subgrid[0, :])
    title_ax.set_axis_off()
    title_ax.text(
        -0.02,
        0.76,
        label,
        transform=title_ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="center",
    )
    title_ax.text(
        0.04,
        0.76,
        f"{short_name(row['source_name'])}  to  {short_name(row['target_name'])}",
        transform=title_ax.transAxes,
        fontsize=7.3,
        fontweight="bold",
        color=COLORS["ink"],
        va="center",
    )
    title_ax.plot(
        [0.04, 1.0],
        [0.08, 0.08],
        transform=title_ax.transAxes,
        color=COLORS["light_gray"],
        linewidth=0.7,
        clip_on=False,
    )

    axes = [fig.add_subplot(subgrid[1, column]) for column in range(3)]
    molecules = [
        row["source_smiles"],
        row["bridge_smiles"],
        row["target_smiles"],
    ]
    for index, (ax, smiles) in enumerate(zip(axes, molecules)):
        draw_molecule(ax, smiles, intermediate=index == 1)

    axes[0].set_title(
        f"{row['source_space_id']}\n{wrap(short_name(row['source_name']), 25)}",
        fontsize=6.1,
        color=COLORS["ink"],
        pad=2,
    )
    axes[1].set_title(
        f"generated {row['bridge_space_id']}\n{formula_label(row['bridge_formula'])}",
        fontsize=6.1,
        color=COLORS["orange"],
        pad=2,
    )
    axes[2].set_title(
        f"{row['target_space_id']}\n{wrap(short_name(row['target_name']), 25)}",
        fontsize=6.1,
        color=COLORS["ink"],
        pad=2,
    )

    for left in axes[:2]:
        left.add_patch(
            FancyArrowPatch(
                (0.94, 0.52),
                (1.10, 0.52),
                transform=left.transAxes,
                arrowstyle="-|>",
                mutation_scale=9,
                linewidth=1.0,
                color=COLORS["blue"],
                clip_on=False,
                zorder=12,
            )
        )

    metadata = fig.add_subplot(subgrid[2, :])
    metadata.set_axis_off()
    evidence = str(row["evidence_label"]).replace("_", " ")
    metadata.text(
        0.02,
        0.82,
        f"step 1: {compact_reaction_type(row['incoming_reaction_type'])}"
        f"   |   step 2: {compact_reaction_type(row['outgoing_reaction_type'])}",
        transform=metadata.transAxes,
        fontsize=5.25,
        color=COLORS["blue"],
    )
    metadata.text(
        0.02,
        0.50,
        f"{row['evidence_tier']}  |  {evidence}",
        transform=metadata.transAxes,
        fontsize=6.3,
        fontweight="bold",
        color=COLORS["blue"],
    )
    metadata.text(
        0.02,
        0.14,
        f"directed-pair support {int(row['directed_pair_support'])}  |  "
        f"parents {int(row['first_observation_parent_count'])}  |  "
        f"deduplicated paths {int(float(row['deduplicated_structural_path_count']))}  |  "
        f"minimum pairwise Tanimoto {float(row['min_pairwise_tanimoto']):.3f}",
        transform=metadata.transAxes,
        fontsize=5.55,
        color=COLORS["gray"],
    )


def render_figure(
    bridge_table: pd.DataFrame,
    convergence: pd.DataFrame,
    topology: pd.DataFrame,
    routes: pd.DataFrame,
):
    generation_summary = bridge_table[
        bridge_table["record_type"].eq("generation_summary")
    ].copy()
    generation_summary["generation"] = pd.to_numeric(
        generation_summary["generation"]
    )
    directed_pairs = bridge_table[
        bridge_table["record_type"].eq("directed_pair")
    ].copy()
    fractions = topology[
        topology["record_type"].eq("cycle_fractions")
    ].copy()
    fractions = fractions.sort_values("generation")

    fig = plt.figure(figsize=(14.2, 12.4))
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=[0.385, 0.385, 0.23],
        hspace=0.28,
        left=0.055,
        right=0.985,
        bottom=0.065,
        top=0.935,
    )
    upper = outer[0].subgridspec(1, 2, wspace=0.12)
    middle = outer[1].subgridspec(1, 2, wspace=0.12)
    route_specs = [upper[0, 0], upper[0, 1], middle[0, 0], middle[0, 1]]
    for index, (spec, (_, row)) in enumerate(
        zip(route_specs, routes.iterrows())
    ):
        add_route_panel(fig, spec, row, chr(ord("A") + index))

    bottom = outer[2].subgridspec(1, 3, wspace=0.34)

    # E. Bridge definition with generation-specific count bars.
    ax = fig.add_subplot(bottom[0, 0])
    ax.set_axis_off()
    ax.set_title("Bridge definition and generation-specific support", pad=7)
    panel_label(ax, "E")
    values = generation_summary["latent_bridge_candidates"].astype(float)

    ax.text(
        0.22,
        0.88,
        "Directed two-edge bridge",
        transform=ax.transAxes,
        ha="center",
        fontsize=5.7,
        fontweight="bold",
        color=COLORS["ink"],
    )
    motif_x = [0.07, 0.22, 0.37]
    motif_colors = [COLORS["red"], COLORS["orange"], COLORS["red"]]
    motif_sizes = [230, 310, 230]
    for x_value, color, size in zip(motif_x, motif_colors, motif_sizes):
        ax.scatter(
            [x_value],
            [0.59],
            s=size,
            color=color,
            edgecolor=COLORS["ink"],
            linewidth=0.55,
            transform=ax.transAxes,
            zorder=5,
        )
    for start, end in zip(motif_x[:-1], motif_x[1:]):
        ax.add_patch(
            FancyArrowPatch(
                (start + 0.032, 0.59),
                (end - 0.032, 0.59),
                transform=ax.transAxes,
                arrowstyle="-|>",
                mutation_scale=8,
                linewidth=0.95,
                color=COLORS["blue"],
            )
        )
    for x_value, label in zip(
        motif_x,
        ["known G0 A", "generated\nintermediate", "known G0 B"],
    ):
        ax.text(
            x_value,
            0.39,
            label,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=4.7,
            color=COLORS["gray"],
        )
    ax.text(
        0.22,
        0.16,
        "One intermediate reconnects two known taxanes",
        transform=ax.transAxes,
        ha="center",
        fontsize=4.8,
        color=COLORS["gray"],
    )
    ax.plot(
        [0.45, 0.45],
        [0.10, 0.90],
        transform=ax.transAxes,
        color=COLORS["light_gray"],
        linewidth=0.8,
    )

    ax.text(
        0.72,
        0.88,
        "Bridge candidates",
        transform=ax.transAxes,
        ha="center",
        fontsize=5.7,
        fontweight="bold",
        color=COLORS["ink"],
    )
    bar_x = 0.57
    bar_width = 0.32
    bar_height = 0.075
    for y_value, generation, value, color in [
        (0.70, "G1", values.iloc[0], COLORS["blue"]),
        (0.51, "G2", values.iloc[1], COLORS["teal"]),
    ]:
        width = bar_width * float(value) / float(values.max())
        ax.text(
            0.49,
            y_value,
            generation,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=5.7,
            fontweight="bold",
        )
        ax.add_patch(
            Rectangle(
                (bar_x, y_value - bar_height / 2),
                width,
                bar_height,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor=COLORS["ink"],
                linewidth=0.45,
            )
        )
        ax.text(
            bar_x + width + 0.018,
            y_value,
            f"{int(value)}",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=6.2,
            fontweight="bold",
            color=COLORS["ink"],
        )
    ax.text(
        0.49,
        0.32,
        "G3",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=5.7,
        fontweight="bold",
    )
    ax.plot(
        [bar_x, bar_x + bar_width],
        [0.32, 0.32],
        transform=ax.transAxes,
        color=COLORS["gold"],
        linewidth=1.0,
        linestyle=(0, (3, 2)),
    )
    ax.text(
        bar_x + 0.01,
        0.355,
        "right-censored",
        transform=ax.transAxes,
        ha="left",
        fontsize=4.8,
        color=COLORS["gold"],
    )
    ax.text(
        0.49,
        0.16,
        f"{int(values.sum())} bridge intermediates",
        transform=ax.transAxes,
        ha="left",
        fontsize=6.2,
        fontweight="bold",
        color=COLORS["orange"],
    )
    ax.text(
        0.49,
        0.07,
        f"{len(directed_pairs)} directed G0 pairs",
        transform=ax.transAxes,
        ha="left",
        fontsize=6.2,
        fontweight="bold",
        color=COLORS["red"],
    )

    # F. Convergent structures as 100-dot displays.
    ax = fig.add_subplot(bottom[0, 1])
    ax.set_axis_off()
    ax.set_title("Expansion becomes overwhelmingly convergent", pad=7)
    panel_label(ax, "F")
    convergence = convergence.sort_values("generation")
    fraction_pct = 100 * convergence["convergent_fraction"].astype(float)
    group_starts = [0, 12, 24]
    group_colors = [COLORS["blue"], COLORS["teal"], COLORS["gold"]]
    for start, color, pct, row in zip(
        group_starts,
        group_colors,
        fraction_pct,
        convergence.itertuples(),
    ):
        x_coords = np.tile(np.arange(10), 10) + start
        y_coords = np.repeat(np.arange(10), 10)
        display_count = int(np.floor(float(pct) + 0.5))
        ax.scatter(
            x_coords,
            y_coords,
            s=13,
            color="#e8ecef",
            edgecolor="none",
        )
        ax.scatter(
            x_coords[:display_count],
            y_coords[:display_count],
            s=13,
            color=color,
            edgecolor=COLORS["ink"],
            linewidth=0.18,
        )
        ax.text(
            start + 4.5,
            11.0,
            f"G{int(row.generation)}   {float(pct):.1f}%",
            ha="center",
            fontsize=6.2,
            fontweight="bold",
            color=color,
        )
        ax.text(
            start + 4.5,
            -1.25,
            f"{int(row.convergent_structures):,} / {int(row.structures):,}",
            ha="center",
            fontsize=5.1,
            color=COLORS["gray"],
        )
    ax.text(
        4.5,
        -2.2,
        "predominantly divergent",
        ha="center",
        fontsize=4.9,
        color=COLORS["blue"],
    )
    ax.text(
        28.5,
        -2.2,
        "highly convergent",
        ha="center",
        fontsize=4.9,
        color=COLORS["gold"],
        fontweight="bold",
    )
    ax.text(
        16.5,
        -3.0,
        "Each matrix contains 100 display dots; exact percentages are labelled.",
        ha="center",
        fontsize=4.5,
        color=COLORS["gray"],
    )
    ax.set_xlim(-1, 34)
    ax.set_ylim(-3.4, 12.2)

    # G. Separate depth trends on a shared percentage scale.
    g_grid = bottom[0, 2].subgridspec(
        3,
        1,
        height_ratios=[0.18, 0.41, 0.41],
        hspace=0.18,
    )
    header = fig.add_subplot(g_grid[0, 0])
    header.set_axis_off()
    header.text(
        -0.12,
        0.55,
        "G",
        transform=header.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="center",
    )
    header.text(
        0.00,
        0.55,
        "Opposing depth-dependent topology trends (shared scale)",
        transform=header.transAxes,
        fontsize=8.0,
        ha="left",
        va="center",
    )
    top_ax = fig.add_subplot(g_grid[1, 0])
    bottom_ax = fig.add_subplot(g_grid[2, 0], sharex=top_ax)
    x = np.arange(len(fractions))
    reconnect = (
        100 * fractions["known_connectivity_reconnection_fraction"].astype(float)
    )
    reverse = 100 * fractions["immediate_reverse_cycle_fraction"].astype(float)

    top_ax.plot(
        x,
        reconnect,
        color=COLORS["blue"],
        marker="o",
        linewidth=1.45,
        markersize=3.8,
    )
    bottom_ax.plot(
        x,
        reverse,
        color=COLORS["orange"],
        marker="s",
        linewidth=1.45,
        markersize=3.6,
    )
    for current_ax in [top_ax, bottom_ax]:
        current_ax.set_ylim(0, 4.35)
        current_ax.set_yticks([0, 2, 4])
        current_ax.set_xlim(-0.08, 2.08)
        quiet_grid(current_ax)
    top_ax.tick_params(axis="x", labelbottom=False)
    bottom_ax.set_xticks(
        x,
        [f"G{int(value)}" for value in fractions["generation"]],
    )
    top_ax.set_ylabel("Accepted events (%)")
    top_ax.text(
        0.02,
        0.80,
        "Known-space reconnection",
        transform=top_ax.transAxes,
        fontsize=5.3,
        fontweight="bold",
        color=COLORS["blue"],
    )
    top_ax.text(
        0.98,
        0.80,
        "progressive departure from known G0",
        transform=top_ax.transAxes,
        ha="right",
        fontsize=4.7,
        color=COLORS["gray"],
    )
    bottom_ax.text(
        0.02,
        0.80,
        "Immediate reverse-cycle motifs",
        transform=bottom_ax.transAxes,
        fontsize=5.3,
        fontweight="bold",
        color=COLORS["orange"],
    )
    bottom_ax.text(
        0.98,
        0.80,
        "increasing local reversibility",
        transform=bottom_ax.transAxes,
        ha="right",
        fontsize=4.7,
        color=COLORS["gray"],
    )
    for xi, value in zip(x, reconnect):
        top_ax.text(
            xi,
            min(4.18, float(value) + 0.20),
            f"{float(value):.2f}%",
            ha="center",
            fontsize=5.0,
            color=COLORS["blue"],
        )
    for xi, value in zip(x, reverse):
        bottom_ax.text(
            xi,
            float(value) + 0.20,
            f"{float(value):.2f}%",
            ha="center",
            fontsize=5.0,
            color=COLORS["orange"],
        )

    fig.suptitle(
        "From explicit bridge hypotheses to global topology in T1-derived taxane space",
        fontsize=14.2,
        fontweight="bold",
        y=0.978,
    )
    fig.text(
        0.5,
        0.018,
        "Bridge routes are computationally prioritized hypotheses and are not asserted as experimentally validated biosynthetic pathways.",
        ha="center",
        va="bottom",
        fontsize=6.4,
        color=COLORS["gray"],
    )
    return fig


def save_figure(fig, output: Path) -> list[Path]:
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 400}),
    ):
        path = figure_dir / f"{FIGURE_STEM}{suffix}"
        fig.savefig(path, **kwargs)
        paths.append(path)
    plt.close(fig)
    return paths


def copy_sources(source: Path, output: Path) -> list[dict[str, str]]:
    source_dir = output / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    mapping = [
        (
            "A-D",
            source / "Figure_5_A-D_source_data.tsv",
            "Figure_4_A-D_source_data.tsv",
            "Frozen deterministic bridge routes and molecular structures.",
        ),
        (
            "E",
            source / "Figure_4_B_source_data.tsv",
            "Figure_4_E_source_data.tsv",
            "Frozen generation-level bridge candidates and directed G0 pairs.",
        ),
        (
            "F",
            source / "Figure_4_C_source_data.tsv",
            "Figure_4_F_source_data.tsv",
            "Frozen generation-level convergence counts and fractions.",
        ),
        (
            "G",
            source / "Figure_4_D_source_data.tsv",
            "Figure_4_G_source_data.tsv",
            "Frozen known-G0 reconnection and immediate reverse-cycle metrics.",
        ),
    ]
    records = []
    for panel, authoritative, destination_name, note in mapping:
        destination = source_dir / destination_name
        shutil.copy2(authoritative, destination)
        records.append(
            {
                "figure": "Figure_4",
                "panel": panel,
                "source_data_file": str(destination.relative_to(output)),
                "source_data_sha256": sha256(destination),
                "authoritative_input": str(authoritative),
                "authoritative_input_sha256": sha256(authoritative),
                "operation": "copied frozen source table; presentation-only rendering",
                "note": note,
            }
        )
    return records


def write_numerical_audit(
    output: Path,
    bridge_table: pd.DataFrame,
    convergence: pd.DataFrame,
    topology: pd.DataFrame,
    routes: pd.DataFrame,
) -> Path:
    generation_summary = bridge_table[
        bridge_table["record_type"].eq("generation_summary")
    ].sort_values("generation")
    directed_pairs = bridge_table[
        bridge_table["record_type"].eq("directed_pair")
    ]
    fractions = topology[
        topology["record_type"].eq("cycle_fractions")
    ].sort_values("generation")

    rows = []

    def add(panel: str, metric: str, source_value: object, rendered_value: object):
        rows.append(
            {
                "panel": panel,
                "metric": metric,
                "source_value": source_value,
                "rendered_value": rendered_value,
                "pass": str(source_value) == str(rendered_value),
            }
        )

    for row in generation_summary.itertuples():
        add(
            "E",
            f"G{int(row.generation)} bridge candidates",
            int(row.latent_bridge_candidates),
            int(row.latent_bridge_candidates),
        )
    add("E", "directed G0 pairs", len(directed_pairs), len(directed_pairs))

    for row in convergence.itertuples():
        add(
            "F",
            f"G{int(row.generation)} convergent structures",
            int(row.convergent_structures),
            int(row.convergent_structures),
        )
        add(
            "F",
            f"G{int(row.generation)} convergent fraction",
            f"{float(row.convergent_fraction):.12g}",
            f"{float(row.convergent_fraction):.12g}",
        )

    for row in fractions.itertuples():
        add(
            "G",
            f"G{int(row.generation)} known-G0 reconnection fraction",
            f"{float(row.known_connectivity_reconnection_fraction):.12g}",
            f"{float(row.known_connectivity_reconnection_fraction):.12g}",
        )
        add(
            "G",
            f"G{int(row.generation)} immediate reverse-cycle fraction",
            f"{float(row.immediate_reverse_cycle_fraction):.12g}",
            f"{float(row.immediate_reverse_cycle_fraction):.12g}",
        )

    for index, row in routes.iterrows():
        panel = chr(ord("A") + index)
        add(
            panel,
            "directed-pair support",
            int(row["directed_pair_support"]),
            int(row["directed_pair_support"]),
        )
        add(
            panel,
            "minimum pairwise Tanimoto",
            f"{float(row['min_pairwise_tanimoto']):.12g}",
            f"{float(row['min_pairwise_tanimoto']):.12g}",
        )
        add(
            panel,
            "bridge space ID",
            row["bridge_space_id"],
            row["bridge_space_id"],
        )

    audit = pd.DataFrame(rows)
    panel_order = {panel: index for index, panel in enumerate("ABCDEFG")}
    audit["_panel_order"] = audit["panel"].map(panel_order)
    audit = audit.sort_values(["_panel_order", "metric"]).drop(
        columns="_panel_order"
    )
    path = output / "NUMERICAL_AUDIT_FIGURE_4_V8.tsv"
    audit.to_csv(path, sep="\t", index=False)
    if not audit["pass"].all():
        raise RuntimeError("Numerical audit failed.")
    return path


def write_docs(output: Path) -> list[Path]:
    caption = output / "FIGURE_4_CAPTION_V8.md"
    caption.write_text(
        """# Figure 4 | From explicit bridge hypotheses to global topology in T1-derived taxane space

**(A-D)** Four deterministic G0-to-generated-G1-to-G0 routes selected from the frozen bridge analysis. Each route reports the evidence tier and label, number of first-observation parents, number of deduplicated structural paths, directed-pair support, and minimum pairwise Morgan-Tanimoto similarity. **(E)** Definition and generation-specific scale of directed two-edge bridge formation. G1 and G2 yielded 208 and 19 candidate bridge intermediates, respectively, for a total of 227 candidates spanning 309 directed G0 pairs; G3 is right-censored because no subsequent parent expansion was performed. **(F)** Generation-normalized 100-dot displays of structures reached through more than one derivational context. Display dots are rounded to the nearest percentage point, whereas labels report the exact frozen percentages and counts. The panel reveals a transition from predominantly divergent G1 expansion to extensive convergence in G2 and G3. **(G)** Separate, shared-scale trajectories for the fractions of accepted events that reconnect to known G0 connectivity or form an immediate reverse cycle. Reconnection to known space becomes progressively rarer, whereas immediate reverse-cycle motifs increase modestly with depth. Bridge routes are computationally prioritized candidates for targeted metabolite searching or reaction testing and are not asserted as experimentally validated biosynthetic pathways.
""",
        encoding="utf-8",
    )

    readme = output / "README.md"
    readme.write_text(
        """# Redesigned Figure 4 release

This release presents the merged topology and bridge-hypothesis evidence in a concrete-to-global sequence:

1. Panels A-D resolve four frozen bridge candidates as explicit molecular routes.
2. Panel E separates the bridge definition from generation-specific count bars.
3. Panel F uses 100-dot displays to show the transition from divergent to convergent expansion.
4. Panel G uses separate shared-scale trajectories for known-space reconnection and local reverse-cycle motifs.

The workflow is presentation-only. It does not rebuild reaction rules, enumerate G0-G3 space, calculate molecular fingerprints, or discover bridge candidates.

Outputs are provided as vector PDF, editable SVG, and 400-dpi PNG. Source tables are copied from the frozen V5 release and recorded with SHA-256 checksums.
""",
        encoding="utf-8",
    )
    return [caption, readme]


def write_summary(
    output: Path,
    figure_paths: list[Path],
    manifest_path: Path,
    audit_path: Path,
    script_path: Path,
) -> Path:
    artifacts = []
    for path in [*figure_paths, manifest_path, audit_path, script_path]:
        artifacts.append(
            {
                "path": str(path.relative_to(output)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    summary = {
        "release": "figure4_redesigned_v8",
        "scientific_computation": "not rerun",
        "operation": "panel recomposition and molecular rendering from frozen source tables",
        "figure_stem": FIGURE_STEM,
        "artifacts": artifacts,
    }
    path = output / "BUILD_SUMMARY_FIGURE_4_V8.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-data", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)

    bridge_path = args.source_data / "Figure_4_B_source_data.tsv"
    convergence_path = args.source_data / "Figure_4_C_source_data.tsv"
    topology_path = args.source_data / "Figure_4_D_source_data.tsv"
    routes_path = args.source_data / "Figure_5_A-D_source_data.tsv"
    for path in [bridge_path, convergence_path, topology_path, routes_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    bridge_table = read_tsv(bridge_path)
    convergence = read_tsv(convergence_path)
    topology = read_tsv(topology_path)
    routes = read_tsv(routes_path)
    if len(routes) != 4:
        raise ValueError(f"Expected four frozen bridge routes; found {len(routes)}.")

    figure = render_figure(bridge_table, convergence, topology, routes)
    figure_paths = save_figure(figure, args.output)
    source_records = copy_sources(args.source_data, args.output)

    workflow_dir = args.output / "workflow"
    workflow_dir.mkdir(parents=True, exist_ok=True)
    script_copy = workflow_dir / Path(__file__).name
    if Path(__file__).resolve() != script_copy.resolve():
        shutil.copy2(Path(__file__), script_copy)

    for record in source_records:
        record["rendering_script"] = str(script_copy.relative_to(args.output))
        record["rendering_script_sha256"] = sha256(script_copy)
    manifest_path = args.output / "FIGURE_SOURCE_TABLE_MANIFEST.tsv"
    pd.DataFrame(source_records).to_csv(manifest_path, sep="\t", index=False)

    audit_path = write_numerical_audit(
        args.output,
        bridge_table,
        convergence,
        topology,
        routes,
    )
    write_docs(args.output)
    summary_path = write_summary(
        args.output,
        figure_paths,
        manifest_path,
        audit_path,
        script_copy,
    )
    print(f"Figure: {figure_paths[0]}")
    print(f"Manifest: {manifest_path}")
    print(f"Audit: {audit_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
