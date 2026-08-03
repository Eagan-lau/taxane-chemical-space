#!/usr/bin/env python3
"""Recompose Figures 2 and 3 from frozen released source data.

This presentation-only workflow does not rebuild reaction rules, enumerate
chemical space, calculate molecular descriptors, or alter frozen scientific
results. It enriches the released network figure with representative molecular
anchors and reframes Figure 3 around cross-substrate grammar transferability.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, Patch, Rectangle
from PIL import Image, ImageChops
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D


DEFAULT_WORK = Path(".")

COLORS = {
    "ink": "#20252b",
    "gray": "#7d8995",
    "grid": "#dfe4e8",
    "blue": "#3f83bd",
    "blue_light": "#d9e8f4",
    "teal": "#2f9183",
    "teal_light": "#d6ebe6",
    "orange": "#df7b19",
    "orange_light": "#f7e2cd",
    "red": "#d44d4f",
    "gold": "#d6a514",
    "g2": "#008c83",
    "g3": "#d95f02",
    "white": "#ffffff",
    "paper": "#fbfcfd",
}

NETWORK_NODE_IDS = [
    "G0_00124",  # taxadiene
    "G0_00620",  # 10-deacetylbaccatin III
    "G0_00429",  # paclitaxel
    "G0_00057",  # 7-(beta-xylopyranozyl)baccatin III
]

NETWORK_ROLES = {
    "G0_00124": "biosynthetic entry scaffold",
    "G0_00620": "oxygenated taxane core",
    "G0_00429": "late-stage taxane product",
    "G0_00057": "glycosylated taxane derivative",
}


def configure_matplotlib() -> None:
    mpl.use("Agg")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.0,
            "axes.titlesize": 10.0,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
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


def formula_label(value: object) -> str:
    import re

    tokens = []
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", str(value)):
        suffix = rf"_{{{count}}}" if count and count != "1" else ""
        tokens.append(rf"\mathrm{{{element}}}{suffix}")
    return "$" + "".join(tokens) + "$" if tokens else str(value)


def elemental_delta_label(value: object) -> str:
    try:
        delta = json.loads(str(value))
    except (json.JSONDecodeError, TypeError):
        return str(value)
    if not delta:
        return "No elemental change"
    order = ("C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I")
    tokens = []
    for element in order:
        amount = int(delta.get(element, 0))
        if amount:
            tokens.append(rf"\Delta {element}={amount:+d}")
    for element in sorted(set(delta) - set(order)):
        amount = int(delta[element])
        if amount:
            tokens.append(rf"\Delta {element}={amount:+d}")
    return "$" + r";\ ".join(tokens) + "$"


def quiet_axes(ax, axis: str = "y") -> None:
    ax.grid(axis=axis, color=COLORS["grid"], linewidth=0.55, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax, label: str, x: float = -0.10, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def molecule_image(
    smiles: str,
    *,
    highlighted_atoms: set[int] | None = None,
    highlighted_bonds: set[tuple[int, int]] | None = None,
    width: int = 1000,
    height: int = 700,
) -> Image.Image:
    highlighted_atoms = highlighted_atoms or set()
    highlighted_bonds = highlighted_bonds or set()
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol, canonOrient=True, clearConfs=True)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    options = drawer.drawOptions()
    options.useBWAtomPalette()
    options.padding = 0.045
    options.bondLineWidth = 1.8
    options.minFontSize = 18
    options.maxFontSize = 30
    options.fixedBondLength = 29
    options.highlightRadius = 0.31
    options.fillHighlights = True
    options.continuousHighlight = False
    options.addStereoAnnotation = False
    highlighted_bond_ids = []
    for bond in mol.GetBonds():
        key = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
        if key in highlighted_bonds:
            highlighted_bond_ids.append(bond.GetIdx())
    color = mpl.colors.to_rgb(COLORS["orange"])
    drawer.DrawMolecule(
        mol,
        highlightAtoms=sorted(highlighted_atoms),
        highlightBonds=highlighted_bond_ids,
        highlightAtomColors={idx: color for idx in highlighted_atoms},
        highlightBondColors={idx: color for idx in highlighted_bond_ids},
    )
    drawer.FinishDrawing()
    image = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")
    background = Image.new("RGB", image.size, "white")
    bbox = ImageChops.difference(image, background).getbbox()
    if bbox:
        margin = 14
        image = image.crop(
            (
                max(0, bbox[0] - margin),
                max(0, bbox[1] - margin),
                min(image.width, bbox[2] + margin),
                min(image.height, bbox[3] + margin),
            )
        )
    return image


def draw_molecule(
    ax,
    smiles: str,
    *,
    highlighted_atoms: set[int] | None = None,
    highlighted_bonds: set[tuple[int, int]] | None = None,
) -> None:
    ax.set_axis_off()
    image = molecule_image(
        smiles,
        highlighted_atoms=highlighted_atoms,
        highlighted_bonds=highlighted_bonds,
    )
    ax.imshow(image, interpolation="lanczos", aspect="equal")


def _bond_order(bond) -> float:
    return 1.5 if bond.GetIsAromatic() else float(bond.GetBondTypeAsDouble())


def molecule_difference(
    source_smiles: str, product_smiles: str
) -> tuple[set[int], set[int], set[tuple[int, int]], set[tuple[int, int]]]:
    """Return display-only MCS differences for reaction-center highlighting."""
    source = Chem.MolFromSmiles(source_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    if source is None or product is None:
        return set(), set(), set(), set()
    result = rdFMCS.FindMCS(
        [source, product],
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchChiralTag=False,
        timeout=10,
    )
    query = Chem.MolFromSmarts(result.smartsString) if result.smartsString else None
    if query is None:
        return (
            set(range(source.GetNumAtoms())),
            set(range(product.GetNumAtoms())),
            set(),
            set(),
        )
    source_match = source.GetSubstructMatch(query)
    product_match = product.GetSubstructMatch(query)
    mapping = dict(zip(source_match, product_match))
    inverse = {value: key for key, value in mapping.items()}
    changed_source = set(range(source.GetNumAtoms())) - set(source_match)
    changed_product = set(range(product.GetNumAtoms())) - set(product_match)
    changed_source_bonds: set[tuple[int, int]] = set()
    changed_product_bonds: set[tuple[int, int]] = set()
    for source_idx, product_idx in mapping.items():
        left = source.GetAtomWithIdx(source_idx)
        right = product.GetAtomWithIdx(product_idx)
        if (
            left.GetAtomicNum(),
            left.GetFormalCharge(),
            left.GetIsAromatic(),
            left.GetDegree(),
        ) != (
            right.GetAtomicNum(),
            right.GetFormalCharge(),
            right.GetIsAromatic(),
            right.GetDegree(),
        ):
            changed_source.add(source_idx)
            changed_product.add(product_idx)
    for bond in source.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a not in mapping or b not in mapping:
            changed_source_bonds.add(tuple(sorted((a, b))))
            changed_source.update((a, b))
            continue
        product_bond = product.GetBondBetweenAtoms(mapping[a], mapping[b])
        if product_bond is None or _bond_order(product_bond) != _bond_order(bond):
            changed_source_bonds.add(tuple(sorted((a, b))))
            changed_source.update((a, b))
            changed_product.update((mapping[a], mapping[b]))
            if product_bond is not None:
                changed_product_bonds.add(
                    tuple(sorted((mapping[a], mapping[b])))
                )
    for bond in product.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a not in inverse or b not in inverse:
            changed_product_bonds.add(tuple(sorted((a, b))))
            changed_product.update((a, b))
    return (
        changed_source,
        changed_product,
        changed_source_bonds,
        changed_product_bonds,
    )


def intensity_rgba(
    intensity: np.ndarray, color: str, maximum_alpha: float
) -> np.ndarray:
    rgba = np.zeros((*intensity.shape, 4), dtype=np.float32)
    rgba[..., :3] = mpl.colors.to_rgb(color)
    rgba[..., 3] = np.clip(intensity, 0, 1) * maximum_alpha
    return rgba


def load_inputs(work: Path) -> dict[str, object]:
    v5 = (
        work
        / "release/manuscript_v5_argument_driven_figures_20260730"
    )
    network_path = (
        v5 / "source_data/Figure_2_A_complete_network_source_data.tsv"
    )
    density_path = v5 / "source_data/Figure_2_A_density_layers.npz"
    figure1_path = (
        work
        / "release/figure1_redesign_v6_20260730"
        / "source_data/Figure_1C_grammar_examples.tsv"
    )
    projected_path = v5 / "source_data/Figure_3_A-C_source_data.tsv"
    locality_path = v5 / "source_data/Figure_3_D_source_data.tsv"
    functional_path = v5 / "source_data/Figure_3_E_source_data.tsv"
    elemental_path = v5 / "source_data/Figure_3_F_source_data.tsv"

    network = read_tsv(network_path)
    density = np.load(density_path, allow_pickle=False)
    anchors = read_tsv(figure1_path)
    projected = read_tsv(projected_path)
    locality = read_tsv(locality_path)
    functional = read_tsv(functional_path)
    elemental = read_tsv(elemental_path)

    nodes = network[network["record_type"].eq("node")].copy()
    nodes["generation"] = pd.to_numeric(nodes["generation"]).astype(int)
    g0 = nodes[nodes["generation"].eq(0)].copy()
    representative = g0[g0["space_id"].isin(NETWORK_NODE_IDS)].copy()
    representative["formula"] = representative["smiles"].map(
        lambda value: rdMolDescriptors.CalcMolFormula(
            Chem.MolFromSmiles(str(value))
        )
    )
    representative["interpretive_role"] = representative["space_id"].map(
        NETWORK_ROLES
    )
    representative["display_order"] = representative["space_id"].map(
        {space_id: index for index, space_id in enumerate(NETWORK_NODE_IDS)}
    )
    representative = representative.sort_values("display_order")

    shared_rules = [
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000001",
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000003",
    ]
    anchor_rows = anchors[
        anchors["smarts_rule_id"].isin(shared_rules)
    ].rename(
        columns={
            "display_name": "anchor_display_name",
            "local_edit": "anchor_local_edit",
            "reaction_smarts": "anchor_reaction_smarts",
            "enzyme": "anchor_enzyme",
            "ec": "anchor_ec",
            "substrate_smiles": "anchor_substrate_smiles",
            "product_smiles": "anchor_product_smiles",
        }
    )
    projected_rows = projected[
        projected["grammar_rule_id"].isin(shared_rules)
    ].copy()
    migration = anchor_rows.merge(
        projected_rows,
        left_on="smarts_rule_id",
        right_on="grammar_rule_id",
        how="inner",
        validate="one_to_one",
    )
    migration["display_order"] = migration["smarts_rule_id"].map(
        {rule_id: index for index, rule_id in enumerate(shared_rules)}
    )
    migration = migration.sort_values("display_order")

    return {
        "network": network,
        "density": density,
        "representative": representative,
        "migration": migration,
        "locality": locality,
        "functional": functional,
        "elemental": elemental,
        "paths": {
            "network": network_path,
            "density": density_path,
            "anchors": figure1_path,
            "projected": projected_path,
            "locality": locality_path,
            "functional": functional_path,
            "elemental": elemental_path,
        },
    }


def split_network(network: pd.DataFrame) -> dict[str, pd.DataFrame]:
    nodes = network[network["record_type"].eq("node")].copy()
    edges = network[network["record_type"].eq("edge")].copy()
    nodes["generation"] = pd.to_numeric(nodes["generation"]).astype(int)
    edges["generation"] = pd.to_numeric(edges["generation"]).astype(int)
    g0 = nodes[nodes["generation"].eq(0)].copy()
    g1 = nodes[nodes["generation"].eq(1)].copy()
    g0_edges = edges[edges["generation"].eq(0)].copy()
    g1_edges = edges[edges["generation"].eq(1)].copy()
    for frame in (g0, g1, g0_edges, g1_edges):
        for column in (
            "component_index",
            "layout_component_index",
            "x",
            "y",
            "layout_x",
            "layout_y",
            "source_x",
            "source_y",
            "target_x",
            "target_y",
        ):
            if column in frame:
                frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return {
        "g0": g0,
        "g1": g1,
        "g0_edges": g0_edges,
        "g1_edges": g1_edges,
    }


def draw_network_axis(
    ax,
    network: pd.DataFrame,
    density,
    selected: pd.DataFrame,
) -> None:
    parts = split_network(network)
    g0 = parts["g0"]
    g1 = parts["g1"]
    g0_edges = parts["g0_edges"]
    g1_edges = parts["g1_edges"]
    bounds = np.asarray(density["bounds"], dtype=float)
    extent = (bounds[0], bounds[1], bounds[2], bounds[3])
    ax.set_facecolor(COLORS["paper"])
    ax.imshow(
        intensity_rgba(
            density["g3_component_normalized_intensity"],
            COLORS["g3"],
            0.27,
        ),
        extent=extent,
        origin="lower",
        interpolation="bilinear",
        zorder=0,
    )
    ax.imshow(
        intensity_rgba(
            density["g2_component_normalized_intensity"],
            COLORS["g2"],
            0.34,
        ),
        extent=extent,
        origin="lower",
        interpolation="bilinear",
        zorder=1,
    )

    core_component = int(g0.groupby("component_index").size().idxmax())
    core_g0 = g0[g0["component_index"].eq(core_component)]
    embedded_g0 = g0[~g0["component_index"].eq(core_component)]
    core_g1 = g1[g1["layout_component_index"].eq(core_component)]
    embedded_g1 = g1[~g1["layout_component_index"].eq(core_component)]
    core_g0_edges = g0_edges[
        g0_edges["layout_component_index"].eq(core_component)
    ]
    embedded_g0_edges = g0_edges[
        ~g0_edges["layout_component_index"].eq(core_component)
    ]
    core_g1_edges = g1_edges[
        g1_edges["layout_component_index"].eq(core_component)
    ]
    embedded_g1_edges = g1_edges[
        ~g1_edges["layout_component_index"].eq(core_component)
    ]

    def segments(frame: pd.DataFrame) -> list[list[tuple[float, float]]]:
        return [
            [(row.source_x, row.source_y), (row.target_x, row.target_y)]
            for row in frame.itertuples(index=False)
        ]

    ax.add_collection(
        LineCollection(
            segments(embedded_g1_edges),
            colors=COLORS["blue"],
            linewidths=0.08,
            alpha=0.012,
            zorder=2,
        )
    )
    ax.add_collection(
        LineCollection(
            segments(embedded_g0_edges),
            colors=COLORS["ink"],
            linewidths=0.18,
            alpha=0.045,
            zorder=3,
        )
    )
    ax.add_collection(
        LineCollection(
            segments(core_g1_edges),
            colors=COLORS["blue"],
            linewidths=0.14,
            alpha=0.06,
            zorder=3,
        )
    )
    ax.add_collection(
        LineCollection(
            segments(core_g0_edges),
            colors=COLORS["ink"],
            linewidths=0.34,
            alpha=0.23,
            zorder=4,
        )
    )
    marker_size = 1.55
    ax.scatter(
        embedded_g1["layout_x"],
        embedded_g1["layout_y"],
        s=marker_size,
        color=COLORS["blue"],
        alpha=0.25,
        edgecolor="none",
        zorder=4,
    )
    ax.scatter(
        embedded_g0["x"],
        embedded_g0["y"],
        s=marker_size,
        color="#e33438",
        alpha=0.38,
        edgecolor="none",
        zorder=5,
    )
    ax.scatter(
        core_g1["layout_x"],
        core_g1["layout_y"],
        s=marker_size,
        color=COLORS["blue"],
        alpha=0.76,
        edgecolor="none",
        zorder=5,
    )
    ax.scatter(
        core_g0["x"],
        core_g0["y"],
        s=marker_size,
        color="#e33438",
        alpha=0.97,
        edgecolor="white",
        linewidth=0.06,
        zorder=6,
    )

    for number, row in enumerate(selected.itertuples(index=False), start=1):
        ax.scatter(
            [row.x],
            [row.y],
            s=46,
            facecolor="none",
            edgecolor=COLORS["ink"],
            linewidth=0.8,
            zorder=12,
        )
        ax.text(
            row.x,
            row.y,
            str(number),
            ha="center",
            va="center",
            fontsize=4.5,
            fontweight="bold",
            color=COLORS["ink"],
            zorder=13,
        )

    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[3], bounds[2])
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#cfd6dc")
        spine.set_linewidth(0.65)


def draw_molecule_card(ax, row, number: int) -> None:
    ax.set_axis_off()
    ax.add_patch(
        Rectangle(
            (0.01, 0.01),
            0.98,
            0.98,
            transform=ax.transAxes,
            facecolor=COLORS["white"],
            edgecolor="#ccd4da",
            linewidth=0.8,
        )
    )
    ax.text(
        0.06,
        0.93,
        f"{number}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.0,
        fontweight="bold",
        color=COLORS["white"],
        bbox={
            "boxstyle": "circle,pad=0.28",
            "facecolor": COLORS["red"],
            "edgecolor": "none",
        },
    )
    ax.text(
        0.18,
        0.94,
        wrap(row.molecule_names, 25),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.7,
        fontweight="bold",
        color=COLORS["ink"],
    )
    structure_ax = ax.inset_axes([0.06, 0.27, 0.88, 0.54])
    draw_molecule(structure_ax, row.smiles)
    ax.text(
        0.50,
        0.18,
        formula_label(row.formula),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.5,
        color=COLORS["ink"],
    )
    ax.text(
        0.50,
        0.07,
        row.interpretive_role,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.3,
        color=COLORS["gray"],
    )


def build_figure2(data: dict[str, object], output: Path) -> list[Path]:
    network = data["network"]
    density = data["density"]
    selected = data["representative"]
    fig = plt.figure(figsize=(15.0, 9.8))
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=[2.25, 10.5, 2.25],
        hspace=0.14,
        wspace=0.05,
    )
    left_axes = [fig.add_subplot(grid[index, 0]) for index in range(2)]
    network_ax = fig.add_subplot(grid[:, 1])
    right_axes = [fig.add_subplot(grid[index, 2]) for index in range(2)]
    card_axes = [left_axes[0], left_axes[1], right_axes[0], right_axes[1]]

    draw_network_axis(network_ax, network, density, selected)
    for index, (card_ax, row) in enumerate(
        zip(card_axes, selected.itertuples(index=False)), start=1
    ):
        draw_molecule_card(card_ax, row, index)
        card_anchor = (0.99, 0.50) if index <= 2 else (0.01, 0.50)
        connection = ConnectionPatch(
            xyA=(row.x, row.y),
            coordsA=network_ax.transData,
            xyB=card_anchor,
            coordsB=card_ax.transAxes,
            arrowstyle="-",
            linewidth=0.65,
            color=COLORS["gray"],
            alpha=0.72,
            zorder=10,
        )
        fig.add_artist(connection)

    summary = (
        network[network["record_type"].eq("generation_summary")]
        .set_index("generation")["structure_count"]
        .astype(int)
    )
    edge_count = int(network["record_type"].eq("edge").sum())
    network_ax.text(
        0.015,
        0.985,
        f"{int(summary.sum()):,} unique structures  |  "
        f"{edge_count:,} explicit G0/G1 edge records",
        transform=network_ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        color=COLORS["ink"],
        bbox={
            "boxstyle": "round,pad=0.40",
            "facecolor": COLORS["white"],
            "edgecolor": "#d3d9de",
            "linewidth": 0.65,
            "alpha": 0.96,
        },
        zorder=20,
    )
    network_ax.text(
        0.015,
        0.942,
        "G0 648  |  G1 15,801  |  G2 223,823  |  "
        "G3 2,362,766 (exploratory)",
        transform=network_ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.9,
        color=COLORS["gray"],
        bbox={
            "boxstyle": "round,pad=0.34",
            "facecolor": COLORS["white"],
            "edgecolor": "#e0e4e8",
            "linewidth": 0.5,
            "alpha": 0.94,
        },
        zorder=20,
    )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="#e33438",
            markeredgecolor="none",
            markersize=5.3,
            label="Known taxane seed (G0)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=COLORS["blue"],
            markeredgecolor="none",
            markersize=5.3,
            label="One-rule intermediate (G1)",
        ),
        Line2D(
            [0],
            [0],
            color=COLORS["ink"],
            lw=0.9,
            alpha=0.55,
            label="Established G0 edge",
        ),
        Line2D(
            [0],
            [0],
            color=COLORS["blue"],
            lw=0.9,
            alpha=0.7,
            label="G0 to G1 derivation",
        ),
        Patch(
            facecolor=COLORS["g2"],
            alpha=0.62,
            label="Two-rule descendant density (G2)",
        ),
        Patch(
            facecolor=COLORS["g3"],
            alpha=0.54,
            label="Exploratory three-rule density (G3)",
        ),
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        fontsize=7.0,
        columnspacing=1.6,
        handletextpad=0.55,
    )
    fig.suptitle(
        "Multiscale organization of T1-derived taxane chemical space",
        fontsize=15.0,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.948,
        "Node-resolved G0/G1 topology, complete G2/G3 descendant density, "
        "and representative molecular anchors",
        ha="center",
        va="center",
        fontsize=8.0,
        color=COLORS["gray"],
    )
    fig.subplots_adjust(
        left=0.018,
        right=0.982,
        top=0.915,
        bottom=0.095,
    )
    return save_figure(
        fig,
        output,
        "Figure_2_T1_Taxane_Chemical_Space_V6",
    )


def add_axis_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    linewidth: float = 1.0,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=linewidth,
            color=color,
            clip_on=False,
        )
    )


def draw_migration_row(
    ax,
    row,
    *,
    y: float,
    color: str,
    rule_number: int,
) -> None:
    row_height = 0.34
    ax.add_patch(
        Rectangle(
            (0.015, y),
            0.97,
            row_height,
            transform=ax.transAxes,
            facecolor="#ffffff",
            edgecolor="#d8dde2",
            linewidth=0.65,
        )
    )
    ax.add_patch(
        Rectangle(
            (0.015, y),
            0.010,
            row_height,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
        )
    )

    anchor_x, anchor_w = 0.045, 0.17
    smarts_x, smarts_w = 0.245, 0.21
    source_x, molecule_w = 0.485, 0.17
    product_x = 0.695
    label_x = 0.885
    center_y = y + row_height / 2

    ax.text(
        anchor_x,
        y + row_height * 0.76,
        "CURATED PATHWAY ANCHOR",
        transform=ax.transAxes,
        fontsize=5.2,
        fontweight="bold",
        color=COLORS["gray"],
        ha="left",
        va="center",
    )
    ax.text(
        anchor_x,
        y + row_height * 0.50,
        row.anchor_local_edit.replace("  ", " "),
        transform=ax.transAxes,
        fontsize=7.5,
        fontweight="bold",
        color=color,
        ha="left",
        va="center",
    )
    ax.text(
        anchor_x,
        y + row_height * 0.25,
        f"{row.anchor_enzyme}  |  EC {row.anchor_ec}",
        transform=ax.transAxes,
        fontsize=5.1,
        color=COLORS["gray"],
        ha="left",
        va="center",
    )

    ax.add_patch(
        Rectangle(
            (smarts_x, y + 0.055),
            smarts_w,
            row_height - 0.11,
            transform=ax.transAxes,
            facecolor="#f7f9fa",
            edgecolor=color,
            linewidth=0.85,
        )
    )
    ax.text(
        smarts_x + smarts_w / 2,
        y + row_height * 0.72,
        f"GENERALIZED T1 RULE {rule_number}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.3,
        fontweight="bold",
        color=color,
    )
    ax.text(
        smarts_x + smarts_w / 2,
        y + row_height * 0.44,
        wrap(row.anchor_reaction_smarts, 35),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.0,
        family="DejaVu Sans Mono",
        color=COLORS["ink"],
    )
    ax.text(
        smarts_x + smarts_w / 2,
        y + row_height * 0.18,
        "local context abstracted",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=4.8,
        color=COLORS["gray"],
    )

    source_atoms, target_atoms, source_bonds, target_bonds = (
        molecule_difference(row.source_smiles, row.target_smiles)
    )
    source_ax = ax.inset_axes(
        [source_x, y + 0.045, molecule_w, row_height - 0.09]
    )
    product_ax = ax.inset_axes(
        [product_x, y + 0.045, molecule_w, row_height - 0.09]
    )
    draw_molecule(
        source_ax,
        row.source_smiles,
        highlighted_atoms=source_atoms,
        highlighted_bonds=source_bonds,
    )
    draw_molecule(
        product_ax,
        row.target_smiles,
        highlighted_atoms=target_atoms,
        highlighted_bonds=target_bonds,
    )
    ax.text(
        source_x + molecule_w / 2,
        y + row_height * 0.92,
        f"G0 | {wrap(row.source_name, 25)}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.3,
        fontweight="bold",
        color=COLORS["ink"],
    )
    ax.text(
        product_x + molecule_w / 2,
        y + row_height * 0.92,
        "GENERATED G1 PRODUCT",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.2,
        fontweight="bold",
        color=COLORS["ink"],
    )
    ax.text(
        source_x + molecule_w / 2,
        y + row_height * 0.08,
        formula_label(row.source_formula),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.5,
    )
    ax.text(
        product_x + molecule_w / 2,
        y + row_height * 0.08,
        formula_label(row.target_formula),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.5,
    )

    ax.text(
        label_x,
        y + row_height * 0.64,
        "TRANSFERRED EDIT",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=4.9,
        fontweight="bold",
        color=COLORS["gray"],
    )
    ax.text(
        label_x,
        y + row_height * 0.43,
        row.anchor_local_edit.replace("  ", " "),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.0,
        fontweight="bold",
        color=color,
    )
    ax.text(
        label_x,
        y + row_height * 0.20,
        f"Tanimoto {float(row.source_product_tanimoto):.3f}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=4.9,
        color=COLORS["gray"],
    )

    add_axis_arrow(
        ax,
        (anchor_x + anchor_w, center_y),
        (smarts_x - 0.012, center_y),
        color=color,
    )
    add_axis_arrow(
        ax,
        (smarts_x + smarts_w + 0.008, center_y),
        (source_x - 0.012, center_y),
        color=color,
    )
    add_axis_arrow(
        ax,
        (source_x + molecule_w + 0.008, center_y),
        (product_x - 0.012, center_y),
        color=color,
    )


def build_figure3(data: dict[str, object], output: Path) -> list[Path]:
    migration = data["migration"]
    locality = data["locality"]
    functional = data["functional"]
    elemental = data["elemental"]

    fig = plt.figure(figsize=(14.2, 8.9))
    grid = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.22, 1.0],
        hspace=0.22,
        wspace=0.52,
    )

    ax_a = fig.add_subplot(grid[0, :])
    ax_a.set_axis_off()
    ax_a.set_title(
        "Cross-substrate transferability of curated T1 reaction grammar",
        fontsize=10.5,
        pad=6,
    )
    panel_label(ax_a, "A", x=-0.025, y=1.00)
    row_colors = [COLORS["blue"], COLORS["orange"]]
    row_positions = [0.56, 0.13]
    for index, (row, y, color) in enumerate(
        zip(
            migration.itertuples(index=False),
            row_positions,
            row_colors,
        ),
        start=1,
    ):
        draw_migration_row(
            ax_a,
            row,
            y=y,
            color=color,
            rule_number=index,
        )

    ax_b = fig.add_subplot(grid[1, 0])
    pivot = locality.pivot(
        index="generation",
        columns="edit_locality",
        values="derivation_event_count",
    ).fillna(0)
    pivot = pivot.reindex(columns=["0", "1", "2-3", ">3"], fill_value=0)
    fractions = pivot.div(pivot.sum(axis=1), axis=0)
    bottom = np.zeros(len(fractions))
    categories = [
        ("0", COLORS["gray"]),
        ("1", COLORS["blue"]),
        ("2-3", COLORS["teal"]),
        (">3", COLORS["orange"]),
    ]
    for category, color in categories:
        ax_b.bar(
            [f"G{int(value)}" for value in fractions.index],
            fractions[category],
            bottom=bottom,
            label=category,
            color=color,
            edgecolor=COLORS["ink"],
            linewidth=0.25,
        )
        bottom += fractions[category].to_numpy()
    totals = pivot.sum(axis=1).astype(int)
    for index, generation in enumerate(fractions.index):
        ax_b.text(
            index,
            0.885,
            f"{100 * fractions.loc[generation, '1']:.1f}%\n1-atom",
            ha="center",
            va="center",
            fontsize=5.6,
            color=COLORS["white"],
            fontweight="bold",
        )
        ax_b.text(
            index,
            -0.075,
            f"n={int(totals.loc[generation]):,}",
            ha="center",
            va="top",
            fontsize=5.1,
            color=COLORS["gray"],
        )
    ax_b.set_ylim(0, 1)
    ax_b.set_ylabel("Fraction of accepted events")
    ax_b.set_title("Reaction-edit locality")
    ax_b.legend(
        title="changed source atoms",
        frameon=False,
        ncol=2,
        fontsize=5.1,
        title_fontsize=5.2,
        loc="lower left",
    )
    quiet_axes(ax_b, "y")
    panel_label(ax_b, "B")

    ax_c = fig.add_subplot(grid[1, 1])
    functional_summary = functional.sort_values(
        "derivation_event_count", ascending=True
    )
    functional_labels = {
        "free_hydroxyl:+1": "free hydroxyl +1",
        "free_hydroxyl:-1;ester:+1;carboxylic_acid_or_carboxylate:+1;ether:+1": (
            "free hydroxyl -1; ester +1;\n"
            "carboxylic acid/carboxylate +1; ether +1"
        ),
        "free_hydroxyl:-1;ester:+1;ether:+1": (
            "free hydroxyl -1; ester +1; ether +1"
        ),
        "free_hydroxyl:+1;ester:-1;ether:-1": (
            "free hydroxyl +1; ester -1; ether -1"
        ),
        "free_hydroxyl:-1;ketone_or_aldehyde:+1": (
            "free hydroxyl -1; ketone/aldehyde +1"
        ),
        "no_counted_functional_state_change": "no counted functional-state change",
        "free_hydroxyl:-2;carboxylic_acid_or_carboxylate:+1": (
            "free hydroxyl -2; carboxylic acid/carboxylate +1"
        ),
    }
    labels = [
        functional_labels.get(value, str(value).replace("_", " "))
        for value in functional_summary["functional_state_transition"]
    ]
    ax_c.barh(
        labels,
        functional_summary["derivation_event_count"],
        color=COLORS["teal"],
    )
    ax_c.set_xscale("log")
    ax_c.set_xlabel("Accepted derivation events (log scale)")
    ax_c.set_title("Dominant functional-state transitions")
    ax_c.tick_params(axis="y", labelsize=5.0)
    quiet_axes(ax_c, "x")
    panel_label(ax_c, "C")

    ax_d = fig.add_subplot(grid[1, 2])
    elemental_summary = elemental.sort_values("derivation_event_count")
    ax_d.barh(
        [
            elemental_delta_label(value)
            for value in elemental_summary["observed_element_delta"]
        ],
        elemental_summary["derivation_event_count"],
        color=COLORS["orange"],
    )
    ax_d.set_xscale("log")
    ax_d.set_xlabel("Accepted derivation events (log scale)")
    ax_d.set_title("Elemental edit signatures")
    ax_d.tick_params(axis="y", labelsize=5.1)
    quiet_axes(ax_d, "x")
    panel_label(ax_d, "D")

    generation_totals = (
        locality.groupby("generation")["derivation_event_count"]
        .sum()
        .astype(int)
    )
    fig.suptitle(
        "Transferable local chemistry underlying T1 space expansion",
        fontsize=14.5,
        fontweight="bold",
        y=0.986,
    )
    fig.text(
        0.5,
        0.947,
        f"{int(generation_totals.sum()):,} accepted derivation events  |  "
        f"G1 {int(generation_totals.loc[1]):,}  |  "
        f"G2 {int(generation_totals.loc[2]):,}  |  "
        f"G3 {int(generation_totals.loc[3]):,}",
        ha="center",
        va="center",
        fontsize=7.2,
        color=COLORS["gray"],
        bbox={
            "boxstyle": "round,pad=0.34",
            "facecolor": COLORS["white"],
            "edgecolor": "#d8dde2",
            "linewidth": 0.55,
        },
    )
    fig.subplots_adjust(
        left=0.075,
        right=0.985,
        top=0.900,
        bottom=0.085,
    )
    return save_figure(
        fig,
        output,
        "Figure_3_T1_Chemical_Regularities_V6",
    )


def save_figure(fig, output: Path, stem: str) -> list[Path]:
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 400}),
    ):
        path = figure_dir / f"{stem}{suffix}"
        fig.savefig(path, **kwargs)
        paths.append(path)
    plt.close(fig)
    return paths


def validate(data: dict[str, object]) -> pd.DataFrame:
    network = data["network"]
    locality = data["locality"]
    summary = (
        network[network["record_type"].eq("generation_summary")]
        .set_index("generation")["structure_count"]
        .astype(int)
    )
    events = (
        locality.groupby("generation")["derivation_event_count"]
        .sum()
        .astype(int)
    )
    observations = {
        "G0_unique_structures": int(summary.loc["G0"]),
        "G1_unique_structures": int(summary.loc["G1"]),
        "G2_unique_structures": int(summary.loc["G2"]),
        "G3_unique_structures": int(summary.loc["G3"]),
        "total_unique_structures": int(summary.sum()),
        "explicit_G0_G1_edge_records": int(
            network["record_type"].eq("edge").sum()
        ),
        "G1_accepted_derivation_events": int(events.loc[1]),
        "G2_accepted_derivation_events": int(events.loc[2]),
        "G3_accepted_derivation_events": int(events.loc[3]),
        "total_accepted_derivation_events": int(events.sum()),
        "network_molecular_anchor_count": len(data["representative"]),
        "grammar_migration_chain_count": len(data["migration"]),
    }
    expected = {
        "G0_unique_structures": 648,
        "G1_unique_structures": 15801,
        "G2_unique_structures": 223823,
        "G3_unique_structures": 2362766,
        "total_unique_structures": 2603038,
        "explicit_G0_G1_edge_records": 17661,
        "G1_accepted_derivation_events": 18504,
        "G2_accepted_derivation_events": 490855,
        "G3_accepted_derivation_events": 7418850,
        "total_accepted_derivation_events": 7928209,
        "network_molecular_anchor_count": 4,
        "grammar_migration_chain_count": 2,
    }
    rows = []
    for claim, expected_value in expected.items():
        observed = observations[claim]
        rows.append(
            {
                "claim_id": claim,
                "expected_frozen_value": expected_value,
                "observed_value": observed,
                "status": "PASS" if observed == expected_value else "FAIL",
            }
        )
    audit = pd.DataFrame(rows)
    if not audit["status"].eq("PASS").all():
        raise RuntimeError(
            "Numerical audit failed:\n"
            + audit[audit["status"].ne("PASS")].to_string(index=False)
        )
    return audit


def write_release(
    data: dict[str, object],
    audit: pd.DataFrame,
    output: Path,
    script_path: Path,
    figure_paths: list[Path],
) -> None:
    source_dir = output / "source_data"
    workflow_dir = output / "workflow"
    source_dir.mkdir(parents=True, exist_ok=True)
    workflow_dir.mkdir(parents=True, exist_ok=True)

    source_tables = {
        "Figure_2_molecular_anchors.tsv": data["representative"],
        "Figure_3_A_grammar_migration_chains.tsv": data["migration"],
        "Figure_3_B_edit_locality.tsv": data["locality"],
        "Figure_3_C_functional_state_transitions.tsv": data["functional"],
        "Figure_3_D_elemental_edit_signatures.tsv": data["elemental"],
    }
    for name, frame in source_tables.items():
        frame.to_csv(source_dir / name, sep="\t", index=False)

    network_source = source_dir / "Figure_2_complete_network_source.tsv"
    density_source = source_dir / "Figure_2_complete_density_layers.npz"
    shutil.copy2(data["paths"]["network"], network_source)
    shutil.copy2(data["paths"]["density"], density_source)

    audit.to_csv(output / "NUMERICAL_AUDIT_V6.tsv", sep="\t", index=False)
    shutil.copy2(script_path, workflow_dir / script_path.name)

    caption = """# Figure 2 | Multiscale organization of T1-derived taxane chemical space

Complete integrated representation of the frozen T1-derived taxane reaction-grammar space. All 648 known taxanes (G0, red) and 15,801 unique one-rule products (G1, blue) are represented as individual nodes. Established G0 edges and accepted G0-to-G1 derivations are drawn explicitly. The 223,823 G2 structures and 2,362,766 exploratory G3 structures are represented as component-normalized descendant-density layers. Four representative G0 structures are linked to their actual network coordinates: taxadiene, 10-deacetylbaccatin III, paclitaxel, and 7-(beta-xylopyranozyl)baccatin III. Molecular callouts provide chemical context but do not define layout axes or component membership. Teal and orange density encode G2 and G3, respectively. Layout coordinates organize network topology and are not chemical-distance axes.

# Figure 3 | Transferable local chemistry underlying T1 space expansion

**(A)** Two representative cross-substrate grammar-transfer chains. Curated Taxol-pathway reaction anchors are abstracted into directional T1 reaction SMARTS and applied to different known G0 taxanes to generate accepted G1 products. Orange highlighting is a depiction aid for changed molecular regions. These products are reaction-grammar-accessible hypotheses rather than experimentally observed metabolites. **(B)** Distribution of reaction-edit locality across generations; percentages identify events with one changed source atom and labels report generation-specific accepted-event denominators. **(C)** Dominant structure-derived functional-state transitions. **(D)** Most frequent elemental edit signatures. The header reports the complete 7,928,209 accepted derivation events used for the generation-level chemical summaries.
"""
    (output / "FIGURE_CAPTIONS_V6.md").write_text(caption, encoding="utf-8")

    manifest_rows = []
    for figure, panel, local_file, authoritative in [
        (
            "Figure_2",
            "network",
            network_source,
            data["paths"]["network"],
        ),
        (
            "Figure_2",
            "density",
            density_source,
            data["paths"]["density"],
        ),
        (
            "Figure_2",
            "molecular anchors",
            source_dir / "Figure_2_molecular_anchors.tsv",
            data["paths"]["network"],
        ),
        (
            "Figure_3",
            "A",
            source_dir / "Figure_3_A_grammar_migration_chains.tsv",
            str(data["paths"]["anchors"])
            + ";"
            + str(data["paths"]["projected"]),
        ),
        (
            "Figure_3",
            "B",
            source_dir / "Figure_3_B_edit_locality.tsv",
            data["paths"]["locality"],
        ),
        (
            "Figure_3",
            "C",
            source_dir / "Figure_3_C_functional_state_transitions.tsv",
            data["paths"]["functional"],
        ),
        (
            "Figure_3",
            "D",
            source_dir / "Figure_3_D_elemental_edit_signatures.tsv",
            data["paths"]["elemental"],
        ),
    ]:
        manifest_rows.append(
            {
                "figure": figure,
                "panel": panel,
                "source_data_file": str(local_file.relative_to(output)),
                "source_data_sha256": sha256(local_file),
                "authoritative_input": str(authoritative),
                "operation": "display-only extraction or recomposition",
                "rendering_script": f"workflow/{script_path.name}",
            }
        )
    pd.DataFrame(manifest_rows).to_csv(
        output / "FIGURE_SOURCE_MANIFEST_V6.tsv",
        sep="\t",
        index=False,
    )

    summary = {
        "release": "Figures 2 and 3 V6",
        "scientific_recalculation_performed": False,
        "audit_status": "PASS",
        "figures": [
            {
                "path": str(path.relative_to(output)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in figure_paths
        ],
    }
    (output / "BUILD_SUMMARY_V6.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    readme = """# Figures 2 and 3 V6

This presentation-only release enriches the frozen complete T1 chemical-space
network with four representative molecular anchors and reframes Figure 3
around two cross-substrate reaction-grammar transfer chains. No reaction-rule
construction, molecular enumeration, fingerprint calculation, descriptor
calculation, or scientific benchmark was rerun.

Figure 2 retains the complete released G0/G1 node-resolved topology and G2/G3
component-normalized density arrays. Figure 3 retains the released locality,
functional-state, and elemental-edit counts. All displayed molecular structures
and rule-transfer examples are extracted from released source tables.
"""
    (output / "README.md").write_text(readme, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    data = load_inputs(args.work)
    audit = validate(data)
    figure2 = build_figure2(data, args.output)
    figure3 = build_figure3(data, args.output)
    write_release(
        data,
        audit,
        args.output,
        Path(__file__).resolve(),
        figure2 + figure3,
    )
    print(f"Built Figures 2 and 3 V6 at {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
