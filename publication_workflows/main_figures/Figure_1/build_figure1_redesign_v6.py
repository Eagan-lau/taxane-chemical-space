#!/usr/bin/env python3
"""Render the publication-grade Figure 1 from frozen evidence tables."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageChops
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


DEFAULT_WORK = Path(".")

COLORS = {
    "ink": "#20262e",
    "gray": "#7f8a96",
    "gray_light": "#eef1f4",
    "grid": "#dce2e7",
    "blue": "#347dbb",
    "blue_light": "#a9c9e3",
    "teal": "#2f8f7f",
    "teal_light": "#b9d9d2",
    "orange": "#dc7b1f",
    "orange_light": "#f3d4b7",
    "gold": "#d9a51f",
    "red": "#c64d4d",
    "white": "#ffffff",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.05) -> None:
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


def quiet_axes(ax, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(
            axis=grid_axis,
            color=COLORS["grid"],
            linewidth=0.55,
            alpha=0.9,
        )
        ax.set_axisbelow(True)


def add_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    facecolor: str,
    edgecolor: str,
    title: str,
    subtitle: str = "",
    title_color: str | None = None,
    linewidth: float = 1.0,
    title_size: float = 7.2,
    subtitle_size: float = 5.7,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height * 0.61,
        title,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=title_size,
        fontweight="bold",
        color=title_color or COLORS["ink"],
    )
    if subtitle:
        ax.text(
            x + width / 2,
            y + height * 0.27,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=subtitle_size,
            color=COLORS["gray"],
        )


def add_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = COLORS["gray"],
    mutation_scale: float = 9,
    linewidth: float = 1.0,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=linewidth,
            color=color,
            transform=ax.transAxes,
            clip_on=False,
        )
    )


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.3,
            "axes.titlesize": 9.0,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.4,
            "legend.fontsize": 5.8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.7,
        }
    )


def load_frozen_data(work: Path) -> dict[str, object]:
    v5 = (
        work
        / "release/manuscript_v5_argument_driven_figures_20260730/source_data"
    )
    provenance = read_tsv(v5 / "Figure_1_B_source_data.tsv")
    stage = read_tsv(v5 / "Figure_1_A_source_data.tsv")
    tier_qc = read_tsv(v5 / "Figure_1_C_source_data.tsv")
    behavior = read_tsv(v5 / "Figure_1_D-F_tier_behavior_source_data.tsv")
    overlap = read_tsv(v5 / "Figure_1_E_source_data.tsv")

    provenance_roles = {
        "TaxolKnownPathway_Curated": "Curated domain anchor",
        "KEGG": "Biochemical resource",
        "Rhea": "Biochemical resource",
        "BioNaviNP_BioChem": "Biochemical resource",
        "MetaNetX": "Integrated reaction resource",
        "RetroRules": "Integrated rule resource",
        "BioNaviNP_USPTO_NPL": "Reaction-chemistry provenance",
    }
    display_names = {
        "TaxolKnownPathway_Curated": "Curated Taxol pathway",
        "BioNaviNP_BioChem": "BioNavi-NP BioChem",
        "BioNaviNP_USPTO_NPL": "BioNavi-NP USPTO_NPL",
    }
    provenance = provenance.copy()
    provenance["provenance_role"] = provenance["source_database"].map(
        provenance_roles
    )
    provenance["display_name"] = provenance["source_database"].map(
        display_names
    ).fillna(provenance["source_database"])

    stage_counts = dict(zip(stage["stage"], stage["count"]))
    tier_qc = tier_qc.set_index("tier")
    behavior = behavior.set_index("tier")
    pipeline_rows = [
        {
            "stage": "Normalized reactions",
            "tier": "all",
            "count": int(stage_counts["Normalized reactions"]),
            "unit": "reaction records",
            "role": "evidence capture",
        },
        {
            "stage": "Deduplicated templates / anchors",
            "tier": "all",
            "count": int(stage_counts["Deduplicated templates / anchors"]),
            "unit": "templates or anchors",
            "role": "build-level deduplication",
        },
        {
            "stage": "Generalized reaction SMARTS",
            "tier": "all",
            "count": int(stage_counts["Generalized reaction SMARTS"]),
            "unit": "SMARTS rows",
            "role": "mapped-centre abstraction and release QC",
        },
    ]
    tier_roles = {
        "T1": "highest-evidence biochemical; primary",
        "T2": "extended biochemical; sensitivity",
        "T3": "exploratory reaction chemistry",
    }
    for tier in ("T1", "T2", "T3"):
        pipeline_rows.extend(
            [
                {
                    "stage": "Exclusive release",
                    "tier": tier,
                    "count": int(tier_qc.loc[tier, "input_release_rows"]),
                    "unit": "SMARTS rows",
                    "role": tier_roles[tier],
                },
                {
                    "stage": "Grammar-QC eligible",
                    "tier": tier,
                    "count": int(
                        tier_qc.loc[tier, "grammar_qc_eligible_rows"]
                    ),
                    "unit": "SMARTS rows",
                    "role": tier_roles[tier],
                },
                {
                    "stage": "Initial executable representatives",
                    "tier": tier,
                    "count": int(
                        tier_qc.loc[
                            tier, "initial_executable_representatives"
                        ]
                    ),
                    "unit": "executable representatives",
                    "role": tier_roles[tier],
                },
                {
                    "stage": "Compiled comparison grammar",
                    "tier": tier,
                    "count": int(behavior.loc[tier, "compiled_rules"]),
                    "unit": "executable productions",
                    "role": tier_roles[tier],
                },
            ]
        )
    pipeline = pd.DataFrame(pipeline_rows)

    grammar_path = (
        work
        / "inputs/G0_G3_primary_release/02_primary_grammar/"
        "02d_T1_final_primary_grammar/taxane_reaction_grammar.primary.tsv"
    )
    grammar_columns = [
        "smarts_rule_id",
        "reaction_smarts",
        "reaction_type",
        "reaction_subtype",
        "normalized_direction",
        "direction_qc_status",
        "evidence_layer_best",
        "template_sources",
        "consensus_qc_status",
        "final_rule_confidence",
        "exclusive_release_tier",
        "grammar_provenance_scope",
    ]
    grammar = read_tsv(grammar_path, usecols=grammar_columns)
    example_ids = [
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000001",
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000003",
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000009",
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000004",
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000005",
        "SMRT_TAXANE_DOMAIN_CONSENSUS_000000013",
    ]
    examples = (
        grammar[grammar["smarts_rule_id"].isin(example_ids)]
        .set_index("smarts_rule_id")
        .loc[example_ids]
        .reset_index()
    )
    if len(examples) != len(example_ids):
        raise RuntimeError("The frozen T1 grammar examples were not uniquely found.")
    examples["display_name"] = [
        "C-H hydroxylation",
        "O-acetyl transfer",
        "Alcohol oxidation",
        "O-deacetylation",
        "O-benzoyl transfer",
        "N-benzoyl transfer",
    ]
    examples["local_edit"] = [
        "C-H  ->  C-OH",
        "O-H  ->  O-acetyl",
        "C-OH  ->  C=O",
        "O-acetyl  ->  O-H",
        "O-H  ->  O-benzoyl",
        "N-H  ->  N-benzoyl",
    ]

    pathway_path = work / "reaction_databases/taxol_pathway.csv"
    pathway = pd.read_csv(pathway_path, encoding="utf-8-sig")
    # Frozen, one-based source rows selected only to illustrate the three
    # corresponding local edits. No pathway reaction is recalculated here.
    pathway_rows = [25, 24, 17, 16, 11, 21]
    pathway_examples = []
    for (_, rule), row_number in zip(examples.iterrows(), pathway_rows):
        reaction = pathway.iloc[row_number - 1]
        pathway_examples.append(
            {
                "display_name": rule["display_name"],
                "local_edit": rule["local_edit"],
                "smarts_rule_id": rule["smarts_rule_id"],
                "reaction_smarts": rule["reaction_smarts"],
                "reaction_type": rule["reaction_type"],
                "pathway_row_number": row_number,
                "enzyme": reaction["Enzyme"],
                "ec": reaction["EC"],
                "substrate_smiles": reaction["Substrate"],
                "product_smiles": reaction["Product"],
            }
        )
    pathway_examples = pd.DataFrame(pathway_examples)

    return {
        "provenance": provenance,
        "pipeline": pipeline,
        "tier_qc": tier_qc.reset_index(),
        "behavior": behavior.reset_index(),
        "overlap": overlap,
        "grammar_examples": pathway_examples,
        "grammar_source": grammar_path,
        "authoritative_sources": {
            "provenance": v5 / "Figure_1_B_source_data.tsv",
            "stage": v5 / "Figure_1_A_source_data.tsv",
            "tier_qc": v5 / "Figure_1_C_source_data.tsv",
            "behavior": v5 / "Figure_1_D-F_tier_behavior_source_data.tsv",
            "overlap": v5 / "Figure_1_E_source_data.tsv",
            "grammar": grammar_path,
            "taxol_pathway": pathway_path,
        },
    }


def validate(data: dict[str, object]) -> pd.DataFrame:
    provenance = data["provenance"]
    pipeline = data["pipeline"]
    tier_qc = data["tier_qc"].set_index("tier")
    behavior = data["behavior"].set_index("tier")
    overlap = data["overlap"]
    grammar_examples = data["grammar_examples"]

    checks = [
        ("normalized_reactions", 630280, int(
            pipeline.loc[
                pipeline["stage"] == "Normalized reactions", "count"
            ].iloc[0]
        )),
        ("generalized_reaction_smarts", 353524, int(
            pipeline.loc[
                pipeline["stage"] == "Generalized reaction SMARTS", "count"
            ].iloc[0]
        )),
        ("exclusive_partition_sum", 353524, int(
            pipeline.loc[
                pipeline["stage"] == "Exclusive release", "count"
            ].sum()
        )),
        ("T1_compiled_rules", 74, int(behavior.loc["T1", "compiled_rules"])),
        ("T2_compiled_rules", 17, int(behavior.loc["T2", "compiled_rules"])),
        ("T3_compiled_rules", 1568, int(behavior.loc["T3", "compiled_rules"])),
        (
            "T1_initial_executable_representatives",
            11214,
            int(tier_qc.loc["T1", "initial_executable_representatives"]),
        ),
        (
            "T2_initial_executable_representatives",
            6593,
            int(tier_qc.loc["T2", "initial_executable_representatives"]),
        ),
        (
            "T3_initial_executable_representatives",
            2647,
            int(tier_qc.loc["T3", "initial_executable_representatives"]),
        ),
        (
            "T1_unique_G1",
            15801,
            int(behavior.loc["T1", "unique_G1_structures"]),
        ),
        (
            "T2_unique_G1",
            3069,
            int(behavior.loc["T2", "unique_G1_structures"]),
        ),
        (
            "T3_unique_G1",
            2973,
            int(behavior.loc["T3", "unique_G1_structures"]),
        ),
        ("provenance_records", 630280, int(
            provenance["normalized_source_reaction_records"].sum()
        )),
        ("overlap_pair_count", 3, int(len(overlap))),
        ("grammar_example_count", 6, int(len(grammar_examples))),
    ]
    rows = []
    for claim, expected, observed in checks:
        rows.append(
            {
                "claim_id": claim,
                "expected_frozen_value": expected,
                "observed_value": observed,
                "status": "PASS" if expected == observed else "FAIL",
            }
        )
    audit = pd.DataFrame(rows)
    if (audit["status"] != "PASS").any():
        raise RuntimeError("Figure 1 numerical audit failed.")
    return audit


def draw_panel_a(ax, provenance: pd.DataFrame) -> None:
    role_colors = {
        "Curated domain anchor": COLORS["red"],
        "Biochemical resource": COLORS["blue"],
        "Integrated reaction resource": COLORS["teal"],
        "Integrated rule resource": COLORS["teal"],
        "Reaction-chemistry provenance": COLORS["orange"],
    }
    ordered = provenance.sort_values(
        "normalized_source_reaction_records", ascending=True
    )
    y = np.arange(len(ordered))
    values = ordered["normalized_source_reaction_records"].to_numpy()
    colors = [role_colors[value] for value in ordered["provenance_role"]]
    ax.barh(
        y,
        values,
        color=colors,
        edgecolor=COLORS["ink"],
        linewidth=0.25,
        height=0.72,
    )
    ax.set_yticks(y, ordered["display_name"])
    ax.set_xscale("log")
    ax.set_xlim(20, 7.5e5)
    ax.set_xlabel("Normalized source reaction records (log scale)")
    ax.set_title("Database provenance and evidence roles", pad=8)
    for yi, value in zip(y, values):
        ax.text(
            value * 1.06,
            yi,
            f"{int(value):,}",
            va="center",
            ha="left",
            fontsize=5.8,
            color=COLORS["ink"],
        )
    handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor=color,
            markeredgecolor=COLORS["ink"],
            markeredgewidth=0.25,
            markersize=6,
            label=label,
        )
        for label, color in [
            ("curated domain anchor", COLORS["red"]),
            ("biochemical resource", COLORS["blue"]),
            ("integrated resource", COLORS["teal"]),
            ("reaction chemistry", COLORS["orange"]),
        ]
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        ncol=2,
        loc="lower right",
        fontsize=5.2,
        columnspacing=0.8,
        handletextpad=0.4,
    )
    ax.text(
        0.0,
        -0.24,
        "Colour denotes provenance role, not evidence confidence.",
        transform=ax.transAxes,
        fontsize=5.5,
        color=COLORS["gray"],
    )
    quiet_axes(ax, "x")
    panel_label(ax, "A")


def draw_panel_b(
    ax,
    pipeline: pd.DataFrame,
    tier_qc: pd.DataFrame,
    behavior: pd.DataFrame,
) -> None:
    ax.set_axis_off()
    ax.set_title(
        "Evidence normalization, tier quality gates and grammar selection",
        pad=8,
    )
    panel_label(ax, "B", x=-0.035, y=1.04)

    common = pipeline[pipeline["tier"] == "all"].set_index("stage")
    common_boxes = [
        (
            0.02,
            0.20,
            "Normalized reactions",
            COLORS["gray_light"],
            COLORS["gray"],
        ),
        (
            0.28,
            0.20,
            "Deduplicated templates / anchors",
            COLORS["blue_light"],
            COLORS["blue"],
        ),
        (
            0.54,
            0.20,
            "Generalized reaction SMARTS",
            COLORS["teal_light"],
            COLORS["teal"],
        ),
    ]
    for x, width, stage, face, edge in common_boxes:
        count = int(common.loc[stage, "count"])
        add_box(
            ax,
            x,
            0.79,
            width,
            0.14,
            facecolor=face,
            edgecolor=edge,
            title=f"{count:,}",
            subtitle=stage,
            title_size=8.5,
            subtitle_size=5.0,
        )
    add_arrow(ax, (0.225, 0.86), (0.272, 0.86))
    add_arrow(ax, (0.485, 0.86), (0.532, 0.86))
    ax.text(
        0.25,
        0.938,
        "principal-pair projection\nand deduplication",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=4.8,
        color=COLORS["gray"],
    )
    ax.text(
        0.51,
        0.938,
        "atom mapping and\nreaction-centre abstraction",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=4.8,
        color=COLORS["gray"],
    )
    add_arrow(ax, (0.745, 0.86), (0.792, 0.86))
    partition = tier_qc.set_index("tier")
    partition_text = (
        f"T1 {int(partition.loc['T1', 'input_release_rows']):,}  |  "
        f"T2 {int(partition.loc['T2', 'input_release_rows']):,}  |  "
        f"T3 {int(partition.loc['T3', 'input_release_rows']):,}"
    )
    add_box(
        ax,
        0.80,
        0.79,
        0.18,
        0.14,
        facecolor=COLORS["white"],
        edgecolor=COLORS["ink"],
        title="Exclusive tier assignment",
        subtitle=partition_text,
        title_size=6.6,
        subtitle_size=4.7,
    )
    ax.text(
        0.77,
        0.938,
        "evidence stratification",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=4.8,
        color=COLORS["gray"],
    )
    add_arrow(
        ax,
        (0.89, 0.785),
        (0.89, 0.705),
        color=COLORS["gray"],
        mutation_scale=8,
    )

    quality_ax = ax.inset_axes([0.02, 0.08, 0.53, 0.57])
    tier_order = ["T1", "T2", "T3"]
    quality = tier_qc.set_index("tier").loc[tier_order]
    x = np.arange(len(tier_order))
    width = 0.24
    quality_columns = [
        ("input_release_rows", "release rows", COLORS["gray"]),
        ("grammar_qc_eligible_rows", "QC eligible", COLORS["blue"]),
        (
            "initial_executable_representatives",
            "initial executable representatives",
            COLORS["orange"],
        ),
    ]
    for index, (column, label, color) in enumerate(quality_columns):
        quality_ax.bar(
            x + (index - 1) * width,
            quality[column],
            width=width,
            color=color,
            edgecolor=COLORS["ink"],
            linewidth=0.25,
            label=label,
        )
    quality_ax.set_yscale("log")
    quality_ax.set_ylim(1.8e3, 3.0e5)
    quality_ax.set_xticks(x, tier_order)
    quality_ax.set_ylabel("Records (log scale)", labelpad=2)
    quality_ax.set_title(
        "Tier-specific evidence and executability gates",
        fontsize=7.0,
        pad=5,
    )
    quality_ax.legend(
        frameon=False,
        fontsize=4.7,
        loc="upper right",
        handlelength=1.2,
    )
    quiet_axes(quality_ax, "y")

    behavior_index = behavior.set_index("tier")
    role_specs = [
        (
            "T1",
            COLORS["blue"],
            COLORS["blue_light"],
            "primary grammar",
            "highest-evidence biochemical",
        ),
        (
            "T2",
            COLORS["teal"],
            COLORS["teal_light"],
            "sensitivity grammar",
            "extended biochemical",
        ),
        (
            "T3",
            COLORS["orange"],
            COLORS["orange_light"],
            "exploratory grammar",
            "reaction-chemistry exploration",
        ),
    ]
    ax.text(
        0.77,
        0.65,
        "Tier-specific semantic selection\nand compilation",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.7,
        fontweight="bold",
        color=COLORS["ink"],
    )
    for index, (tier, edge, fill, role, evidence) in enumerate(role_specs):
        y = 0.46 - index * 0.16
        count = int(behavior_index.loc[tier, "compiled_rules"])
        add_arrow(
            ax,
            (0.56, y + 0.055),
            (0.605, y + 0.055),
            color=edge,
            mutation_scale=7,
            linewidth=0.8,
        )
        add_box(
            ax,
            0.61,
            y,
            0.34,
            0.11,
            facecolor=fill,
            edgecolor=edge,
            title=f"{tier}  |  {count:,} executable productions",
            subtitle=f"{role}; {evidence}",
            title_color=edge,
            linewidth=1.25 if tier == "T1" else 0.9,
            title_size=6.6,
            subtitle_size=4.8,
        )
    ax.text(
        0.77,
        0.03,
        "Final boxes report the frozen tier-comparison grammars. "
        "T1 includes separately traceable curated taxane-domain rules.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=4.8,
        color=COLORS["gray"],
    )


def draw_panel_c(ax, examples: pd.DataFrame) -> None:
    ax.set_axis_off()
    ax.set_title(
        "How reaction grammar encodes transferable molecular edits",
        pad=8,
    )
    panel_label(ax, "C", x=-0.025, y=1.04)
    ax.text(
        0.5,
        0.91,
        "A production couples a matched local molecular pattern to a "
        "directional product rewrite; the surrounding taxane context is "
        "retained.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.2,
        color=COLORS["ink"],
    )

    card_specs = [
        (0.02, COLORS["blue"], COLORS["blue_light"]),
        (0.345, COLORS["orange"], COLORS["orange_light"]),
        (0.67, COLORS["teal"], COLORS["teal_light"]),
    ]
    for (_, row), (x, edge, fill) in zip(examples.iterrows(), card_specs):
        width = 0.305
        card = FancyBboxPatch(
            (x, 0.13),
            width,
            0.68,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            facecolor=COLORS["white"],
            edgecolor=COLORS["grid"],
            linewidth=0.9,
            transform=ax.transAxes,
            clip_on=False,
        )
        ax.add_patch(card)
        ax.add_patch(
            FancyBboxPatch(
                (x + 0.012, 0.705),
                width - 0.024,
                0.08,
                boxstyle="round,pad=0.004,rounding_size=0.009",
                facecolor=fill,
                edgecolor=edge,
                linewidth=0.8,
                transform=ax.transAxes,
                clip_on=False,
            )
        )
        ax.text(
            x + width / 2,
            0.745,
            str(row["display_name"]),
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7.1,
            fontweight="bold",
            color=edge,
        )
        ax.text(
            x + 0.018,
            0.665,
            "Schematic molecular context",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=4.8,
            color=COLORS["gray"],
            fontweight="bold",
        )
        ax.text(
            x + 0.078,
            0.575,
            str(row["substrate_context"]),
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8.8,
            color=COLORS["ink"],
        )
        add_arrow(
            ax,
            (x + 0.135, 0.575),
            (x + 0.177, 0.575),
            color=edge,
            mutation_scale=9,
            linewidth=1.0,
        )
        ax.text(
            x + 0.236,
            0.575,
            str(row["product_context"]),
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8.8,
            color=edge,
        )
        ax.text(
            x + width / 2,
            0.485,
            f"Local edit: {row['local_edit']}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=6.0,
            fontweight="bold",
            color=edge,
        )
        ax.plot(
            [x + 0.018, x + width - 0.018],
            [0.435, 0.435],
            transform=ax.transAxes,
            color=COLORS["grid"],
            linewidth=0.7,
        )
        ax.text(
            x + 0.018,
            0.395,
            "Executable reaction-SMARTS production",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=4.8,
            color=COLORS["gray"],
            fontweight="bold",
        )
        reaction_smarts = str(row["reaction_smarts"]).replace(
            ">>", "\n>> "
        )
        ax.text(
            x + 0.018,
            0.31,
            reaction_smarts,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=5.2,
            family="DejaVu Sans Mono",
            color=COLORS["ink"],
            linespacing=1.25,
        )
        ax.text(
            x + 0.018,
            0.19,
            f"T1 | direction-aware | {row['smarts_rule_id']}",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=4.5,
            color=COLORS["gray"],
        )
    ax.text(
        0.5,
        0.045,
        "Reaction grammar is a transferable local rewrite, not a memorized "
        "whole-molecule reaction. Each example is a real frozen T1 rule.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.6,
        fontweight="bold",
        color=COLORS["ink"],
    )


def draw_panel_b_v4(
    ax,
    pipeline: pd.DataFrame,
    tier_qc: pd.DataFrame,
    behavior: pd.DataFrame,
) -> None:
    """Draw a compact construction flow and tier-retention matrix."""
    ax.set_axis_off()
    ax.set_title(
        "From heterogeneous reaction evidence to evidence-stratified grammars",
        pad=8,
    )
    panel_label(ax, "B", x=-0.035, y=1.04)

    common = pipeline[pipeline["tier"] == "all"].set_index("stage")
    stages = [
        ("Normalized\nreactions", "Normalized reactions", COLORS["gray"]),
        (
            "Deduplicated\ntemplates",
            "Deduplicated templates / anchors",
            COLORS["blue"],
        ),
        (
            "Generalized\nreaction SMARTS",
            "Generalized reaction SMARTS",
            COLORS["teal"],
        ),
    ]
    x_positions = [0.015, 0.235, 0.455]
    block_width = 0.185
    for (label, stage, color), x in zip(stages, x_positions):
        count = int(common.loc[stage, "count"])
        ax.add_patch(
            Rectangle(
                (x, 0.735),
                block_width,
                0.19,
                transform=ax.transAxes,
                facecolor=COLORS["white"],
                edgecolor=COLORS["grid"],
                linewidth=0.8,
            )
        )
        ax.add_patch(
            Rectangle(
                (x, 0.735),
                0.008,
                0.19,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="none",
            )
        )
        ax.text(
            x + 0.018,
            0.855,
            f"{count:,}",
            transform=ax.transAxes,
            fontsize=10.0,
            fontweight="bold",
            color=COLORS["ink"],
            ha="left",
            va="center",
        )
        ax.text(
            x + 0.018,
            0.785,
            label,
            transform=ax.transAxes,
            fontsize=5.8,
            color=COLORS["gray"],
            ha="left",
            va="center",
            linespacing=1.15,
        )
    add_arrow(ax, (0.202, 0.83), (0.228, 0.83), mutation_scale=7)
    add_arrow(ax, (0.422, 0.83), (0.448, 0.83), mutation_scale=7)
    add_arrow(ax, (0.642, 0.83), (0.668, 0.83), mutation_scale=7)

    partition = tier_qc.set_index("tier")
    ax.add_patch(
        Rectangle(
            (0.675, 0.735),
            0.31,
            0.19,
            transform=ax.transAxes,
            facecolor="#fafbfd",
            edgecolor=COLORS["grid"],
            linewidth=0.8,
        )
    )
    ax.text(
        0.69,
        0.883,
        "Exclusive evidence tiers",
        transform=ax.transAxes,
        fontsize=6.2,
        fontweight="bold",
        ha="left",
        va="center",
    )
    tier_colors = {
        "T1": COLORS["blue"],
        "T2": COLORS["teal"],
        "T3": COLORS["orange"],
    }
    for i, tier in enumerate(("T1", "T2", "T3")):
        x = 0.69 + i * 0.095
        ax.text(
            x,
            0.814,
            tier,
            transform=ax.transAxes,
            fontsize=6.0,
            fontweight="bold",
            color=tier_colors[tier],
            ha="left",
            va="center",
        )
        ax.text(
            x,
            0.772,
            f"{int(partition.loc[tier, 'input_release_rows']):,}",
            transform=ax.transAxes,
            fontsize=5.5,
            color=COLORS["ink"],
            ha="left",
            va="center",
        )

    ax.text(
        0.02,
        0.64,
        "Tier retention through grammar quality control",
        transform=ax.transAxes,
        fontsize=6.5,
        fontweight="bold",
        ha="left",
        va="center",
    )
    headers = [
        ("EVIDENCE", 0.02, "left"),
        ("EXCLUSIVE RELEASE", 0.265, "center"),
        ("QC ELIGIBLE", 0.455, "center"),
        ("EXECUTABLE", 0.635, "center"),
        ("FINAL GRAMMAR", 0.815, "center"),
        ("ROLE", 0.985, "right"),
    ]
    for label, x, align in headers:
        ax.text(
            x,
            0.57,
            label,
            transform=ax.transAxes,
            fontsize=4.7,
            fontweight="bold",
            color=COLORS["gray"],
            ha=align,
            va="center",
        )
    ax.plot(
        [0.02, 0.985],
        [0.535, 0.535],
        transform=ax.transAxes,
        color=COLORS["grid"],
        linewidth=0.8,
    )

    quality = tier_qc.set_index("tier")
    final = behavior.set_index("tier")
    roles = {"T1": "PRIMARY", "T2": "SENSITIVITY", "T3": "EXPLORATORY"}
    y_positions = {"T1": 0.435, "T2": 0.29, "T3": 0.145}
    for tier in ("T1", "T2", "T3"):
        y = y_positions[tier]
        color = tier_colors[tier]
        if tier == "T1":
            ax.add_patch(
                Rectangle(
                    (0.012, y - 0.058),
                    0.976,
                    0.116,
                    transform=ax.transAxes,
                    facecolor="#f1f7fc",
                    edgecolor="none",
                )
            )
        ax.add_patch(
            Rectangle(
                (0.02, y - 0.052),
                0.008,
                0.104,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="none",
            )
        )
        ax.text(
            0.045,
            y,
            tier,
            transform=ax.transAxes,
            fontsize=8.2,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
        )
        values = [
            int(quality.loc[tier, "input_release_rows"]),
            int(quality.loc[tier, "grammar_qc_eligible_rows"]),
            int(quality.loc[tier, "initial_executable_representatives"]),
            int(final.loc[tier, "compiled_rules"]),
        ]
        for x, value in zip((0.265, 0.455, 0.635, 0.815), values):
            ax.text(
                x,
                y,
                f"{value:,}",
                transform=ax.transAxes,
                fontsize=6.8 if x == 0.815 else 6.2,
                fontweight="bold" if x == 0.815 else "normal",
                color=color if x == 0.815 else COLORS["ink"],
                ha="center",
                va="center",
            )
        ax.text(
            0.985,
            y,
            roles[tier],
            transform=ax.transAxes,
            fontsize=5.0,
            fontweight="bold",
            color=color,
            ha="right",
            va="center",
        )
        if tier != "T3":
            ax.plot(
                [0.02, 0.985],
                [y - 0.073, y - 0.073],
                transform=ax.transAxes,
                color=COLORS["grid"],
                linewidth=0.55,
            )
    ax.text(
        0.02,
        0.025,
        "T1 contains the highest-evidence biochemical layer and traceable "
        "curated taxane-domain productions.",
        transform=ax.transAxes,
        fontsize=5.0,
        color=COLORS["gray"],
        ha="left",
        va="bottom",
    )


def draw_panel_b_v5(
    ax,
    pipeline: pd.DataFrame,
    tier_qc: pd.DataFrame,
    behavior: pd.DataFrame,
) -> None:
    """Draw the evidence hierarchy as a branched, contracting flow."""
    ax.set_axis_off()
    ax.set_title(
        "Evidence stratification and contraction to executable grammars",
        pad=8,
    )
    panel_label(ax, "B", x=-0.035, y=1.04)

    common = pipeline[pipeline["tier"] == "all"].set_index("stage")
    stages = [
        (
            0.02,
            "Normalized reactions",
            int(common.loc["Normalized reactions", "count"]),
            COLORS["gray"],
        ),
        (
            0.23,
            "Deduplicated templates",
            int(common.loc["Deduplicated templates / anchors", "count"]),
            COLORS["blue"],
        ),
        (
            0.44,
            "Generalized reaction SMARTS",
            int(common.loc["Generalized reaction SMARTS", "count"]),
            COLORS["teal"],
        ),
    ]
    for x, label, count, color in stages:
        ax.add_patch(
            Rectangle(
                (x, 0.785),
                0.17,
                0.15,
                transform=ax.transAxes,
                facecolor=COLORS["white"],
                edgecolor=COLORS["grid"],
                linewidth=0.8,
            )
        )
        ax.add_patch(
            Rectangle(
                (x, 0.785),
                0.007,
                0.15,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="none",
            )
        )
        ax.text(
            x + 0.016,
            0.877,
            f"{count:,}",
            transform=ax.transAxes,
            fontsize=9.0,
            fontweight="bold",
            color=COLORS["ink"],
            ha="left",
            va="center",
        )
        ax.text(
            x + 0.016,
            0.818,
            label,
            transform=ax.transAxes,
            fontsize=5.1,
            color=COLORS["gray"],
            ha="left",
            va="center",
        )
    add_arrow(ax, (0.192, 0.86), (0.223, 0.86), mutation_scale=7)
    add_arrow(ax, (0.402, 0.86), (0.433, 0.86), mutation_scale=7)

    hub_x, hub_y = 0.66, 0.86
    add_arrow(ax, (0.612, 0.86), (hub_x - 0.012, hub_y), mutation_scale=7)
    ax.scatter(
        [hub_x],
        [hub_y],
        s=66,
        marker="D",
        color=COLORS["ink"],
        transform=ax.transAxes,
        clip_on=False,
        zorder=5,
    )
    ax.text(
        hub_x,
        0.955,
        "evidence\nstratification",
        transform=ax.transAxes,
        fontsize=4.6,
        color=COLORS["gray"],
        ha="center",
        va="bottom",
        linespacing=1.05,
    )

    column_x = [0.72, 0.81, 0.90, 0.975]
    column_labels = ["EXCLUSIVE", "QC", "EXECUTABLE", "COMPILED"]
    for x, label in zip(column_x, column_labels):
        ax.text(
            x,
            0.685,
            label,
            transform=ax.transAxes,
            fontsize=4.6,
            fontweight="bold",
            color=COLORS["gray"],
            ha="center",
            va="center",
        )

    quality = tier_qc.set_index("tier")
    final = behavior.set_index("tier")
    tier_specs = [
        ("T1", 0.52, COLORS["blue"], "PRIMARY"),
        ("T2", 0.315, COLORS["teal"], "SENSITIVITY"),
        ("T3", 0.11, COLORS["orange"], "EXPLORATORY"),
    ]

    def node_size(value: int) -> float:
        return 165 + 110 * max(0.0, np.log10(max(value, 1)) - 1.0)

    for tier, y, color, role in tier_specs:
        values = [
            int(quality.loc[tier, "input_release_rows"]),
            int(quality.loc[tier, "grammar_qc_eligible_rows"]),
            int(quality.loc[tier, "initial_executable_representatives"]),
            int(final.loc[tier, "compiled_rules"]),
        ]
        branch = FancyArrowPatch(
            (hub_x, hub_y - 0.012),
            (column_x[0] - 0.012, y),
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.10",
            mutation_scale=7,
            linewidth=1.1,
            color=color,
            alpha=0.72,
            transform=ax.transAxes,
            clip_on=False,
            zorder=1,
        )
        ax.add_patch(branch)
        for index in range(3):
            linewidth = 1.0 + 0.48 * (
                np.log10(max(values[index + 1], 1)) - 1.0
            )
            connector = FancyArrowPatch(
                (column_x[index] + 0.012, y),
                (column_x[index + 1] - 0.012, y),
                arrowstyle="-|>",
                mutation_scale=7,
                linewidth=max(0.9, linewidth),
                color=color,
                alpha=0.58,
                transform=ax.transAxes,
                clip_on=False,
                zorder=1,
            )
            ax.add_patch(connector)

        for index, (x, value) in enumerate(zip(column_x, values)):
            alpha = [0.24, 0.38, 0.62, 0.95][index]
            ax.scatter(
                [x],
                [y],
                s=node_size(value),
                color=color,
                alpha=alpha,
                edgecolors=color,
                linewidths=0.8,
                transform=ax.transAxes,
                clip_on=False,
                zorder=3,
            )
            ax.text(
                x,
                y,
                f"{value:,}",
                transform=ax.transAxes,
                fontsize=5.2 if value < 100000 else 4.8,
                fontweight="bold" if index == 3 else "normal",
                color=COLORS["white"] if index >= 2 else COLORS["ink"],
                ha="center",
                va="center",
                zorder=6,
            )
        ax.text(
            0.675,
            y,
            tier,
            transform=ax.transAxes,
            fontsize=7.2,
            fontweight="bold",
            color=color,
            ha="right",
            va="center",
        )
        ax.text(
            0.998,
            y - 0.073,
            role,
            transform=ax.transAxes,
            fontsize=4.7,
            fontweight="bold",
            color=color,
            ha="right",
            va="center",
        )
def reaction_center_highlights(
    substrate: Chem.Mol,
    product: Chem.Mol,
) -> tuple[list[int], list[int], list[int], list[int], dict[int, int]]:
    """Identify display-only changed atoms and bonds from a maximum common graph."""
    result = rdFMCS.FindMCS(
        [substrate, product],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=10,
    )
    pattern = Chem.MolFromSmarts(result.smartsString)
    substrate_match = substrate.GetSubstructMatch(pattern)
    product_match = product.GetSubstructMatch(pattern)
    atom_map = dict(zip(substrate_match, product_match))
    substrate_atoms = set(range(substrate.GetNumAtoms())) - set(substrate_match)
    product_atoms = set(range(product.GetNumAtoms())) - set(product_match)
    substrate_bonds: set[int] = set()
    product_bonds: set[int] = set()

    for atom_index in list(substrate_atoms):
        substrate_atoms.update(
            neighbor.GetIdx()
            for neighbor in substrate.GetAtomWithIdx(atom_index).GetNeighbors()
        )
    for atom_index in list(product_atoms):
        product_atoms.update(
            neighbor.GetIdx()
            for neighbor in product.GetAtomWithIdx(atom_index).GetNeighbors()
        )

    for bond in substrate.GetBonds():
        left = bond.GetBeginAtomIdx()
        right = bond.GetEndAtomIdx()
        if left not in atom_map or right not in atom_map:
            continue
        product_bond = product.GetBondBetweenAtoms(
            atom_map[left],
            atom_map[right],
        )
        if product_bond is None or product_bond.GetBondType() != bond.GetBondType():
            substrate_atoms.update((left, right))
            product_atoms.update((atom_map[left], atom_map[right]))
            substrate_bonds.add(bond.GetIdx())
            if product_bond is not None:
                product_bonds.add(product_bond.GetIdx())

    for substrate_index, product_index in atom_map.items():
        left = substrate.GetAtomWithIdx(substrate_index)
        right = product.GetAtomWithIdx(product_index)
        if (
            left.GetFormalCharge() != right.GetFormalCharge()
            or left.GetTotalNumHs() != right.GetTotalNumHs()
        ):
            substrate_atoms.add(substrate_index)
            product_atoms.add(product_index)

    return (
        sorted(substrate_atoms),
        sorted(product_atoms),
        sorted(substrate_bonds),
        sorted(product_bonds),
        atom_map,
    )


def render_reaction_pair(
    substrate: Chem.Mol,
    product: Chem.Mol,
    highlights: tuple[
        list[int],
        list[int],
        list[int],
        list[int],
        dict[int, int],
    ],
    *,
    width: int = 1800,
    height: int = 560,
) -> np.ndarray:
    """Render an aligned substrate-product pair at a common scale."""
    rdDepictor.Compute2DCoords(substrate)
    rdDepictor.GenerateDepictionMatching2DStructure(
        product,
        substrate,
        [(left, right) for left, right in highlights[4].items()],
    )
    panel_width = width // 2
    drawer = rdMolDraw2D.MolDraw2DCairo(
        width,
        height,
        panel_width,
        height,
    )
    options = drawer.drawOptions()
    options.padding = 0.055
    options.bondLineWidth = 1.8
    options.highlightBondWidthMultiplier = 12
    options.highlightRadius = 0.24
    options.addStereoAnnotation = False
    highlight_color = (0.95, 0.55, 0.16)
    drawer.DrawMolecules(
        [substrate, product],
        highlightAtoms=[highlights[0], highlights[1]],
        highlightBonds=[highlights[2], highlights[3]],
        highlightAtomColors=[
            {atom_index: highlight_color for atom_index in highlights[0]},
            {atom_index: highlight_color for atom_index in highlights[1]},
        ],
        highlightBondColors=[
            {bond_index: highlight_color for bond_index in highlights[2]},
            {bond_index: highlight_color for bond_index in highlights[3]},
        ],
    )
    drawer.FinishDrawing()
    image = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")
    background = Image.new("RGB", image.size, "white")
    difference = ImageChops.difference(image, background).convert("L")
    box = difference.point(lambda value: 255 if value > 8 else 0).getbbox()
    if box is not None:
        left, upper, right, lower = box
        pad = 18
        image = image.crop(
            (
                max(0, left - pad),
                max(0, upper - pad),
                min(image.width, right + pad),
                min(image.height, lower + pad),
            )
        )
    return np.asarray(image)


def draw_panel_c_v4(ax, examples: pd.DataFrame) -> None:
    """Use complete curated molecular pairs to define reaction grammar visually."""
    ax.set_axis_off()
    ax.set_title(
        "Reaction grammar links whole-molecule context to transferable local edits",
        pad=8,
    )
    panel_label(ax, "C", x=-0.018, y=1.03)

    ax.text(
        0.21,
        0.968,
        "CURATED SUBSTRATE",
        transform=ax.transAxes,
        fontsize=5.0,
        fontweight="bold",
        color=COLORS["gray"],
        ha="center",
        va="center",
    )
    ax.text(
        0.565,
        0.968,
        "CURATED PRODUCT",
        transform=ax.transAxes,
        fontsize=5.0,
        fontweight="bold",
        color=COLORS["gray"],
        ha="center",
        va="center",
    )
    ax.text(
        0.865,
        0.968,
        "GENERALIZED T1 PRODUCTION",
        transform=ax.transAxes,
        fontsize=5.0,
        fontweight="bold",
        color=COLORS["gray"],
        ha="center",
        va="center",
    )
    row_bases = [0.66, 0.355, 0.05]
    accents = [COLORS["blue"], COLORS["orange"], COLORS["teal"]]
    for row_index, ((_, row), y, accent) in enumerate(
        zip(examples.iterrows(), row_bases, accents)
    ):
        substrate = Chem.MolFromSmiles(str(row["substrate_smiles"]))
        product = Chem.MolFromSmiles(str(row["product_smiles"]))
        if substrate is None or product is None:
            raise ValueError("A curated pathway molecule could not be parsed.")
        highlights = reaction_center_highlights(substrate, product)
        pair_image = render_reaction_pair(
            substrate,
            product,
            highlights,
        )

        if row_index:
            ax.plot(
                [0.012, 0.988],
                [y + 0.295, y + 0.295],
                transform=ax.transAxes,
                color=COLORS["grid"],
                linewidth=0.7,
            )
        ax.add_patch(
            Rectangle(
                (0.012, y + 0.022),
                0.007,
                0.235,
                transform=ax.transAxes,
                facecolor=accent,
                edgecolor="none",
            )
        )
        ax.text(
            0.03,
            y + 0.235,
            str(row["display_name"]),
            transform=ax.transAxes,
            fontsize=7.0,
            fontweight="bold",
            color=accent,
            ha="left",
            va="top",
        )

        pair_ax = ax.inset_axes([0.07, y + 0.015, 0.64, 0.245])
        pair_ax.imshow(pair_image)
        pair_ax.annotate(
            "",
            xy=(0.535, 0.50),
            xytext=(0.465, 0.50),
            xycoords="axes fraction",
            arrowprops={
                "arrowstyle": "-|>",
                "color": accent,
                "linewidth": 1.5,
                "mutation_scale": 11,
            },
        )
        pair_ax.set_axis_off()

        ax.plot(
            [0.735, 0.735],
            [y + 0.03, y + 0.255],
            transform=ax.transAxes,
            color=accent,
            linewidth=1.5,
        )
        ax.text(
            0.755,
            y + 0.225,
            str(row["local_edit"]),
            transform=ax.transAxes,
            fontsize=7.1,
            fontweight="bold",
            color=accent,
            ha="left",
            va="center",
        )
        left, right = str(row["reaction_smarts"]).split(">>", 1)
        ax.text(
            0.755,
            y + 0.148,
            f"{left}\n>> {right}",
            transform=ax.transAxes,
            fontsize=5.3,
            family="DejaVu Sans Mono",
            color=COLORS["ink"],
            ha="left",
            va="center",
            linespacing=1.3,
        )
        ax.text(
            0.755,
            y + 0.058,
            f"T1  |  {row['enzyme']}  |  EC {row['ec']}",
            transform=ax.transAxes,
            fontsize=4.8,
            color=COLORS["gray"],
            ha="left",
            va="center",
        )
    ax.scatter(
        [0.02],
        [0.012],
        s=28,
        color=COLORS["orange"],
        edgecolors="none",
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(
        0.032,
        0.012,
        "highlighted reaction centre",
        transform=ax.transAxes,
        fontsize=4.7,
        color=COLORS["gray"],
        ha="left",
        va="center",
    )


def draw_panel_c_v5(ax, examples: pd.DataFrame) -> None:
    """Render six aligned curated transformations in a molecular gallery."""
    ax.set_axis_off()
    ax.set_title(
        "Whole-molecule examples of transferable T1 reaction grammar",
        pad=8,
    )
    panel_label(ax, "C", x=-0.018, y=1.03)
    ax.scatter(
        [0.805],
        [0.982],
        s=25,
        color=COLORS["orange"],
        edgecolors="none",
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(
        0.816,
        0.982,
        "reaction centre",
        transform=ax.transAxes,
        fontsize=4.7,
        color=COLORS["gray"],
        ha="left",
        va="center",
    )
    ax.text(
        0.985,
        0.982,
        "SUBSTRATE  ->  PRODUCT",
        transform=ax.transAxes,
        fontsize=4.7,
        fontweight="bold",
        color=COLORS["gray"],
        ha="right",
        va="center",
    )

    x_positions = [0.012, 0.508]
    y_positions = [0.655, 0.345, 0.035]
    accents = [
        COLORS["blue"],
        COLORS["orange"],
        COLORS["teal"],
        COLORS["red"],
        COLORS["gold"],
        COLORS["ink"],
    ]
    for index, (_, row) in enumerate(examples.iterrows()):
        column = 0 if index < 3 else 1
        row_index = index if index < 3 else index - 3
        x = x_positions[column]
        y = y_positions[row_index]
        accent = accents[index]

        substrate = Chem.MolFromSmiles(str(row["substrate_smiles"]))
        product = Chem.MolFromSmiles(str(row["product_smiles"]))
        if substrate is None or product is None:
            raise ValueError("A curated pathway molecule could not be parsed.")
        highlights = reaction_center_highlights(substrate, product)
        pair_image = render_reaction_pair(substrate, product, highlights)

        if row_index:
            ax.plot(
                [x, x + 0.476],
                [y + 0.302, y + 0.302],
                transform=ax.transAxes,
                color=COLORS["grid"],
                linewidth=0.65,
            )
        if column:
            ax.plot(
                [0.495, 0.495],
                [0.02, 0.955],
                transform=ax.transAxes,
                color=COLORS["grid"],
                linewidth=0.8,
            )
        ax.add_patch(
            Rectangle(
                (x, y + 0.02),
                0.006,
                0.255,
                transform=ax.transAxes,
                facecolor=accent,
                edgecolor="none",
            )
        )
        ax.text(
            x + 0.016,
            y + 0.265,
            str(row["display_name"]),
            transform=ax.transAxes,
            fontsize=6.8,
            fontweight="bold",
            color=accent,
            ha="left",
            va="top",
        )

        pair_ax = ax.inset_axes([x + 0.012, y + 0.025, 0.320, 0.225])
        pair_ax.imshow(pair_image)
        pair_ax.annotate(
            "",
            xy=(0.535, 0.50),
            xytext=(0.465, 0.50),
            xycoords="axes fraction",
            arrowprops={
                "arrowstyle": "-|>",
                "color": accent,
                "linewidth": 1.3,
                "mutation_scale": 10,
            },
        )
        pair_ax.set_axis_off()

        text_x = x + 0.345
        ax.text(
            text_x,
            y + 0.218,
            str(row["local_edit"]),
            transform=ax.transAxes,
            fontsize=6.3,
            fontweight="bold",
            color=accent,
            ha="left",
            va="center",
        )
        left, right = str(row["reaction_smarts"]).split(">>", 1)
        ax.text(
            text_x,
            y + 0.135,
            f"{left}\n>> {right}",
            transform=ax.transAxes,
            fontsize=4.15,
            family="DejaVu Sans Mono",
            color=COLORS["ink"],
            ha="left",
            va="center",
            linespacing=1.22,
        )
        ax.text(
            text_x,
            y + 0.052,
            f"T1  |  {row['enzyme']}  |  EC {row['ec']}",
            transform=ax.transAxes,
            fontsize=4.25,
            color=COLORS["gray"],
            ha="left",
            va="center",
        )


def draw_panel_d(ax, behavior: pd.DataFrame) -> None:
    colors = [COLORS["blue"], COLORS["teal"], COLORS["gold"]]
    tiers = behavior["tier"].tolist()
    x = np.arange(len(tiers))

    bars = ax.bar(
        x,
        behavior["unique_G1_structures"],
        color=colors,
        edgecolor=COLORS["ink"],
        linewidth=0.4,
        width=0.63,
    )
    ax.set_xticks(x, tiers)
    ax.set_ylabel("Unique one-step structures")
    ax.set_ylim(0, 16600)
    ax.set_title("One-step coverage and product acceptance", pad=7)
    for bar, value in zip(bars, behavior["unique_G1_structures"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 260,
            f"{int(value):,}",
            ha="center",
            va="bottom",
            fontsize=6.3,
            fontweight="bold",
        )
    twin = ax.twinx()
    fractions = (
        behavior["accepted_event_fraction_of_raw_products"].to_numpy() * 100
    )
    twin.plot(
        x,
        fractions,
        color=COLORS["orange"],
        marker="o",
        linewidth=1.3,
        markersize=5,
    )
    twin.set_ylabel("Accepted raw products (%)", color=COLORS["orange"])
    twin.set_ylim(0, 60)
    twin.tick_params(axis="y", colors=COLORS["ink"])
    for xi, value in zip(x, fractions):
        twin.text(
            xi,
            value + 1.9,
            f"{value:.1f}%",
            ha="center",
            va="bottom",
            fontsize=5.8,
            color=COLORS["orange"],
        )
    quiet_axes(ax, "y")
    panel_label(ax, "D", x=-0.10, y=1.05)


def draw_panel_e(ax, overlap: pd.DataFrame) -> None:
    tier_order = ["T1", "T2", "T3"]
    matrix = np.eye(3)
    lookup = {tier: index for index, tier in enumerate(tier_order)}
    for _, row in overlap.iterrows():
        i = lookup[str(row["left_tier"])]
        j = lookup[str(row["right_tier"])]
        matrix[i, j] = matrix[j, i] = float(row["jaccard_similarity"])
    image = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(3), tier_order)
    ax.set_yticks(range(3), tier_order)
    ax.set_title("Pairwise overlap of one-step structure sets", pad=8)
    for i in range(3):
        for j in range(3):
            ax.text(
                j,
                i,
                f"{matrix[i, j]:.3f}",
                ha="center",
                va="center",
                color="white" if matrix[i, j] > 0.55 else COLORS["ink"],
                fontsize=6.2,
            )
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard similarity")
    panel_label(ax, "E", x=-0.29, y=1.10)


def build_figure(data: dict[str, object], output: Path) -> list[Path]:
    fig = plt.figure(figsize=(15.2, 13.2))
    grid = fig.add_gridspec(
        3,
        12,
        height_ratios=[1.02, 1.72, 1.0],
        hspace=0.40,
        wspace=0.78,
    )

    ax_a = fig.add_subplot(grid[0, :4])
    draw_panel_a(ax_a, data["provenance"])

    ax_b = fig.add_subplot(grid[0, 4:])
    draw_panel_b(
        ax_b,
        data["pipeline"],
        data["tier_qc"],
        data["behavior"],
    )

    ax_c = fig.add_subplot(grid[1, :])
    draw_panel_c_v5(ax_c, data["grammar_examples"])

    ax_d = fig.add_subplot(grid[2, :6])
    draw_panel_d(ax_d, data["behavior"])

    ax_e = fig.add_subplot(grid[2, 6:])
    draw_panel_e(ax_e, data["overlap"])

    fig.suptitle(
        "Evidence-stratified construction and selection of the T1 reaction grammar",
        fontsize=14.2,
        fontweight="bold",
        y=0.988,
    )
    fig.text(
        0.5,
        0.018,
        "T1 combines the highest evidence tier with the broadest accepted "
        "one-step taxane coverage; T2 and T3 remain sensitivity and "
        "exploratory layers, respectively.",
        ha="center",
        va="bottom",
        fontsize=6.6,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.subplots_adjust(
        left=0.115,
        right=0.985,
        top=0.935,
        bottom=0.068,
    )

    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    stem = figure_dir / "Figure_1_Evidence_Stratified_Grammar_V6"
    outputs = []
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 400}),
    ):
        path = stem.with_suffix(suffix)
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def write_release(
    data: dict[str, object],
    audit: pd.DataFrame,
    output: Path,
    script_path: Path,
    figure_paths: list[Path],
) -> list[Path]:
    source_dir = output / "source_data"
    workflow_dir = output / "workflow"
    source_dir.mkdir(parents=True, exist_ok=True)
    workflow_dir.mkdir(parents=True, exist_ok=True)

    source_tables = {
        "Figure_1A_database_provenance.tsv": data["provenance"],
        "Figure_1B_grammar_construction.tsv": data["pipeline"],
        "Figure_1B_tier_quality_gates.tsv": data["tier_qc"],
        "Figure_1C_grammar_examples.tsv": data["grammar_examples"],
        "Figure_1D_tier_behavior.tsv": data["behavior"],
        "Figure_1E_pairwise_overlap.tsv": data["overlap"],
    }
    written: list[Path] = []
    for name, frame in source_tables.items():
        path = source_dir / name
        frame.to_csv(path, sep="\t", index=False)
        written.append(path)

    audit_path = output / "FIGURE_1_NUMERICAL_AUDIT.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    written.append(audit_path)

    caption_path = output / "FIGURE_1_CAPTION.md"
    caption_path.write_text(
        """# Figure 1 | Evidence-stratified construction and selection of the T1 reaction grammar

**(A)** Database-level provenance of 630,280 normalized reaction records. Colours distinguish curated domain anchors, biochemical resources, integrated reaction or rule resources, and reaction-chemistry provenance; colour does not encode confidence. **(B)** Audited construction from normalized reactions through deduplicated templates or anchors and 353,524 generalized reaction SMARTS. The upper workflow records the common normalization, principal-pair projection, deduplication, atom mapping, reaction-centre abstraction, and mutually exclusive T1/T2/T3 evidence assignment. The lower-left chart summarizes tier-specific release, quality-control, and initial-executability gates, whereas the lower-right boxes report the frozen compiled grammars and their intended analytical roles. T1 contains the highest-evidence biochemical layer and separately traceable curated taxane-domain productions. **(C)** Twelve complete molecular structures from six curated Taxol-pathway substrate-product pairs illustrate C-H hydroxylation, O-acetyl transfer, alcohol oxidation, O-deacetylation, O-benzoyl transfer, and N-benzoyl transfer. Substrate and product depictions share a common core orientation, and orange highlighting identifies the displayed reaction centre. The corresponding frozen T1 domain-consensus SMARTS shows how each whole-molecule transformation is abstracted into a directional, transferable local rewrite. These pathway pairs illustrate existing molecular edits and are not newly inferred reactions. **(D)** Unique one-step structures and accepted fractions of raw generated products for the three mutually exclusive evidence layers. **(E)** Pairwise Jaccard overlap of the corresponding G1 structure sets. T1 was selected for primary inference because it combines the highest evidence tier with the broadest accepted one-step coverage; T2 and T3 were retained as sensitivity and exploratory layers.
""",
        encoding="utf-8",
    )
    written.append(caption_path)

    readme_path = output / "README.md"
    readme_path.write_text(
        """# Figure 1 Redesign V6

This release recomposes Figure 1 from frozen source tables only. It does not
rebuild reaction rules, reassign evidence tiers, rerun molecular generation,
or change any headline numerical result.

The visual hierarchy is:

1. database provenance;
2. the original V3 evidence-normalization, tier-quality-gate, and grammar-
   selection panel;
3. twelve complete molecular structures from six curated Taxol-pathway
   reactions, paired with the corresponding generalized T1 productions;
4. one-step coverage and product acceptance;
5. pairwise overlap of evidence-layer outputs.

Panel B is retained directly from Figure 1 V3. Panel C uses the frozen pathway
molecules for illustration only and computes reaction-centre highlighting
solely for figure rendering. Panels D and E retain the original dual-axis and
Jaccard-matrix designs.
Previous Figure 1 releases are preserved unchanged.
""",
        encoding="utf-8",
    )
    written.append(readme_path)

    workflow_copy = workflow_dir / script_path.name
    shutil.copy2(script_path, workflow_copy)
    written.append(workflow_copy)

    manifest_rows = []
    authoritative = data["authoritative_sources"]
    panel_sources = {
        "A": [authoritative["provenance"]],
        "B": [
            authoritative["stage"],
            authoritative["tier_qc"],
            authoritative["behavior"],
        ],
        "C": [authoritative["grammar"]],
        "D": [authoritative["behavior"]],
        "E": [authoritative["overlap"]],
    }
    panel_sources["C"].append(authoritative["taxol_pathway"])
    panel_files = {
        "A": source_dir / "Figure_1A_database_provenance.tsv",
        "B": source_dir / "Figure_1B_grammar_construction.tsv",
        "C": source_dir / "Figure_1C_grammar_examples.tsv",
        "D": source_dir / "Figure_1D_tier_behavior.tsv",
        "E": source_dir / "Figure_1E_pairwise_overlap.tsv",
    }
    for panel in ("A", "B", "C", "D", "E"):
        manifest_rows.append(
            {
                "panel": panel,
                "source_data_file": str(
                    panel_files[panel].relative_to(output)
                ),
                "source_data_sha256": sha256(panel_files[panel]),
                "authoritative_inputs": ";".join(
                    str(path) for path in panel_sources[panel]
                ),
                "authoritative_input_sha256": ";".join(
                    sha256(path) for path in panel_sources[panel]
                ),
                "operation": "display-only extraction or annotation",
                "rendering_script": f"workflow/{script_path.name}",
                "rendering_script_sha256": sha256(script_path),
            }
        )
    manifest_path = output / "FIGURE_1_SOURCE_MANIFEST.tsv"
    pd.DataFrame(manifest_rows).to_csv(
        manifest_path, sep="\t", index=False
    )
    written.append(manifest_path)

    summary_path = output / "BUILD_SUMMARY.json"
    payload = {
        "release": "Figure 1 redesign V6",
        "scientific_recalculation_performed": False,
        "figures": [
            {
                "path": str(path.relative_to(output)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in figure_paths
        ],
        "audit_status": "PASS",
        "source_manifest": str(manifest_path.relative_to(output)),
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    written.append(summary_path)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    data = load_frozen_data(args.work)
    audit = validate(data)
    figures = build_figure(data, args.output)
    written = write_release(
        data,
        audit,
        args.output,
        Path(__file__).resolve(),
        figures,
    )
    print(f"Figure 1 redesign built at {args.output}")
    for path in figures + written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
