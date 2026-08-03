#!/usr/bin/env python3
"""Recompose Figure 3A as frozen G0-to-G3 derivation chains.

This presentation-only workflow reads released nodes and derivation events.
It does not rebuild rules, enumerate products, or recalculate scientific
statistics. Panels B-D are rendered from the same frozen tables used in V6.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle

import build_figure2_figure3_v6 as base


DEFAULT_WORK = Path(".")
DEFAULT_OUTPUT = Path("figure3_output")

SPACE_DB = (
    "inputs/G0_G3_primary_release/03_primary_G0_G3/"
    "taxane_reaction_grammar_space.sqlite"
)

CHAIN_DEFINITIONS = [
    {
        "chain_id": "chain_1",
        "space_ids": [
            "G0_00029",
            "G1_00000695",
            "G2_00009939",
            "G3_00105258",
        ],
        "accent": base.COLORS["blue"],
        "description": (
            "Anchor-derived hydroxylation initiates iterative "
            "oxygenation and acyl-state diversification"
        ),
    },
    {
        "chain_id": "chain_2",
        "space_ids": [
            "G0_00129",
            "G1_00003186",
            "G2_00045680",
            "G3_00485701",
        ],
        "accent": base.COLORS["orange"],
        "description": (
            "Anchor-derived O-acetyl transfer initiates iterative "
            "oxygenation and acyl-state diversification"
        ),
    },
]

GENERATION_COLORS = {
    0: base.COLORS["red"],
    1: base.COLORS["blue"],
    2: base.COLORS["teal"],
    3: base.COLORS["orange"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def query_frame(
    connection: sqlite3.Connection,
    query: str,
    parameters: tuple[object, ...],
) -> pd.DataFrame:
    return pd.read_sql_query(query, connection, params=parameters)


def load_derivation_chains(
    work: Path,
    migration: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, Path]:
    database = work / SPACE_DB
    if not database.exists():
        raise FileNotFoundError(database)

    local_edit_by_rule = (
        migration.set_index("grammar_rule_id")["anchor_local_edit"]
        .astype(str)
        .str.replace("  ", " ", regex=False)
        .to_dict()
    )
    anchor_by_rule = (
        migration.set_index("grammar_rule_id")[
            [
                "anchor_enzyme",
                "anchor_ec",
                "anchor_reaction_smarts",
            ]
        ]
        .to_dict("index")
    )

    structure_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    with sqlite3.connect(database) as connection:
        for chain_order, definition in enumerate(CHAIN_DEFINITIONS, start=1):
            space_ids = definition["space_ids"]
            placeholders = ",".join("?" for _ in space_ids)
            nodes = query_frame(
                connection,
                f"""
                SELECT
                    space_id,
                    generation_first,
                    molecule_names,
                    smiles,
                    formula
                FROM nodes
                WHERE space_id IN ({placeholders})
                """,
                tuple(space_ids),
            ).set_index("space_id")
            missing = [value for value in space_ids if value not in nodes.index]
            if missing:
                raise RuntimeError(
                    f"Missing frozen nodes for {definition['chain_id']}: {missing}"
                )

            for generation, space_id in enumerate(space_ids):
                node = nodes.loc[space_id]
                if int(node["generation_first"]) != generation:
                    raise RuntimeError(
                        f"Generation mismatch for {space_id}: "
                        f"{node['generation_first']} != {generation}"
                    )
                structure_rows.append(
                    {
                        "chain_id": definition["chain_id"],
                        "chain_order": chain_order,
                        "chain_description": definition["description"],
                        "generation": generation,
                        "space_id": space_id,
                        "molecule_name": str(node["molecule_names"] or ""),
                        "smiles": str(node["smiles"]),
                        "formula": str(node["formula"]),
                    }
                )

            for generation in range(1, 4):
                source_id = space_ids[generation - 1]
                target_id = space_ids[generation]
                event = query_frame(
                    connection,
                    """
                    SELECT
                        event_id,
                        generation,
                        source_space_id,
                        target_space_id,
                        target_is_new,
                        grammar_rule_id,
                        smarts_rule_id,
                        semantic_group_id,
                        reaction_type,
                        evidence_layer,
                        final_rule_confidence,
                        observed_element_delta,
                        source_atom_retention,
                        observed_changed_source_atoms,
                        source_product_tanimoto,
                        immediate_reverse_cycle
                    FROM derivation_events
                    WHERE source_space_id = ?
                      AND target_space_id = ?
                      AND generation = ?
                    """,
                    (source_id, target_id, generation),
                )
                if len(event) != 1:
                    raise RuntimeError(
                        "Expected one frozen derivation event for "
                        f"{source_id} -> {target_id}; observed {len(event)}"
                    )
                row = event.iloc[0].to_dict()
                rule_id = str(row["grammar_rule_id"])
                anchor = anchor_by_rule.get(rule_id, {})
                row.update(
                    {
                        "chain_id": definition["chain_id"],
                        "chain_order": chain_order,
                        "step_order": generation,
                        "source_generation": generation - 1,
                        "target_generation": generation,
                        "local_edit": local_edit_by_rule.get(
                            rule_id,
                            str(row["reaction_type"]).replace("_", " "),
                        ),
                        "anchor_enzyme": anchor.get("anchor_enzyme", ""),
                        "anchor_ec": anchor.get("anchor_ec", ""),
                        "reaction_smarts": anchor.get(
                            "anchor_reaction_smarts",
                            "",
                        ),
                    }
                )
                edge_rows.append(row)

    structures = pd.DataFrame(structure_rows).sort_values(
        ["chain_order", "generation"]
    )
    edges = pd.DataFrame(edge_rows).sort_values(
        ["chain_order", "step_order"]
    )
    return structures, edges, database


def draw_arrow(
    ax,
    x0: float,
    x1: float,
    y: float,
    color: str,
) -> None:
    arrow = FancyArrowPatch(
        (x0, y),
        (x1, y),
        transform=ax.transAxes,
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=1.15,
        color=color,
        shrinkA=0,
        shrinkB=0,
        clip_on=False,
    )
    ax.add_patch(arrow)


def draw_chain_row(
    ax,
    structures: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    y: float,
    accent: str,
) -> None:
    row_height = 0.39
    ax.add_patch(
        Rectangle(
            (0.012, y),
            0.976,
            row_height,
            transform=ax.transAxes,
            facecolor=base.COLORS["white"],
            edgecolor="#d6dde2",
            linewidth=0.7,
        )
    )
    ax.add_patch(
        Rectangle(
            (0.012, y),
            0.009,
            row_height,
            transform=ax.transAxes,
            facecolor=accent,
            edgecolor="none",
        )
    )

    description = str(structures.iloc[0]["chain_description"])
    ax.text(
        0.032,
        y + row_height * 0.925,
        description,
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=5.5,
        fontweight="bold",
        color=base.COLORS["ink"],
    )
    ax.text(
        0.972,
        y + row_height * 0.925,
        "three sequential accepted T1 derivations",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=4.6,
        color=base.COLORS["gray"],
    )

    molecule_x = [0.050, 0.292, 0.534, 0.776]
    molecule_w = 0.165
    molecule_y = y + 0.055
    molecule_h = row_height * 0.55
    edge_by_step = {
        int(row.step_order): row
        for row in edges.itertuples(index=False)
    }

    previous_smiles: str | None = None
    for row in structures.itertuples(index=False):
        generation = int(row.generation)
        x = molecule_x[generation]
        highlighted_atoms: set[int] = set()
        highlighted_bonds: set[tuple[int, int]] = set()
        if previous_smiles is not None:
            _, highlighted_atoms, _, highlighted_bonds = (
                base.molecule_difference(previous_smiles, row.smiles)
            )
        molecule_ax = ax.inset_axes(
            [x, molecule_y, molecule_w, molecule_h]
        )
        base.draw_molecule(
            molecule_ax,
            row.smiles,
            highlighted_atoms=highlighted_atoms,
            highlighted_bonds=highlighted_bonds,
        )

        generation_label = (
            "G0 | known taxane"
            if generation == 0
            else f"G{generation} | generated product"
        )
        ax.text(
            x + molecule_w / 2,
            y + row_height * 0.815,
            generation_label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=5.0,
            fontweight="bold",
            color=GENERATION_COLORS[generation],
            bbox={
                "boxstyle": "round,pad=0.22",
                "facecolor": base.COLORS["white"],
                "edgecolor": GENERATION_COLORS[generation],
                "linewidth": 0.6,
            },
        )
        if generation == 0:
            ax.text(
                x + molecule_w / 2,
                y + row_height * 0.735,
                str(row.molecule_name),
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=4.1,
                color=base.COLORS["ink"],
            )
        ax.text(
            x + molecule_w / 2,
            y + row_height * 0.055,
            base.formula_label(row.formula),
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=5.2,
            color=base.COLORS["ink"],
        )
        previous_smiles = str(row.smiles)

    for step in range(1, 4):
        edge = edge_by_step[step]
        left = molecule_x[step - 1] + molecule_w + 0.004
        right = molecule_x[step] - 0.004
        middle = (left + right) / 2
        color = GENERATION_COLORS[step]
        draw_arrow(
            ax,
            left,
            right,
            y + row_height * 0.43,
            color,
        )
        local_edit = str(edge.local_edit).replace(" -> ", " → ")
        ax.text(
            middle,
            y + row_height * 0.60,
            local_edit,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=4.6,
            fontweight="bold",
            color=color,
        )
        ax.text(
            middle,
            y + row_height * 0.27,
            f"T1 grammar\nTanimoto {float(edge.source_product_tanimoto):.3f}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=4.1,
            color=base.COLORS["gray"],
            linespacing=1.15,
        )


def draw_panel_a(
    ax,
    structures: pd.DataFrame,
    edges: pd.DataFrame,
) -> None:
    ax.set_axis_off()
    ax.set_title(
        "Representative iterative T1 derivation chains across G0-G3",
        fontsize=10.5,
        pad=6,
    )
    base.panel_label(ax, "A", x=-0.025, y=1.00)
    for y, definition in zip([0.53, 0.08], CHAIN_DEFINITIONS):
        chain_id = definition["chain_id"]
        draw_chain_row(
            ax,
            structures[structures["chain_id"].eq(chain_id)],
            edges[edges["chain_id"].eq(chain_id)],
            y=y,
            accent=str(definition["accent"]),
        )
    ax.text(
        0.5,
        0.015,
        (
            "Orange molecular highlights identify the local structural "
            "difference from the immediately preceding generation."
        ),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=4.7,
        color=base.COLORS["gray"],
    )


def draw_panel_b(ax, locality: pd.DataFrame) -> None:
    pivot = locality.pivot(
        index="generation",
        columns="edit_locality",
        values="derivation_event_count",
    ).fillna(0)
    pivot = pivot.reindex(columns=["0", "1", "2-3", ">3"], fill_value=0)
    fractions = pivot.div(pivot.sum(axis=1), axis=0)
    bottom = np.zeros(len(fractions))
    categories = [
        ("0", base.COLORS["gray"]),
        ("1", base.COLORS["blue"]),
        ("2-3", base.COLORS["teal"]),
        (">3", base.COLORS["orange"]),
    ]
    for category, color in categories:
        ax.bar(
            [f"G{int(value)}" for value in fractions.index],
            fractions[category],
            bottom=bottom,
            label=category,
            color=color,
            edgecolor=base.COLORS["ink"],
            linewidth=0.25,
        )
        bottom += fractions[category].to_numpy()
    totals = pivot.sum(axis=1).astype(int)
    for index, generation in enumerate(fractions.index):
        ax.text(
            index,
            0.885,
            f"{100 * fractions.loc[generation, '1']:.1f}%\n1-atom",
            ha="center",
            va="center",
            fontsize=5.6,
            color=base.COLORS["white"],
            fontweight="bold",
        )
        ax.text(
            index,
            -0.075,
            f"n={int(totals.loc[generation]):,}",
            ha="center",
            va="top",
            fontsize=5.1,
            color=base.COLORS["gray"],
        )
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of accepted events")
    ax.set_title("Reaction-edit locality")
    ax.legend(
        title="changed source atoms",
        frameon=False,
        ncol=2,
        fontsize=5.1,
        title_fontsize=5.2,
        loc="lower left",
    )
    base.quiet_axes(ax, "y")
    base.panel_label(ax, "B")


def draw_panel_c(ax, functional: pd.DataFrame) -> None:
    functional_summary = functional.sort_values(
        "derivation_event_count",
        ascending=True,
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
        "no_counted_functional_state_change": (
            "no counted functional-state change"
        ),
        "free_hydroxyl:-2;carboxylic_acid_or_carboxylate:+1": (
            "free hydroxyl -2; carboxylic acid/carboxylate +1"
        ),
    }
    labels = [
        functional_labels.get(value, str(value).replace("_", " "))
        for value in functional_summary["functional_state_transition"]
    ]
    ax.barh(
        labels,
        functional_summary["derivation_event_count"],
        color=base.COLORS["teal"],
    )
    ax.set_xscale("log")
    ax.set_xlabel("Accepted derivation events (log scale)")
    ax.set_title("Dominant functional-state transitions")
    ax.tick_params(axis="y", labelsize=5.0)
    base.quiet_axes(ax, "x")
    base.panel_label(ax, "C")


def draw_panel_d(ax, elemental: pd.DataFrame) -> None:
    elemental_summary = elemental.sort_values("derivation_event_count")
    ax.barh(
        [
            base.elemental_delta_label(value)
            for value in elemental_summary["observed_element_delta"]
        ],
        elemental_summary["derivation_event_count"],
        color=base.COLORS["orange"],
    )
    ax.set_xscale("log")
    ax.set_xlabel("Accepted derivation events (log scale)")
    ax.set_title("Elemental edit signatures")
    ax.tick_params(axis="y", labelsize=5.1)
    base.quiet_axes(ax, "x")
    base.panel_label(ax, "D")


def build_figure(
    data: dict[str, object],
    structures: pd.DataFrame,
    edges: pd.DataFrame,
    output: Path,
) -> list[Path]:
    locality = data["locality"]
    functional = data["functional"]
    elemental = data["elemental"]

    fig = plt.figure(figsize=(14.2, 9.3))
    grid = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.46, 1.0],
        hspace=0.24,
        wspace=0.52,
    )
    draw_panel_a(fig.add_subplot(grid[0, :]), structures, edges)
    draw_panel_b(fig.add_subplot(grid[1, 0]), locality)
    draw_panel_c(fig.add_subplot(grid[1, 1]), functional)
    draw_panel_d(fig.add_subplot(grid[1, 2]), elemental)

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
        color=base.COLORS["gray"],
        bbox={
            "boxstyle": "round,pad=0.34",
            "facecolor": base.COLORS["white"],
            "edgecolor": "#d8dde2",
            "linewidth": 0.55,
        },
    )
    fig.subplots_adjust(
        left=0.075,
        right=0.985,
        top=0.900,
        bottom=0.082,
    )
    return base.save_figure(
        fig,
        output,
        "Figure_3_T1_Chemical_Regularities_V7",
    )


def validate(
    structures: pd.DataFrame,
    edges: pd.DataFrame,
    locality: pd.DataFrame,
) -> pd.DataFrame:
    expected_events = {
        "chain_1": [814, 40152, 837078],
        "chain_2": [3683, 117446, 2020022],
    }
    rows: list[dict[str, object]] = []
    for chain_id, event_ids in expected_events.items():
        chain_structures = structures[
            structures["chain_id"].eq(chain_id)
        ].sort_values("generation")
        chain_edges = edges[edges["chain_id"].eq(chain_id)].sort_values(
            "step_order"
        )
        checks = {
            "four_generations": chain_structures["generation"].tolist()
            == [0, 1, 2, 3],
            "three_derivation_events": chain_edges["event_id"].astype(int).tolist()
            == event_ids,
            "directional_generation_order": (
                chain_edges["source_generation"].astype(int).tolist()
                == [0, 1, 2]
                and chain_edges["target_generation"].astype(int).tolist()
                == [1, 2, 3]
            ),
            "new_products": chain_edges["target_is_new"].astype(int).eq(1).all(),
            "no_immediate_reverse_cycles": (
                chain_edges["immediate_reverse_cycle"].astype(int).eq(0).all()
            ),
            "single_source_atom_edits": (
                chain_edges["observed_changed_source_atoms"]
                .astype(int)
                .eq(1)
                .all()
            ),
        }
        for check, status in checks.items():
            rows.append(
                {
                    "claim_id": f"{chain_id}_{check}",
                    "observed_value": str(bool(status)),
                    "expected_value": "True",
                    "status": "PASS" if status else "FAIL",
                }
            )

    totals = (
        locality.groupby("generation")["derivation_event_count"]
        .sum()
        .astype(int)
        .to_dict()
    )
    expected_totals = {1: 18504, 2: 490855, 3: 7418850}
    for generation, expected in expected_totals.items():
        observed = int(totals[generation])
        rows.append(
            {
                "claim_id": f"G{generation}_accepted_derivation_events",
                "observed_value": observed,
                "expected_value": expected,
                "status": "PASS" if observed == expected else "FAIL",
            }
        )
    audit = pd.DataFrame(rows)
    if not audit["status"].eq("PASS").all():
        raise RuntimeError(
            "Figure 3 V7 audit failed:\n"
            + audit[audit["status"].ne("PASS")].to_string(index=False)
        )
    return audit


def write_release(
    output: Path,
    script: Path,
    database: Path,
    data: dict[str, object],
    structures: pd.DataFrame,
    edges: pd.DataFrame,
    audit: pd.DataFrame,
    figure_paths: list[Path],
) -> None:
    source_dir = output / "source_data"
    workflow_dir = output / "workflow"
    source_dir.mkdir(parents=True, exist_ok=True)
    workflow_dir.mkdir(parents=True, exist_ok=True)

    structure_path = source_dir / "Figure_3_A_G0_G3_chain_structures.tsv"
    edge_path = source_dir / "Figure_3_A_G0_G3_chain_edges.tsv"
    locality_path = source_dir / "Figure_3_B_edit_locality.tsv"
    functional_path = source_dir / "Figure_3_C_functional_state_transitions.tsv"
    elemental_path = source_dir / "Figure_3_D_elemental_edit_signatures.tsv"
    structures.to_csv(structure_path, sep="\t", index=False)
    edges.to_csv(edge_path, sep="\t", index=False)
    data["locality"].to_csv(locality_path, sep="\t", index=False)
    data["functional"].to_csv(functional_path, sep="\t", index=False)
    data["elemental"].to_csv(elemental_path, sep="\t", index=False)
    audit.to_csv(output / "NUMERICAL_AUDIT_V7.tsv", sep="\t", index=False)
    shutil.copy2(script, workflow_dir / script.name)
    shutil.copy2(
        Path(base.__file__),
        workflow_dir / Path(base.__file__).name,
    )

    caption = """# Figure 3 | Transferable local chemistry underlying T1 space expansion

**(A)** Two representative iterative derivation chains traced through the
frozen reaction-grammar space from known taxanes (G0) to accepted one-rule
(G1), two-rule (G2), and three-rule (G3) products. Every arrow is a recorded
directional derivation event and is labelled by its transferred local edit and
source-product Tanimoto similarity. Orange molecular highlights identify the
structural difference from the immediately preceding generation and are
depiction aids rather than atom-mapped reaction centers. Generated structures
are reaction-grammar-accessible hypotheses, not experimentally confirmed
metabolites. **(B)** Distribution of reaction-edit locality across generations;
percentages identify events with one changed source atom and labels report
generation-specific accepted-event denominators. **(C)** Dominant
structure-derived functional-state transitions. **(D)** Most frequent
elemental edit signatures. The header reports all 7,928,209 accepted derivation
events used for the frozen generation-level summaries.
"""
    (output / "FIGURE_3_CAPTION_V7.md").write_text(
        caption,
        encoding="utf-8",
    )

    manifest = pd.DataFrame(
        [
            {
                "panel": "A structures",
                "source_data_file": str(structure_path.relative_to(output)),
                "source_data_sha256": sha256(structure_path),
                "authoritative_input": str(database),
                "operation": "frozen node extraction",
            },
            {
                "panel": "A derivation edges",
                "source_data_file": str(edge_path.relative_to(output)),
                "source_data_sha256": sha256(edge_path),
                "authoritative_input": str(database),
                "operation": "frozen directional-edge extraction",
            },
            {
                "panel": "B",
                "source_data_file": str(locality_path.relative_to(output)),
                "source_data_sha256": sha256(locality_path),
                "authoritative_input": str(data["paths"]["locality"]),
                "operation": "unchanged frozen summary",
            },
            {
                "panel": "C",
                "source_data_file": str(functional_path.relative_to(output)),
                "source_data_sha256": sha256(functional_path),
                "authoritative_input": str(data["paths"]["functional"]),
                "operation": "unchanged frozen summary",
            },
            {
                "panel": "D",
                "source_data_file": str(elemental_path.relative_to(output)),
                "source_data_sha256": sha256(elemental_path),
                "authoritative_input": str(data["paths"]["elemental"]),
                "operation": "unchanged frozen summary",
            },
        ]
    )
    manifest.to_csv(
        output / "FIGURE_3_SOURCE_MANIFEST_V7.tsv",
        sep="\t",
        index=False,
    )

    summary = {
        "release": "Figure 3 V7",
        "scientific_recalculation_performed": False,
        "panel_A_generation_scope": "G0-G3",
        "panel_A_chain_count": 2,
        "panel_A_structure_count": int(len(structures)),
        "panel_A_derivation_event_count": int(len(edges)),
        "panels_B_to_D": "unchanged frozen V6 source tables",
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
    (output / "BUILD_SUMMARY_V7.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    (output / "README.md").write_text(
        """# Figure 3 V7

This presentation-only release extends panel A from G0-to-G1 examples to two
complete, recorded G0-to-G3 derivation chains. The structures and directional
events are extracted from the frozen SQLite release. Panels B-D use unchanged
V6 source tables. No rule construction, product enumeration, descriptor
calculation, fingerprint calculation, or scientific benchmark was rerun.
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base.configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    data = base.load_inputs(args.work)
    structures, edges, database = load_derivation_chains(
        args.work,
        data["migration"],
    )
    audit = validate(structures, edges, data["locality"])
    figure_paths = build_figure(
        data,
        structures,
        edges,
        args.output,
    )
    write_release(
        args.output,
        Path(__file__).resolve(),
        database,
        data,
        structures,
        edges,
        audit,
        figure_paths,
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figures": [str(path) for path in figure_paths],
                "audit_status": "PASS",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
