from __future__ import annotations

import json
import os
from pathlib import Path
import re
import textwrap
from typing import Any

import numpy as np
import pandas as pd

from .io import ensure_dir, read_table, write_json, write_table


os.environ.setdefault("MPLCONFIGDIR", "/tmp/taxane_reaction_grammar_matplotlib")


PALETTE = {
    "G0": "#3B3B3B",
    "G1": "#0072B2",
    "G2": "#009E73",
    "G3": "#E69F00",
    "T1": "#0072B2",
    "T2": "#009E73",
    "T3": "#D55E00",
    "gray": "#8C8C8C",
    "light": "#D8D8D8",
    "red": "#C44E52",
    "purple": "#7A5195",
}

DESCRIPTOR_LABELS = {
    "exact_mass": "exact mass",
    "clogp": "cLogP",
    "tpsa": "tPSA",
    "hbd": "H-bond donors",
    "hba": "H-bond acceptors",
    "rotatable_bonds": "rotatable bonds",
    "fraction_csp3": "fraction Csp3",
    "ring_count": "ring count",
}

REACTION_TYPE_LABELS = {
    "large_side_chain_transfer_or_loss": "large side-chain transfer/loss",
    "acetylation_or_deacetylation_like_acyl_transfer": "acetyl transfer-like",
    "benzoylation_or_aromatic_acyl_transfer": "benzoyl/aromatic acyl transfer",
    "N_benzoylation_or_aromatic_acyl_transfer": "N-benzoyl/aromatic acyl transfer",
    "cyclization_or_ring_rearrangement": "cyclization/ring rearrangement",
    "deacetylation_or_acetyl_ester_hydrolysis": (
        "deacetylation/ester hydrolysis"
    ),
    "oxidation_or_dehydrogenation": "oxidation/dehydrogenation",
    "hydroxylation_or_oxygenation": "hydroxylation/oxygenation",
}

FUNCTIONAL_STATE_LABELS = {
    "free_hydroxyl": "OH",
    "carboxylic_acid_or_carboxylate": "carboxyl",
    "ketone_or_aldehyde": "carbonyl",
    "no_counted_functional_state_change": "no counted change",
}


def _plt():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.titlelocation": "left",
            "axes.titlepad": 8.0,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.75,
            "ytick.major.width": 0.75,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.transparent": False,
        }
    )
    return plt


def _panel(axis, label: str) -> None:
    axis.text(
        -0.17,
        1.13,
        label,
        transform=axis.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )


def _clean(axis, *, grid: bool = True) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if grid:
        axis.grid(axis="y", color="#E6E6E6", linewidth=0.5, zorder=0)
        axis.set_axisbelow(True)


def _wrap(value: object, width: int = 32) -> str:
    return "\n".join(
        textwrap.wrap(
            str(value).replace("_", " "),
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )


def _descriptor_label(value: object) -> str:
    return DESCRIPTOR_LABELS.get(str(value), str(value).replace("_", " "))


def _functional_state_transition_label(value: object) -> str:
    parts = []
    for token in str(value).split(";"):
        if ":" in token:
            state, delta = token.rsplit(":", 1)
            label = FUNCTIONAL_STATE_LABELS.get(
                state, state.replace("_", " ")
            )
            parts.append(f"{label}({delta})")
        else:
            parts.append(
                FUNCTIONAL_STATE_LABELS.get(
                    token, token.replace("_", " ")
                )
            )
    return _wrap(", ".join(parts), 34)


def _reaction_type_label(value: object) -> str:
    return REACTION_TYPE_LABELS.get(
        str(value), str(value).replace("_", " ")
    )


def _rule_label(row: pd.Series, width: int = 27) -> str:
    rule_id = str(row.get("grammar_rule_id", "rule"))
    suffix_match = re.search(r"(\d+)$", rule_id)
    suffix: object = (
        int(suffix_match.group(1)) if suffix_match else rule_id
    )
    prefix = "TD" if "TAXANE_DOMAIN" in rule_id else "T1"
    reaction_type = _reaction_type_label(row.get("reaction_type", ""))
    return f"{prefix}-{suffix}: {reaction_type}" if reaction_type else f"{prefix}-{suffix}"


def _format_element_delta(value: object) -> str:
    try:
        parsed = json.loads(str(value))
    except (TypeError, ValueError, json.JSONDecodeError):
        return _wrap(value, 30)
    if not parsed:
        return "No formula change"
    preferred = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
    elements = [element for element in preferred if element in parsed]
    elements.extend(sorted(set(parsed) - set(elements)))
    return ", ".join(
        f"Δ{element} {int(parsed[element]):+d}" for element in elements
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _portable_source_reference(source_path: Path | str) -> str:
    references = []
    for raw_reference in str(source_path).split(";"):
        path = Path(raw_reference)
        if not path.is_absolute():
            references.append(str(path))
            continue
        parts = path.parts
        if "taxane_space_study_outputs" in parts:
            root_index = parts.index("taxane_space_study_outputs")
            references.append(str(Path(*parts[root_index + 1 :])))
        else:
            references.append(path.name)
    return ";".join(references)


class SourceRecorder:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.source_dir = ensure_dir(output_dir / "source_data")
        self.rows: list[dict[str, str]] = []

    def add(
        self,
        figure: str,
        panel: str,
        frame: pd.DataFrame,
        *,
        source_path: Path | str,
        description: str,
    ) -> Path:
        scaffold_columns = [
            column for column in frame.columns if "scaffold" in column.lower()
        ]
        if scaffold_columns:
            raise ValueError(
                f"{figure}{panel} contains prohibited scaffold fields: "
                f"{scaffold_columns}"
            )
        path = self.source_dir / f"{figure}_{panel}_source_data.tsv"
        write_table(frame, path)
        self.rows.append(
            {
                "figure": figure,
                "panel": panel,
                "source_data_file": str(path.relative_to(self.output_dir)),
                "authoritative_input": _portable_source_reference(source_path),
                "description": description,
            }
        )
        return path

    def write(self, output_dir: Path) -> Path:
        path = output_dir / "figure_source_data_manifest.tsv"
        write_table(pd.DataFrame(self.rows), path)
        return path


def _save(figure, output_dir: Path, stem: str) -> list[Path]:
    paths = []
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 400}),
    ):
        path = output_dir / f"{stem}{suffix}"
        figure.savefig(path, **kwargs)
        paths.append(path)
    return paths


def _fig1(
    provenance_dir: Path,
    analysis_dir: Path,
    output_dir: Path,
    recorder: SourceRecorder,
) -> list[Path]:
    plt = _plt()
    stages_path = provenance_dir / "rule_library_build_stage_counts.tsv"
    sources_path = provenance_dir / "source_database_rule_contributions.tsv"
    tiers_path = provenance_dir / "prepared_grammar_tier_composition.tsv"
    final_sources_path = provenance_dir / "final_primary_grammar_source_support.tsv"
    reactions_path = provenance_dir / "final_primary_grammar_reaction_types.tsv"
    rule_usage_path = analysis_dir / "reaction_grammar_usage.tsv"
    stages = read_table(stages_path)
    sources = read_table(sources_path)
    tiers = read_table(tiers_path)
    final_sources = read_table(final_sources_path)
    reactions = read_table(reactions_path)
    rule_usage = read_table(rule_usage_path)
    stages["record_count"] = pd.to_numeric(stages["record_count"], errors="coerce")
    numeric_source = [
        column
        for column in sources
        if column.endswith("_records") or "supported" in column
    ]
    for column in numeric_source:
        sources[column] = pd.to_numeric(
            sources[column], errors="coerce"
        ).fillna(0)
    for column in tiers.columns:
        if column.endswith("rows") or column in {
            "semantic_groups",
            "initial_executable_representatives",
        }:
            tiers[column] = pd.to_numeric(tiers[column], errors="coerce")

    figure, axes = plt.subplots(2, 3, figsize=(10.2, 7.4))
    axes = axes.ravel()

    selected_names = [
        "normalized source reactions",
        "main-substrate/main-product projections",
        "deduplicated templates or anchors",
        "predictive generalized SMARTS rules",
    ]
    stage_display = {
        "normalized source reactions": "normalized\nsource reactions",
        "main-substrate/main-product projections": "main-pair\nprojections",
        "deduplicated templates or anchors": "deduplicated\ntemplates/anchors",
        "predictive generalized SMARTS rules": "generalized\nSMARTS",
    }
    panel_a = stages[stages["build_stage"].isin(selected_names)].copy()
    panel_a["build_stage"] = pd.Categorical(
        panel_a["build_stage"], selected_names, ordered=True
    )
    panel_a = panel_a.sort_values("build_stage")
    y = np.arange(len(panel_a))
    axes[0].barh(
        y, panel_a["record_count"], color=PALETTE["gray"], height=0.62
    )
    axes[0].set_xscale("log")
    axes[0].set_yticks(
        y,
        [stage_display[str(value)] for value in panel_a["build_stage"]],
    )
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Records (log scale)")
    axes[0].set_title("Database records to generalized reaction grammar")
    for position, value in zip(y, panel_a["record_count"]):
        axes[0].text(
            value * 1.025,
            position,
            f"{int(value):,}",
            ha="left",
            va="center",
            fontsize=6,
        )
    axes[0].set_xlim(
        left=max(1, panel_a["record_count"].min() / 1.35),
        right=panel_a["record_count"].max() * 1.35,
    )
    _clean(axes[0], grid=False)
    _panel(axes[0], "A")
    recorder.add(
        "Fig1",
        "A",
        panel_a,
        source_path=stages_path,
        description="Selected source-to-generalized-rule construction stages.",
    )

    panel_b = sources.sort_values("normalized_source_reaction_records")
    axes[1].barh(
        [str(value).replace("_", " ") for value in panel_b["source_database"]],
        panel_b["normalized_source_reaction_records"],
        color=PALETTE["T1"],
    )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Normalized reaction records (log scale)")
    axes[1].set_title("Database contributions before rule abstraction")
    _clean(axes[1], grid=False)
    _panel(axes[1], "B")
    recorder.add(
        "Fig1",
        "B",
        panel_b[
            [
                "source_database",
                "normalized_source_reaction_records",
                "main_pair_records",
            ]
        ],
        source_path=sources_path,
        description="Normalized and main-pair reaction records by source database.",
    )

    panel_c = sources.copy().sort_values(
        "T1_rule_rows_supported", ascending=True
    )
    y = np.arange(len(panel_c))
    for offset, (tier, color, marker) in enumerate(
        (
            ("T1", PALETTE["T1"], "o"),
            ("T2", PALETTE["T2"], "s"),
            ("T3", PALETTE["T3"], "^"),
        )
    ):
        values = panel_c[f"{tier}_rule_rows_supported"].to_numpy(dtype=float)
        axes[2].scatter(
            np.where(values > 0, values, np.nan),
            y + (offset - 1) * 0.18,
            color=color,
            label=tier,
            marker=marker,
            s=18,
            zorder=3,
        )
    axes[2].set_yticks(
        y,
        [str(value).replace("_", " ") for value in panel_c["source_database"]],
    )
    axes[2].set_xscale("log")
    axes[2].set_xlim(left=0.8)
    axes[2].set_xlabel("Independently supported rule rows (log scale)")
    axes[2].set_title("Overlapping source support across evidence tiers")
    axes[2].legend(frameon=False, ncol=3)
    axes[2].grid(axis="x", color="#E6E6E6", linewidth=0.5, zorder=0)
    axes[2].set_axisbelow(True)
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)
    _panel(axes[2], "C")
    recorder.add(
        "Fig1",
        "C",
        panel_c[
            [
                "source_database",
                "T1_rule_rows_supported",
                "T2_rule_rows_supported",
                "T3_rule_rows_supported",
                "release_support_counting_scope",
            ]
        ],
        source_path=sources_path,
        description="Overlapping database support for exclusive T1/T2/T3 rule rows.",
    )

    panel_d = tiers.sort_values("tier")
    metrics = [
        ("input_release_rows", "release"),
        ("grammar_qc_eligible_rows", "QC-eligible"),
        ("semantic_groups", "semantic groups"),
        ("initial_executable_representatives", "representatives"),
    ]
    x = np.arange(len(panel_d))
    width = 0.19
    colors = [
        PALETTE["gray"],
        PALETTE["T1"],
        PALETTE["T2"],
        PALETTE["T3"],
    ]
    for index, ((column, label), color) in enumerate(zip(metrics, colors)):
        axes[3].bar(
            x + (index - 1.5) * width,
            panel_d[column],
            width=width,
            color=color,
            label=label,
        )
    axes[3].set_yscale("log")
    axes[3].set_xticks(x, panel_d["tier"])
    axes[3].set_ylabel("Count (log scale)")
    axes[3].set_title("Tier-specific QC and semantic compression")
    axes[3].legend(frameon=False, ncol=2)
    _clean(axes[3])
    _panel(axes[3], "D")
    recorder.add(
        "Fig1",
        "D",
        panel_d,
        source_path=tiers_path,
        description="Exclusive evidence-tier attrition and semantic compression.",
    )

    panel_e = final_sources.sort_values("final_grammar_rule_rows_supported")
    axes[4].barh(
        [str(value).replace("_", " ") for value in panel_e["source_database"]],
        pd.to_numeric(panel_e["final_grammar_rule_rows_supported"]),
        color=PALETTE["purple"],
    )
    axes[4].set_xlabel("Final grammar rule rows supported")
    axes[4].set_title("Evidence sources retained in the primary grammar")
    _clean(axes[4])
    _panel(axes[4], "E")
    recorder.add(
        "Fig1",
        "E",
        panel_e,
        source_path=final_sources_path,
        description="Overlapping source support in the final primary executable grammar.",
    )

    rule_usage["derivation_event_count"] = pd.to_numeric(
        rule_usage["derivation_event_count"], errors="coerce"
    ).fillna(0)
    active_rule_types = (
        rule_usage.groupby("reaction_type", dropna=False)
        .agg(
            G0_G3_active_rule_count=("grammar_rule_id", "nunique"),
            G0_G3_derivation_event_count=("derivation_event_count", "sum"),
        )
        .reset_index()
    )
    panel_f_all = reactions.merge(
        active_rule_types,
        on="reaction_type",
        how="left",
    )
    panel_f_all["G0_G3_active_rule_count"] = pd.to_numeric(
        panel_f_all["G0_G3_active_rule_count"], errors="coerce"
    ).fillna(0).astype(int)
    panel_f_all["G0_G3_derivation_event_count"] = pd.to_numeric(
        panel_f_all["G0_G3_derivation_event_count"], errors="coerce"
    ).fillna(0).astype(int)
    panel_f_all["final_grammar_rule_count"] = pd.to_numeric(
        panel_f_all["final_grammar_rule_count"], errors="coerce"
    ).fillna(0).astype(int)
    active_types = panel_f_all[
        panel_f_all["G0_G3_active_rule_count"] > 0
    ]
    inactive_types = panel_f_all[
        panel_f_all["G0_G3_active_rule_count"] == 0
    ].nlargest(3, "final_grammar_rule_count")
    displayed_types = set(active_types["reaction_type"]).union(
        inactive_types["reaction_type"]
    )
    panel_f_all["displayed_in_panel"] = panel_f_all["reaction_type"].isin(
        displayed_types
    )
    panel_f = panel_f_all[panel_f_all["displayed_in_panel"]].copy()
    panel_f["final_grammar_rule_count"] = pd.to_numeric(
        panel_f["final_grammar_rule_count"], errors="coerce"
    )
    panel_f = panel_f.sort_values(
        ["final_grammar_rule_count", "G0_G3_active_rule_count"]
    )
    y = np.arange(len(panel_f))
    axes[5].barh(
        y - 0.18,
        panel_f["final_grammar_rule_count"],
        height=0.34,
        color=PALETTE["gray"],
        label="candidate grammar",
    )
    axes[5].barh(
        y + 0.18,
        panel_f["G0_G3_active_rule_count"],
        height=0.34,
        color=PALETTE["T2"],
        label="activated in G0-G3",
    )
    axes[5].set_yticks(
        y,
        [
            _wrap(_reaction_type_label(value), 26)
            for value in panel_f["reaction_type"]
        ],
    )
    axes[5].set_xlabel("Primary grammar rules")
    axes[5].set_title("Candidate versus G0-G3-active grammar")
    axes[5].legend(frameon=False, fontsize=5.4)
    _clean(axes[5])
    _panel(axes[5], "F")
    recorder.add(
        "Fig1",
        "F",
        panel_f_all,
        source_path=f"{reactions_path};{rule_usage_path}",
        description=(
            "Reaction-type composition of the candidate primary grammar "
            "and the subset that produced at least one accepted G0-G3 event. "
            "All reaction types are retained in the source table; the panel "
            "shows every active type plus the three largest inactive types."
        ),
    )
    figure.tight_layout(pad=1.55, h_pad=2.1, w_pad=2.0)
    paths = _save(
        figure, output_dir, "Fig1_evidence_stratified_reaction_grammar"
    )
    plt.close(figure)
    return paths


def _benchmark_row(label: str, summary: dict[str, Any]) -> dict[str, Any]:
    positive_total = int(summary.get("positive_reactions", 0))
    return {
        "benchmark_mode": label,
        "connectivity_recovery_rate": float(
            summary.get("positive_connectivity_recovery_rate", 0) or 0
        ),
        "full_stereo_recovery_rate": (
            int(summary.get("positive_full_stereo_recovered", 0))
            / positive_total
            if positive_total
            else 0.0
        ),
        "decoy_connectivity_match_rate": float(
            summary.get("decoy_connectivity_match_rate", 0) or 0
        ),
        "positive_reactions": positive_total,
        "compiled_rules": int(summary.get("compiled_evaluated_rules", 0)),
    }


def _fig2(
    sensitivity_dir: Path,
    external_benchmark_dir: Path,
    domain_benchmark_dir: Path,
    output_dir: Path,
    recorder: SourceRecorder,
) -> list[Path]:
    plt = _plt()
    tier_path = sensitivity_dir / "G1_evidence_layer_sensitivity_summary.tsv"
    overlap_path = sensitivity_dir / "G1_pairwise_structure_overlap.tsv"
    membership_path = sensitivity_dir / "G1_structure_layer_membership.tsv"
    tier = read_table(tier_path)
    overlap = read_table(overlap_path)
    membership = read_table(membership_path)
    external_summary_path = external_benchmark_dir / "benchmark_summary.json"
    domain_summary_path = domain_benchmark_dir / "benchmark_summary.json"
    external_summary = _read_json(external_summary_path)
    domain_summary = _read_json(domain_summary_path)
    domain_positive_path = (
        domain_benchmark_dir / "benchmark_positive_recovery.tsv"
    )
    domain_positive = read_table(domain_positive_path)

    figure, axes = plt.subplots(2, 3, figsize=(10.2, 6.3))
    axes = axes.ravel()

    panel_a = pd.DataFrame(
        [
            _benchmark_row(
                "External-only\nleakage control", external_summary
            ),
            _benchmark_row(
                "Domain-informed\nreplay calibration", domain_summary
            ),
        ]
    )
    x = np.arange(len(panel_a))
    width = 0.24
    axes[0].bar(
        x - width,
        panel_a["connectivity_recovery_rate"],
        width,
        color=PALETTE["T1"],
        label="connectivity recovery",
    )
    axes[0].bar(
        x,
        panel_a["full_stereo_recovery_rate"],
        width,
        color=PALETTE["T2"],
        label="full-stereo recovery",
    )
    axes[0].bar(
        x + width,
        panel_a["decoy_connectivity_match_rate"],
        width,
        color=PALETTE["gray"],
        label="matched decoys",
    )
    axes[0].set_xticks(x, panel_a["benchmark_mode"])
    axes[0].set_ylim(0, 1.02)
    axes[0].set_ylabel("Fraction")
    axes[0].set_title("Specificity control and internal pathway replay")
    axes[0].legend(frameon=False, loc="upper left")
    _clean(axes[0])
    _panel(axes[0], "A")
    recorder.add(
        "Fig2",
        "A",
        panel_a,
        source_path=f"{external_summary_path};{domain_summary_path}",
        description="Leakage-controlled external recovery and non-independent domain replay calibration.",
    )

    panel_b = domain_positive[
        [
            "benchmark_id",
            "step_label",
            "recovered_connectivity",
            "recovered_full_stereo",
        ]
    ].copy()
    for column in ("recovered_connectivity", "recovered_full_stereo"):
        panel_b[column] = (
            panel_b[column].astype(str).str.lower().eq("true")
        )
    matrix = panel_b[
        ["recovered_connectivity", "recovered_full_stereo"]
    ].astype(int).to_numpy()
    axes[1].imshow(
        matrix,
        aspect="auto",
        cmap="Blues",
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    axes[1].set_xticks([0, 1], ["connectivity", "full stereo"])
    axes[1].set_yticks(
        np.arange(len(panel_b)),
        [f"{index:02d}" for index in range(1, len(panel_b) + 1)],
        fontsize=5.0,
    )
    axes[1].set_ylabel("Curated reaction index")
    axes[1].set_title("Per-reaction domain replay calibration")
    _panel(axes[1], "B")
    recorder.add(
        "Fig2",
        "B",
        panel_b,
        source_path=domain_positive_path,
        description="Per-reaction connectivity and full-stereochemistry replay status.",
    )

    for column in tier.columns:
        if column not in {"evidence_scope", "comparison_scope"}:
            converted = pd.to_numeric(tier[column], errors="coerce")
            if converted.notna().any():
                tier[column] = converted
    labels = ["T1 + domain", "T2", "T3"]
    x = np.arange(3)
    axes[2].bar(
        x,
        tier["unique_G1_structures"],
        color=[PALETTE["T1"], PALETTE["T2"], PALETTE["T3"]],
        width=0.62,
    )
    axes[2].set_xticks(x, labels)
    axes[2].set_ylabel("Unique G1 structures")
    twin = axes[2].twinx()
    twin.plot(
        x,
        tier["activated_rules_G1"],
        color=PALETTE["G0"],
        marker="o",
        linewidth=1.2,
    )
    twin.set_ylabel("G1-activated rules")
    axes[2].set_title("One-step space by exclusive evidence layer")
    _clean(axes[2])
    _panel(axes[2], "C")
    recorder.add(
        "Fig2",
        "C",
        tier,
        source_path=tier_path,
        description="G1 structure yield and activated-rule counts by exclusive evidence scope.",
    )

    scopes = list(tier["evidence_scope"])
    jaccard = np.eye(len(scopes), dtype=float)
    for row in overlap.to_dict("records"):
        left = scopes.index(row["left_evidence_scope"])
        right = scopes.index(row["right_evidence_scope"])
        value = float(row["jaccard_similarity"])
        jaccard[left, right] = value
        jaccard[right, left] = value
    image = axes[3].imshow(jaccard, vmin=0, vmax=1, cmap="viridis")
    short_scopes = ["T1 + domain", "T2", "T3"]
    axes[3].set_xticks(range(3), short_scopes)
    axes[3].set_yticks(range(3), short_scopes)
    for row in range(3):
        for column in range(3):
            axes[3].text(
                column,
                row,
                f"{jaccard[row, column]:.3f}",
                ha="center",
                va="center",
                color="white" if jaccard[row, column] < 0.5 else "black",
            )
    figure.colorbar(
        image,
        ax=axes[3],
        fraction=0.046,
        pad=0.04,
        label="Jaccard similarity",
    )
    axes[3].set_title("Full-stereochemistry G1 overlap")
    _panel(axes[3], "D")
    recorder.add(
        "Fig2",
        "D",
        overlap,
        source_path=overlap_path,
        description="Pairwise full-InChIKey G1 overlap between exclusive evidence layers.",
    )

    membership_columns = [
        column for column in membership if column.startswith("present_in_")
    ]
    membership_panel = (
        membership.groupby(membership_columns, dropna=False)
        .size()
        .reset_index(name="structure_count")
    )
    membership_panel["combination"] = membership_panel.apply(
        lambda row: " + ".join(
            label
            for label, column in zip(
                ["T1", "T2", "T3"], membership_columns
            )
            if str(row[column]).lower() == "true"
        ),
        axis=1,
    )
    membership_panel = membership_panel.sort_values(
        "structure_count", ascending=False
    )
    axes[4].bar(
        np.arange(len(membership_panel)),
        membership_panel["structure_count"],
        color=PALETTE["purple"],
    )
    axes[4].set_yscale("log")
    axes[4].set_xticks(
        np.arange(len(membership_panel)),
        membership_panel["combination"],
        rotation=40,
        ha="right",
    )
    axes[4].set_ylabel("G1 structures (log scale)")
    axes[4].set_title("Evidence-layer membership combinations")
    _clean(axes[4])
    _panel(axes[4], "E")
    recorder.add(
        "Fig2",
        "E",
        membership_panel,
        source_path=membership_path,
        description="Counts of G1 structures unique to or shared among evidence scopes.",
    )

    panel_f = tier[
        [
            "evidence_scope",
            "accepted_event_fraction_of_raw_products",
            "known_G0_full_recovery_events",
            "known_G0_connectivity_only_recovery_events",
        ]
    ].copy()
    axes[5].bar(
        x,
        panel_f["accepted_event_fraction_of_raw_products"],
        color=[PALETTE["T1"], PALETTE["T2"], PALETTE["T3"]],
    )
    axes[5].set_xticks(x, labels)
    axes[5].set_ylim(0, 1)
    axes[5].set_ylabel("Accepted events / raw products")
    axes[5].set_title("Product-level QC acceptance")
    _clean(axes[5])
    _panel(axes[5], "F")
    recorder.add(
        "Fig2",
        "F",
        panel_f,
        source_path=tier_path,
        description="Product-level acceptance and known-space recovery events by evidence scope.",
    )

    figure.tight_layout(pad=1.55, h_pad=2.1, w_pad=2.0)
    paths = _save(
        figure, output_dir, "Fig2_validation_and_evidence_sensitivity"
    )
    plt.close(figure)
    return paths


def _fig3(
    analysis_dir: Path,
    output_dir: Path,
    recorder: SourceRecorder,
) -> list[Path]:
    plt = _plt()
    generation_path = analysis_dir / "generation_expansion_summary.tsv"
    nearest_path = analysis_dir / "nearest_G0_similarity.tsv"
    descriptors_path = (
        analysis_dir / "physicochemical_descriptor_summary.tsv"
    )
    projection_path = (
        analysis_dir / "chemical_space_fingerprint_projection.tsv"
    )
    generation = read_table(generation_path)
    nearest = read_table(nearest_path)
    descriptors = read_table(descriptors_path)
    projection = read_table(projection_path)
    for column in generation.columns:
        if column != "interpretation_layer":
            generation[column] = pd.to_numeric(
                generation[column], errors="coerce"
            )
    nearest["generation"] = pd.to_numeric(
        nearest["generation"], errors="coerce"
    ).astype(int)
    nearest["nearest_G0_tanimoto"] = pd.to_numeric(
        nearest["nearest_G0_tanimoto"], errors="coerce"
    )
    for column in ("generation", "mean", "standard_deviation"):
        descriptors[column] = pd.to_numeric(
            descriptors[column], errors="coerce"
        )
    for column in ("generation", "molecular_fp_axis_1", "molecular_fp_axis_2"):
        projection[column] = pd.to_numeric(
            projection[column], errors="coerce"
        )

    figure, axes = plt.subplots(2, 3, figsize=(10.2, 6.4))
    axes = axes.ravel()
    x = generation["generation"].astype(int).to_numpy()
    colors = [PALETTE[f"G{value}"] for value in x]

    axes[0].bar(
        x,
        generation["unique_nodes_first_observed"],
        color=colors,
        width=0.68,
    )
    axes[0].set_yscale("log")
    axes[0].set_xticks(x, [f"G{value}" for value in x])
    axes[0].set_ylabel("Unique structures (log scale)")
    axes[0].set_title("Generation-wise expansion")
    _clean(axes[0])
    _panel(axes[0], "A")
    recorder.add(
        "Fig3",
        "A",
        generation,
        source_path=generation_path,
        description="Generation-wise unique structure and event inventory.",
    )

    axes[1].plot(
        x,
        generation["cumulative_unique_nodes"],
        marker="o",
        color=PALETTE["purple"],
    )
    axes[1].set_yscale("log")
    axes[1].set_xticks(x, [f"G{value}" for value in x])
    axes[1].set_ylabel("Cumulative structures (log scale)")
    axes[1].set_title("Cumulative accessible chemical space")
    _clean(axes[1])
    _panel(axes[1], "B")
    recorder.add(
        "Fig3",
        "B",
        generation[
            ["generation", "cumulative_unique_nodes", "interpretation_layer"]
        ],
        source_path=generation_path,
        description="Cumulative full-stereochemistry-aware chemical-space size.",
    )

    nonzero = generation[generation["generation"] > 0]
    gx = nonzero["generation"].astype(int).to_numpy()
    width = 0.23
    axes[2].bar(
        gx - width,
        nonzero["raw_product_tuples"],
        width,
        color=PALETTE["gray"],
        label="raw products",
    )
    axes[2].bar(
        gx,
        nonzero["rejected_product_events"],
        width,
        color=PALETTE["red"],
        label="rejected",
    )
    axes[2].bar(
        gx + width,
        nonzero["derivation_events"],
        width,
        color=PALETTE["T2"],
        label="accepted",
    )
    axes[2].set_yscale("log")
    axes[2].set_xticks(gx, [f"G{value}" for value in gx])
    axes[2].set_ylabel("Events (log scale)")
    axes[2].set_title("Enumeration and product-level QC")
    axes[2].legend(frameon=False)
    _clean(axes[2])
    _panel(axes[2], "C")
    recorder.add(
        "Fig3",
        "C",
        nonzero,
        source_path=generation_path,
        description="Raw, rejected, and accepted products at each generated layer.",
    )

    axes[3].plot(
        gx,
        nonzero["unique_node_yield_per_raw_product"],
        marker="o",
        color=PALETTE["T1"],
        label="unique-node yield",
    )
    axes[3].plot(
        gx,
        nonzero["known_connectivity_reconnection_fraction"],
        marker="s",
        color=PALETTE["T3"],
        label="known-space reconnection",
    )
    axes[3].set_xticks(gx, [f"G{value}" for value in gx])
    axes[3].set_ylim(bottom=0)
    axes[3].set_ylabel("Fraction")
    axes[3].set_title("Novelty yield and known-space reconnection")
    axes[3].legend(frameon=False)
    _clean(axes[3])
    _panel(axes[3], "D")
    recorder.add(
        "Fig3",
        "D",
        nonzero[
            [
                "generation",
                "unique_node_yield_per_raw_product",
                "known_connectivity_reconnection_fraction",
                "interpretation_layer",
            ]
        ],
        source_path=generation_path,
        description="Novel structure yield and connectivity-level G0 reconnection fractions.",
    )

    generations = sorted(nearest["generation"].unique())
    values = [
        nearest.loc[
            nearest["generation"] == generation_value,
            "nearest_G0_tanimoto",
        ]
        .dropna()
        .to_numpy()
        for generation_value in generations
    ]
    violin = axes[4].violinplot(
        values, positions=generations, showmedians=True, widths=0.75
    )
    for body, generation_value in zip(violin["bodies"], generations):
        body.set_facecolor(PALETTE[f"G{generation_value}"])
        body.set_edgecolor("none")
        body.set_alpha(0.75)
    axes[4].set_xticks(
        generations, [f"G{value}" for value in generations]
    )
    axes[4].set_ylim(0, 1.02)
    axes[4].set_ylabel("Nearest-G0 Morgan Tanimoto")
    axes[4].set_title("Displacement from known taxane space")
    _clean(axes[4])
    _panel(axes[4], "E")
    recorder.add(
        "Fig3",
        "E",
        nearest,
        source_path=nearest_path,
        description="Deterministically sampled nearest-G0 fingerprint similarities with sampling metadata.",
    )

    descriptor_names = [
        "exact_mass",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "fraction_csp3",
    ]
    descriptor_panel = descriptors[
        descriptors["descriptor"].isin(descriptor_names)
    ].copy()
    g0_stats = descriptor_panel[
        descriptor_panel["generation"] == 0
    ].set_index("descriptor")
    generation_order = sorted(descriptor_panel["generation"].unique())
    heat = []
    for generation_value in generation_order:
        group = descriptor_panel[
            descriptor_panel["generation"] == generation_value
        ].set_index("descriptor")
        heat.append(
            [
                (group.loc[name, "mean"] - g0_stats.loc[name, "mean"])
                / max(g0_stats.loc[name, "standard_deviation"], 1e-12)
                for name in descriptor_names
            ]
        )
    heat = np.asarray(heat)
    limit = max(1.0, float(np.nanmax(np.abs(heat))))
    image = axes[5].imshow(
        heat, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit
    )
    axes[5].set_xticks(
        range(len(descriptor_names)),
        [_descriptor_label(value) for value in descriptor_names],
        rotation=40,
        ha="right",
    )
    axes[5].set_yticks(
        range(len(heat)),
        [f"G{int(value)}" for value in generation_order],
    )
    axes[5].set_title("Physicochemical displacement from G0")
    figure.colorbar(
        image,
        ax=axes[5],
        fraction=0.046,
        pad=0.04,
        label="Mean shift (G0 SD)",
    )
    _panel(axes[5], "F")
    recorder.add(
        "Fig3",
        "F",
        descriptor_panel,
        source_path=descriptors_path,
        description="Generation-level physicochemical summary used for G0-standardized mean shifts.",
    )

    recorder.add(
        "Fig3",
        "S1",
        projection,
        source_path=projection_path,
        description="Morgan fingerprint TruncatedSVD coordinates and deterministic sampling metadata; source data for supplementary projection figure.",
    )
    figure.tight_layout(pad=1.55, h_pad=2.1, w_pad=2.0)
    paths = _save(
        figure, output_dir, "Fig3_iterative_taxane_space_expansion"
    )
    plt.close(figure)
    return paths


def _fig4(
    analysis_dir: Path,
    output_dir: Path,
    recorder: SourceRecorder,
) -> list[Path]:
    plt = _plt()
    transition_path = (
        analysis_dir / "functional_state_transition_summary.tsv"
    )
    rule_path = analysis_dir / "reaction_grammar_usage.tsv"
    concentration_path = (
        analysis_dir / "reaction_grammar_use_concentration.tsv"
    )
    edit_path = analysis_dir / "reaction_edit_landscape.tsv"
    transition = read_table(transition_path)
    rule = read_table(rule_path)
    concentration = read_table(concentration_path)
    edit = read_table(edit_path)
    for frame in (transition, rule, concentration, edit):
        for column in frame.columns:
            if column.endswith("_count") or column in {
                "generation",
                "generation_min",
                "generation_max",
                "derivation_events",
                "rules_used",
                "semantic_groups_used",
                "top1_rule_event_fraction",
                "top5_rule_event_fraction",
                "top10_rule_event_fraction",
                "effective_rule_number_exp_shannon",
                "observed_changed_source_atoms",
            }:
                frame[column] = pd.to_numeric(
                    frame[column], errors="coerce"
                )

    figure, axes = plt.subplots(2, 3, figsize=(10.2, 6.4))
    axes = axes.ravel()

    top_transitions = (
        transition.groupby("functional_state_transition")[
            "derivation_event_count"
        ]
        .sum()
        .nlargest(8)
        .index
    )
    panel_a = transition[
        transition["functional_state_transition"].isin(top_transitions)
    ].copy()
    pivot = panel_a.pivot_table(
        index="functional_state_transition",
        columns="generation",
        values="derivation_event_count",
        aggfunc="sum",
        fill_value=0,
    )
    pivot = pivot.loc[pivot.sum(axis=1).sort_values().index]
    image = axes[0].imshow(
        np.log10(pivot.to_numpy() + 1), aspect="auto", cmap="magma"
    )
    axes[0].set_yticks(
        range(len(pivot)),
        [
            _functional_state_transition_label(value)
            for value in pivot.index
        ],
    )
    axes[0].set_xticks(
        range(len(pivot.columns)),
        [f"G{int(value)}" for value in pivot.columns],
    )
    axes[0].set_title("Dominant functional-state transitions")
    figure.colorbar(
        image,
        ax=axes[0],
        fraction=0.046,
        pad=0.04,
        label="log10(events + 1)",
    )
    _panel(axes[0], "A")
    recorder.add(
        "Fig4",
        "A",
        panel_a,
        source_path=transition_path,
        description=(
            "Top functional-state transition patterns by generation; raw "
            "transition strings are retained in the source table."
        ),
    )

    top_rules = rule.nlargest(12, "derivation_event_count").copy()
    axes[1].barh(
        np.arange(len(top_rules)),
        top_rules["derivation_event_count"].iloc[::-1],
        color=PALETTE["T1"],
    )
    axes[1].set_yticks(
        np.arange(len(top_rules)),
        [
            _rule_label(row)
            for _, row in top_rules.iloc[::-1].iterrows()
        ],
    )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Derivation events (log scale)")
    axes[1].set_title("Most-used grammar productions")
    _clean(axes[1])
    _panel(axes[1], "B")
    recorder.add(
        "Fig4",
        "B",
        top_rules,
        source_path=rule_path,
        description="Top grammar productions ranked by accepted derivation events.",
    )

    by_type = (
        rule.groupby("reaction_type", dropna=False)["derivation_event_count"]
        .sum()
        .sort_values()
        .tail(12)
        .reset_index()
    )
    axes[2].barh(
        [
            _wrap(_reaction_type_label(value), 25)
            for value in by_type["reaction_type"]
        ],
        by_type["derivation_event_count"],
        color=PALETTE["T2"],
    )
    axes[2].set_xscale("log")
    axes[2].set_xlabel("Derivation events (log scale)")
    axes[2].set_title("Transformation semantics of generated space")
    _clean(axes[2])
    _panel(axes[2], "C")
    recorder.add(
        "Fig4",
        "C",
        by_type,
        source_path=rule_path,
        description="Accepted derivation events aggregated by reaction type.",
    )

    per_generation = concentration[
        concentration["generation_scope"] != "all_generations"
    ].copy()
    gx = (
        per_generation["generation_scope"]
        .str.replace("G", "", regex=False)
        .astype(int)
        .to_numpy()
    )
    axes[3].plot(
        gx,
        per_generation["top1_rule_event_fraction"],
        marker="o",
        label="top 1 rule",
        color=PALETTE["red"],
    )
    axes[3].plot(
        gx,
        per_generation["top5_rule_event_fraction"],
        marker="s",
        label="top 5 rules",
        color=PALETTE["purple"],
    )
    axes[3].plot(
        gx,
        per_generation["top10_rule_event_fraction"],
        marker="^",
        label="top 10 rules",
        color=PALETTE["gray"],
    )
    axes[3].set_ylim(0, 1.02)
    axes[3].set_xticks(gx, [f"G{value}" for value in gx])
    axes[3].set_ylabel("Cumulative event fraction")
    axes[3].set_title("Grammar-use concentration")
    axes[3].legend(frameon=False)
    _clean(axes[3])
    _panel(axes[3], "D")
    recorder.add(
        "Fig4",
        "D",
        concentration,
        source_path=concentration_path,
        description="Shannon, HHI, Gini, and top-k grammar-use concentration metrics.",
    )

    edit_panel = (
        edit.groupby(
            ["generation", "observed_changed_source_atoms"], dropna=False
        )["derivation_event_count"]
        .sum()
        .reset_index()
    )
    edit_panel = edit_panel[
        edit_panel["observed_changed_source_atoms"] <= 12
    ]
    pivot_edit = edit_panel.pivot_table(
        index="observed_changed_source_atoms",
        columns="generation",
        values="derivation_event_count",
        fill_value=0,
    )
    image = axes[4].imshow(
        np.log10(pivot_edit.to_numpy() + 1),
        aspect="auto",
        cmap="viridis",
    )
    axes[4].set_yticks(
        range(len(pivot_edit)),
        [str(int(value)) for value in pivot_edit.index],
    )
    axes[4].set_xticks(
        range(len(pivot_edit.columns)),
        [f"G{int(value)}" for value in pivot_edit.columns],
    )
    axes[4].set_ylabel("Changed source atoms")
    axes[4].set_title("Reaction-edit locality")
    figure.colorbar(
        image,
        ax=axes[4],
        fraction=0.046,
        pad=0.04,
        label="log10(events + 1)",
    )
    _panel(axes[4], "E")
    recorder.add(
        "Fig4",
        "E",
        edit_panel,
        source_path=edit_path,
        description="Accepted events by generation and observed changed source-atom count.",
    )

    delta_panel = (
        edit.groupby(["observed_element_delta"], dropna=False)[
            "derivation_event_count"
        ]
        .sum()
        .nlargest(12)
        .sort_values()
        .reset_index()
    )
    delta_panel["elemental_edit_signature"] = delta_panel[
        "observed_element_delta"
    ].map(_format_element_delta)
    axes[5].barh(
        delta_panel["elemental_edit_signature"],
        delta_panel["derivation_event_count"],
        color=PALETTE["T3"],
    )
    axes[5].set_xscale("log")
    axes[5].set_xlabel("Derivation events (log scale)")
    axes[5].set_title("Observed elemental edit signatures")
    _clean(axes[5])
    _panel(axes[5], "F")
    recorder.add(
        "Fig4",
        "F",
        delta_panel,
        source_path=edit_path,
        description="Most frequent observed molecular formula deltas.",
    )

    figure.tight_layout(pad=1.55, h_pad=2.1, w_pad=2.0)
    paths = _save(
        figure, output_dir, "Fig4_reaction_grammar_transformation_landscape"
    )
    plt.close(figure)
    return paths


def _fig5(
    analysis_dir: Path,
    output_dir: Path,
    recorder: SourceRecorder,
) -> list[Path]:
    plt = _plt()
    convergence_path = (
        analysis_dir / "convergence_and_route_multiplicity.tsv"
    )
    bridges_path = analysis_dir / "latent_bridge_candidates.tsv"
    pair_path = analysis_dir / "known_G0_pair_bridge_summary.tsv"
    convergence = read_table(convergence_path)
    bridges = read_table(bridges_path)
    if pair_path.exists() and pair_path.stat().st_size > 1:
        try:
            pair_frame = read_table(pair_path)
        except pd.errors.EmptyDataError:
            pair_frame = pd.DataFrame()
    else:
        pair_frame = pd.DataFrame()
    for column in (
        "generation",
        "incoming_derivation_events",
        "unique_parent_count",
        "unique_rule_count",
        "unique_semantic_group_count",
        "structural_path_count",
        "log10_structural_path_count",
        "semantic_edge_path_count",
        "log10_semantic_edge_path_count",
        "raw_rule_event_path_count",
        "log10_raw_rule_event_path_count",
    ):
        convergence[column] = pd.to_numeric(
            convergence[column], errors="coerce"
        )
    for column in (
        "generation",
        "known_G0_ancestor_count",
        "known_G0_descendant_count",
        "distinct_G0_pair_bridge_count",
    ):
        bridges[column] = pd.to_numeric(bridges[column], errors="coerce")
    bridges["latent_bridge_candidate"] = (
        bridges["latent_bridge_candidate"].astype(str).str.lower().eq("true")
    )

    figure, axes = plt.subplots(2, 3, figsize=(10.2, 6.4))
    axes = axes.ravel()

    convergence_bool = (
        convergence["is_convergent"].astype(str).str.lower().eq("true")
    )
    generated = convergence[convergence["generation"] > 0].copy()
    generated_convergence_bool = convergence_bool.loc[generated.index]
    panel_a = (
        generated.assign(_convergent=generated_convergence_bool)
        .groupby("generation")
        .agg(
            structures=("space_id", "size"),
            convergent_structures=("_convergent", "sum"),
        )
        .reset_index()
    )
    panel_a["convergent_fraction"] = (
        panel_a["convergent_structures"] / panel_a["structures"]
    )
    axes[0].bar(
        panel_a["generation"],
        panel_a["convergent_fraction"],
        color=[
            PALETTE[f"G{int(value)}"] for value in panel_a["generation"]
        ],
    )
    axes[0].set_xticks(
        panel_a["generation"],
        [f"G{int(value)}" for value in panel_a["generation"]],
    )
    axes[0].set_ylim(
        0, max(0.05, panel_a["convergent_fraction"].max() * 1.2)
    )
    axes[0].set_ylabel("Convergent structures / all structures")
    axes[0].set_title("Convergence by generation")
    _clean(axes[0])
    _panel(axes[0], "A")
    recorder.add(
        "Fig5",
        "A",
        panel_a,
        source_path=convergence_path,
        description=(
            "Generated structures with at least two distinct parent "
            "structures in their first-observed generation."
        ),
    )

    generations = sorted(
        generated["generation"].dropna().astype(int).unique()
    )
    parent_values = [
        generated.loc[
            generated["generation"] == value, "unique_parent_count"
        ]
        .clip(upper=20)
        .to_numpy()
        for value in generations
    ]
    axes[1].boxplot(
        parent_values, positions=generations, showfliers=False
    )
    axes[1].set_xticks(
        generations, [f"G{value}" for value in generations]
    )
    axes[1].set_ylabel("Unique immediate parents (capped at 20)")
    axes[1].set_title("Immediate route multiplicity")
    _clean(axes[1])
    _panel(axes[1], "B")
    recorder.add(
        "Fig5",
        "B",
        generated[
            [
                "space_id",
                "generation",
                "unique_parent_count",
                "unique_rule_count",
                "unique_semantic_group_count",
                "incoming_derivation_events",
            ]
        ],
        source_path=convergence_path,
        description="Per-structure immediate parent multiplicity with semantic-group, rule, and event support retained for audit.",
    )

    path_values = [
        generated.loc[
            generated["generation"] == value,
            "log10_structural_path_count",
        ]
        .dropna()
        .to_numpy()
        for value in generations
    ]
    violin = axes[2].violinplot(
        path_values,
        positions=generations,
        showmedians=True,
        widths=0.75,
    )
    for body, generation_value in zip(violin["bodies"], generations):
        body.set_facecolor(PALETTE[f"G{generation_value}"])
        body.set_alpha(0.75)
    axes[2].set_xticks(
        generations, [f"G{value}" for value in generations]
    )
    axes[2].set_ylabel("log10(distinct structure paths)")
    axes[2].set_title("Multi-generation route multiplicity")
    _clean(axes[2])
    _panel(axes[2], "C")
    recorder.add(
        "Fig5",
        "C",
        generated[
            [
                "space_id",
                "generation",
                "structural_path_count",
                "log10_structural_path_count",
                "semantic_edge_path_count",
                "log10_semantic_edge_path_count",
                "raw_rule_event_path_count",
                "log10_raw_rule_event_path_count",
            ]
        ],
        source_path=convergence_path,
        description="Dynamic-programming path counts over distinct structure edges, with semantic-edge and raw rule-event counts retained for audit.",
    )

    panel_d = (
        bridges.groupby("generation")
        .agg(
            generated_structures=("space_id", "size"),
            latent_bridge_candidates=("latent_bridge_candidate", "sum"),
        )
        .reset_index()
    )
    panel_d["latent_bridge_fraction"] = (
        panel_d["latent_bridge_candidates"] / panel_d["generated_structures"]
    )
    panel_d["bridge_status"] = np.where(
        panel_d["generation"] < panel_d["generation"].max(),
        "evaluated",
        "right_censored_no_subsequent_parent_expansion",
    )
    evaluated_panel_d = panel_d[panel_d["bridge_status"] == "evaluated"]
    axes[3].bar(
        evaluated_panel_d["generation"],
        evaluated_panel_d["latent_bridge_candidates"],
        color=[
            PALETTE[f"G{int(value)}"]
            for value in evaluated_panel_d["generation"]
        ],
    )
    axes[3].set_yscale("symlog", linthresh=1)
    generation_ticks = panel_d["generation"].astype(int).to_list()
    axes[3].set_xticks(
        generation_ticks,
        [f"G{value}" for value in generation_ticks],
    )
    maximum_bridge_count = max(
        1, int(evaluated_panel_d["latent_bridge_candidates"].max())
    )
    axes[3].text(
        max(generation_ticks),
        max(1.2, maximum_bridge_count ** 0.25),
        "right-censored",
        ha="center",
        va="bottom",
        fontsize=6.2,
        color=PALETTE["gray"],
    )
    axes[3].set_xlim(min(generation_ticks) - 0.5, max(generation_ticks) + 0.5)
    axes[3].set_ylabel("Latent bridge candidates (symlog)")
    axes[3].set_title("Latent bridge candidates")
    _clean(axes[3])
    _panel(axes[3], "D")
    recorder.add(
        "Fig5",
        "D",
        panel_d,
        source_path=bridges_path,
        description=(
            "Latent bridge counts and fractions by generation; G3 is "
            "right-censored because it was not used as a parent layer."
        ),
    )

    if not pair_frame.empty:
        pair_frame["supporting_latent_intermediate_count"] = pd.to_numeric(
            pair_frame["supporting_latent_intermediate_count"],
            errors="coerce",
        )
        top_pairs = pair_frame.nlargest(
            20, "supporting_latent_intermediate_count"
        ).copy()
        top_pairs["pair"] = (
            top_pairs["known_G0_source_space_id"]
            + " → "
            + top_pairs["known_G0_target_space_id"]
        )
        axes[4].barh(
            np.arange(len(top_pairs)),
            top_pairs["supporting_latent_intermediate_count"].iloc[::-1],
            color=PALETTE["purple"],
        )
        axes[4].set_yticks(
            np.arange(len(top_pairs)),
            top_pairs["pair"].iloc[::-1],
            fontsize=5.2,
        )
        axes[4].set_xlabel("Supporting latent intermediates")
    else:
        top_pairs = pair_frame
        axes[4].text(
            0.5,
            0.5,
            "No G0-pair bridges at this depth",
            ha="center",
            va="center",
            transform=axes[4].transAxes,
        )
    axes[4].set_title("Most-supported G0 pairs")
    _clean(axes[4])
    _panel(axes[4], "E")
    recorder.add(
        "Fig5",
        "E",
        top_pairs,
        source_path=pair_path,
        description="Top directed known-G0 pairs ranked by supporting latent intermediates.",
    )

    bridge_candidates = bridges[bridges["latent_bridge_candidate"]].copy()
    if not bridge_candidates.empty:
        maximum = int(
            bridge_candidates["distinct_G0_pair_bridge_count"].max()
        )
        bins = np.arange(1, min(50, maximum) + 2)
        for generation_value in sorted(
            bridge_candidates["generation"].unique()
        ):
            values = bridge_candidates.loc[
                bridge_candidates["generation"] == generation_value,
                "distinct_G0_pair_bridge_count",
            ]
            axes[5].hist(
                values,
                bins=bins,
                histtype="step",
                linewidth=1.2,
                color=PALETTE[f"G{int(generation_value)}"],
                label=f"G{int(generation_value)}",
            )
        axes[5].set_yscale("log")
        axes[5].legend(frameon=False)
    axes[5].set_xlabel("Distinct directed G0 pairs per intermediate")
    axes[5].set_ylabel("Intermediate count (log scale)")
    axes[5].set_title("Bridge breadth")
    _clean(axes[5])
    _panel(axes[5], "F")
    recorder.add(
        "Fig5",
        "F",
        bridge_candidates,
        source_path=bridges_path,
        description="Known-G0 pair breadth of every latent bridge candidate.",
    )

    figure.tight_layout(pad=1.55, h_pad=2.1, w_pad=2.0)
    paths = _save(
        figure, output_dir, "Fig5_convergence_and_latent_bridges"
    )
    plt.close(figure)
    return paths


def _supplement_figures(
    analysis_dir: Path,
    t2_analysis_dir: Path,
    t3_analysis_dir: Path,
    output_dir: Path,
    recorder: SourceRecorder,
) -> dict[str, list[Path]]:
    plt = _plt()
    outputs: dict[str, list[Path]] = {}

    projection_path = (
        analysis_dir / "chemical_space_fingerprint_projection.tsv"
    )
    projection = read_table(projection_path)
    for column in ("generation", "molecular_fp_axis_1", "molecular_fp_axis_2"):
        projection[column] = pd.to_numeric(
            projection[column], errors="coerce"
        )
    figure, axis = plt.subplots(figsize=(6.4, 5.2))
    for generation in sorted(
        projection["generation"].dropna().unique(), reverse=True
    ):
        group = projection[projection["generation"] == generation]
        generation = int(generation)
        axis.scatter(
            group["molecular_fp_axis_1"],
            group["molecular_fp_axis_2"],
            s=4 if generation else 8,
            alpha=0.32 if generation else 0.75,
            linewidths=0,
            color=PALETTE[f"G{generation}"],
            label=f"G{generation} (n={len(group):,})",
            rasterized=True,
        )
    variance_1 = float(
        pd.to_numeric(
            projection["axis_1_explained_variance_ratio"],
            errors="coerce",
        ).dropna().iloc[0]
    )
    variance_2 = float(
        pd.to_numeric(
            projection["axis_2_explained_variance_ratio"],
            errors="coerce",
        ).dropna().iloc[0]
    )
    axis.set_xlabel(f"Fingerprint SVD axis 1 ({variance_1:.1%})")
    axis.set_ylabel(f"Fingerprint SVD axis 2 ({variance_2:.1%})")
    axis.set_title("Generation-stratified Morgan-fingerprint projection")
    axis.legend(frameon=False, markerscale=2.5)
    _clean(axis, grid=False)
    recorder.add(
        "FigS1",
        "A",
        projection,
        source_path=projection_path,
        description=(
            "Deterministically sampled Morgan-fingerprint TruncatedSVD "
            "projection with generation population and sample metadata."
        ),
    )
    figure.tight_layout(pad=1.2)
    outputs["FigS1"] = _save(
        figure, output_dir, "FigS1_fingerprint_projection"
    )
    plt.close(figure)

    descriptor_path = (
        analysis_dir / "physicochemical_descriptor_summary.tsv"
    )
    descriptors = read_table(descriptor_path)
    numeric_columns = [
        "generation",
        "median",
        "q1",
        "q3",
        "minimum",
        "maximum",
    ]
    for column in numeric_columns:
        descriptors[column] = pd.to_numeric(
            descriptors[column], errors="coerce"
        )
    descriptor_order = [
        "exact_mass",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "fraction_csp3",
        "ring_count",
    ]
    figure, axes = plt.subplots(2, 4, figsize=(10.2, 5.4))
    for axis, descriptor in zip(axes.ravel(), descriptor_order):
        panel = descriptors[descriptors["descriptor"] == descriptor].copy()
        panel = panel.sort_values("generation")
        x = panel["generation"].to_numpy(dtype=float)
        axis.fill_between(
            x,
            panel["q1"].to_numpy(dtype=float),
            panel["q3"].to_numpy(dtype=float),
            color=PALETTE["light"],
            alpha=0.75,
            label="IQR",
        )
        axis.plot(
            x,
            panel["median"],
            color=PALETTE["T1"],
            marker="o",
            linewidth=1.2,
            label="median",
        )
        axis.set_xticks(x, [f"G{int(value)}" for value in x])
        axis.set_title(_descriptor_label(descriptor))
        _clean(axis)
    axes[0, 0].legend(frameon=False)
    figure.suptitle(
        "Generation-specific physicochemical descriptor summaries",
        fontsize=10,
        y=1.01,
    )
    recorder.add(
        "FigS2",
        "A",
        descriptors[descriptors["descriptor"].isin(descriptor_order)],
        source_path=descriptor_path,
        description=(
            "Absolute medians and interquartile ranges of molecular "
            "descriptors by generation."
        ),
    )
    figure.tight_layout(pad=1.2, h_pad=1.6, w_pad=1.4)
    outputs["FigS2"] = _save(
        figure, output_dir, "FigS2_physicochemical_distributions"
    )
    plt.close(figure)

    rule_generation_path = (
        analysis_dir / "functional_state_transition_rule_summary.tsv"
    )
    rule_generation = read_table(rule_generation_path)
    rule_generation["generation"] = pd.to_numeric(
        rule_generation["generation"], errors="coerce"
    ).astype(int)
    rule_generation["derivation_event_count"] = pd.to_numeric(
        rule_generation["derivation_event_count"], errors="coerce"
    )
    rule_by_generation = (
        rule_generation.groupby(["generation", "grammar_rule_id"])
        ["derivation_event_count"]
        .sum()
        .reset_index()
    )
    active_summary = (
        rule_by_generation.groupby("generation")
        .agg(
            active_rules=("grammar_rule_id", "nunique"),
            derivation_events=("derivation_event_count", "sum"),
        )
        .reset_index()
    )
    top_rule_ids = (
        rule_by_generation.groupby("grammar_rule_id")[
            "derivation_event_count"
        ]
        .sum()
        .nlargest(20)
        .index
    )
    top_rules = rule_generation[
        rule_generation["grammar_rule_id"].isin(top_rule_ids)
    ].copy()
    top_rules = (
        top_rules.groupby(
            ["generation", "grammar_rule_id", "reaction_type"],
            dropna=False,
        )["derivation_event_count"]
        .sum()
        .reset_index()
    )
    rule_labels = {
        row["grammar_rule_id"]: _rule_label(row)
        for _, row in top_rules.drop_duplicates(
            "grammar_rule_id"
        ).iterrows()
    }
    pivot = top_rules.pivot_table(
        index="grammar_rule_id",
        columns="generation",
        values="derivation_event_count",
        fill_value=0,
    )
    pivot = pivot.loc[pivot.sum(axis=1).sort_values().index]
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 7.2))
    x = active_summary["generation"].to_numpy()
    axes[0].bar(
        x,
        active_summary["active_rules"],
        color=[PALETTE[f"G{int(value)}"] for value in x],
    )
    axes[0].set_xticks(x, [f"G{int(value)}" for value in x])
    axes[0].set_ylabel("Activated grammar rules")
    axes[0].set_title("Rule activation by generation")
    _clean(axes[0])
    _panel(axes[0], "A")
    image = axes[1].imshow(
        np.log10(pivot.to_numpy() + 1),
        aspect="auto",
        cmap="magma",
    )
    axes[1].set_yticks(
        range(len(pivot)),
        [rule_labels[value] for value in pivot.index],
        fontsize=4.8,
    )
    axes[1].set_xticks(
        range(len(pivot.columns)),
        [f"G{int(value)}" for value in pivot.columns],
    )
    displayed_rule_count = len(pivot)
    axes[1].set_title(
        f"Generation usage of the {displayed_rule_count} most active rules"
    )
    figure.colorbar(
        image,
        ax=axes[1],
        fraction=0.046,
        pad=0.04,
        label="log10(events + 1)",
    )
    _panel(axes[1], "B")
    recorder.add(
        "FigS3",
        "A",
        active_summary,
        source_path=rule_generation_path,
        description="Activated exact rules and accepted events by generation.",
    )
    recorder.add(
        "FigS3",
        "B",
        top_rules,
        source_path=rule_generation_path,
        description=(
            "Generation-resolved event counts for the "
            f"{displayed_rule_count} most-used rules."
        ),
    )
    figure.tight_layout(pad=1.25, w_pad=2.0)
    outputs["FigS3"] = _save(
        figure, output_dir, "FigS3_rule_activation_across_generations"
    )
    plt.close(figure)

    rejection_path = analysis_dir / "generation_rejection_summary.tsv"
    rejections = read_table(rejection_path)
    rejections["generation"] = pd.to_numeric(
        rejections["generation"], errors="coerce"
    ).astype(int)
    rejections["rejected_product_count"] = pd.to_numeric(
        rejections["rejected_product_count"], errors="coerce"
    )
    rejections["rejection_class"] = rejections["rejection_reason"].map(
        lambda value: (
            "invalid valence"
            if str(value).startswith("Explicit valence")
            else str(value).replace("_", " ")
        )
    )
    rejection_classes = (
        rejections.groupby(["generation", "rejection_class"], as_index=False)[
            "rejected_product_count"
        ].sum()
    )
    rejection_pivot = rejection_classes.pivot_table(
        index="generation",
        columns="rejection_class",
        values="rejected_product_count",
        fill_value=0,
    )
    rejection_fraction = rejection_pivot.div(
        rejection_pivot.sum(axis=1).replace(0, np.nan), axis=0
    ).fillna(0)
    figure, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    rejection_pivot.plot.bar(
        stacked=True, ax=axes[0], colormap="tab20", width=0.72
    )
    rejection_fraction.plot.bar(
        stacked=True, ax=axes[1], colormap="tab20", width=0.72
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Rejected products (log scale)")
    axes[0].set_title("Absolute product rejection counts")
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Fraction of rejected products")
    axes[1].set_title("Composition of product rejection")
    handles, labels = axes[0].get_legend_handles_labels()
    for axis in axes:
        axis.set_xlabel("Generation")
        axis.set_xticklabels(
            [f"G{int(value)}" for value in rejection_pivot.index],
            rotation=0,
        )
        legend = axis.get_legend()
        if legend is not None:
            legend.remove()
        _clean(axis)
    figure.legend(
        handles,
        labels,
        frameon=False,
        ncol=min(3, len(labels)),
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
    )
    _panel(axes[0], "A")
    _panel(axes[1], "B")
    recorder.add(
        "FigS4",
        "A",
        rejections,
        source_path=rejection_path,
        description=(
            "Absolute rejected-product counts by raw reason and normalized "
            "QC class. Atom-specific RDKit valence messages are grouped as "
            "invalid valence for visualization."
        ),
    )
    recorder.add(
        "FigS4",
        "B",
        rejection_fraction.reset_index(),
        source_path=rejection_path,
        description="Within-generation fractions of product rejection reasons.",
    )
    figure.tight_layout(pad=1.2, w_pad=1.8, rect=(0, 0.11, 1, 1))
    outputs["FigS4"] = _save(
        figure, output_dir, "FigS4_product_rejection_reasons"
    )
    plt.close(figure)

    generation_path = analysis_dir / "generation_expansion_summary.tsv"
    generation = read_table(generation_path)
    for column in generation.columns:
        if column not in {"interpretation_layer"}:
            generation[column] = pd.to_numeric(
                generation[column], errors="coerce"
            )
    generated = generation[generation["generation"] > 0].copy()
    gx = generated["generation"].astype(int).to_numpy()
    figure, axes = plt.subplots(1, 2, figsize=(8.6, 3.7))
    axes[0].plot(
        gx,
        generated["immediate_reverse_cycle_events"],
        marker="o",
        color=PALETTE["red"],
        label="immediate reverse cycles",
    )
    axes[0].plot(
        gx,
        generated["known_G0_connectivity_recovery_events"],
        marker="s",
        color=PALETTE["T1"],
        label="G0 connectivity recoveries",
    )
    axes[0].set_yscale("log")
    axes[0].set_xticks(gx, [f"G{value}" for value in gx])
    axes[0].set_ylabel("Accepted events (log scale)")
    axes[0].set_title("Cycles and known-space reconnections")
    axes[0].legend(frameon=False, loc="lower right")
    _clean(axes[0])
    _panel(axes[0], "A")
    reverse_fraction = (
        generated["immediate_reverse_cycle_events"]
        / generated["derivation_events"].replace(0, np.nan)
    )
    axes[1].plot(
        gx,
        reverse_fraction,
        marker="o",
        color=PALETTE["red"],
        label="reverse-cycle fraction",
    )
    axes[1].plot(
        gx,
        generated["known_connectivity_reconnection_fraction"],
        marker="s",
        color=PALETTE["T1"],
        label="G0 reconnection fraction",
    )
    axes[1].set_xticks(gx, [f"G{value}" for value in gx])
    axes[1].set_ylim(bottom=0)
    axes[1].set_ylabel("Fraction of accepted events")
    axes[1].set_title("Generation-normalized reconnection rates")
    axes[1].legend(frameon=False, loc="upper right")
    _clean(axes[1])
    _panel(axes[1], "B")
    recorder.add(
        "FigS5",
        "A",
        generated,
        source_path=generation_path,
        description="Immediate reverse-cycle and known-space recovery events.",
    )
    recorder.add(
        "FigS5",
        "B",
        pd.DataFrame(
            {
                "generation": gx,
                "immediate_reverse_cycle_fraction": reverse_fraction,
                "known_connectivity_reconnection_fraction": generated[
                    "known_connectivity_reconnection_fraction"
                ].to_numpy(),
                "interpretation_layer": generated[
                    "interpretation_layer"
                ].to_numpy(),
            }
        ),
        source_path=generation_path,
        description="Accepted-event-normalized cycle and reconnection rates.",
    )
    figure.tight_layout(pad=1.2, w_pad=1.8)
    outputs["FigS5"] = _save(
        figure, output_dir, "FigS5_cycles_and_known_space_reconnections"
    )
    plt.close(figure)

    transition_path = (
        analysis_dir / "functional_state_transition_summary.tsv"
    )
    transitions = read_table(transition_path)
    transitions["generation"] = pd.to_numeric(
        transitions["generation"], errors="coerce"
    ).astype(int)
    transitions["derivation_event_count"] = pd.to_numeric(
        transitions["derivation_event_count"], errors="coerce"
    )
    top_transition_ids = (
        transitions.groupby("functional_state_transition")[
            "derivation_event_count"
        ]
        .sum()
        .nlargest(40)
        .index
    )
    transition_panel = transitions[
        transitions["functional_state_transition"].isin(top_transition_ids)
    ].copy()
    transition_pivot = transition_panel.pivot_table(
        index="functional_state_transition",
        columns="generation",
        values="derivation_event_count",
        fill_value=0,
    )
    transition_pivot = transition_pivot.loc[
        transition_pivot.sum(axis=1).sort_values().index
    ]
    figure, axis = plt.subplots(figsize=(7.6, 9.0))
    image = axis.imshow(
        np.log10(transition_pivot.to_numpy() + 1),
        aspect="auto",
        cmap="magma",
    )
    axis.set_yticks(
        range(len(transition_pivot)),
        [_wrap(value, 42) for value in transition_pivot.index],
        fontsize=5.3,
    )
    axis.set_xticks(
        range(len(transition_pivot.columns)),
        [f"G{int(value)}" for value in transition_pivot.columns],
    )
    axis.set_title("Extended functional-state transition spectrum")
    figure.colorbar(
        image,
        ax=axis,
        fraction=0.03,
        pad=0.03,
        label="log10(events + 1)",
    )
    recorder.add(
        "FigS6",
        "A",
        transition_panel,
        source_path=transition_path,
        description="Forty most frequent counted functional-state transitions.",
    )
    figure.tight_layout(pad=1.2)
    outputs["FigS6"] = _save(
        figure, output_dir, "FigS6_extended_functional_state_transitions"
    )
    plt.close(figure)

    edit_path = analysis_dir / "reaction_edit_landscape.tsv"
    edits = read_table(edit_path)
    for column in (
        "generation",
        "observed_changed_source_atoms",
        "derivation_event_count",
        "median_source_atom_retention",
        "median_source_product_tanimoto",
    ):
        edits[column] = pd.to_numeric(edits[column], errors="coerce")
    edits = edits.dropna(
        subset=[
            "median_source_atom_retention",
            "median_source_product_tanimoto",
        ]
    )
    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))
    for generation, group in edits.groupby("generation", sort=True):
        axes[0].scatter(
            group["median_source_atom_retention"],
            group["median_source_product_tanimoto"],
            s=np.clip(
                np.sqrt(group["derivation_event_count"]) * 1.8, 5, 90
            ),
            color=PALETTE[f"G{int(generation)}"],
            alpha=0.55,
            linewidths=0,
            label=f"G{int(generation)}",
        )
    axes[0].set_xlim(0, 1.02)
    axes[0].set_ylim(0, 1.02)
    axes[0].set_xlabel("Median source-atom retention")
    axes[0].set_ylabel("Median parent-product Morgan Tanimoto")
    axes[0].set_title("Event-weighted reaction-edit classes")
    axes[0].legend(frameon=False)
    _clean(axes[0])
    _panel(axes[0], "A")
    locality = (
        edits.groupby(["generation", "observed_changed_source_atoms"])[
            "derivation_event_count"
        ]
        .sum()
        .reset_index()
    )
    locality = locality[locality["observed_changed_source_atoms"] <= 20]
    locality_pivot = locality.pivot_table(
        index="observed_changed_source_atoms",
        columns="generation",
        values="derivation_event_count",
        fill_value=0,
    )
    image = axes[1].imshow(
        np.log10(locality_pivot.to_numpy() + 1),
        aspect="auto",
        cmap="viridis",
    )
    axes[1].set_yticks(
        range(len(locality_pivot)),
        [str(int(value)) for value in locality_pivot.index],
    )
    axes[1].set_xticks(
        range(len(locality_pivot.columns)),
        [f"G{int(value)}" for value in locality_pivot.columns],
    )
    axes[1].set_ylabel("Observed changed source atoms")
    axes[1].set_title("Locality spectrum of accepted edits")
    figure.colorbar(
        image,
        ax=axes[1],
        fraction=0.046,
        pad=0.04,
        label="log10(events + 1)",
    )
    _panel(axes[1], "B")
    recorder.add(
        "FigS7",
        "A",
        edits,
        source_path=edit_path,
        description=(
            "Reaction-edit classes positioned by median atom retention and "
            "parent-product similarity; marker area reflects event count."
        ),
    )
    recorder.add(
        "FigS7",
        "B",
        locality,
        source_path=edit_path,
        description="Accepted events by generation and observed edit locality.",
    )
    figure.tight_layout(pad=1.2, w_pad=1.8)
    outputs["FigS7"] = _save(
        figure, output_dir, "FigS7_reaction_edit_similarity_and_locality"
    )
    plt.close(figure)

    scope_directories = {
        "T1 + domain": analysis_dir,
        "T2": t2_analysis_dir,
        "T3": t3_analysis_dir,
    }
    descriptor_frames = []
    similarity_frames = []
    for scope, directory in scope_directories.items():
        descriptor_frame = read_table(
            directory / "physicochemical_descriptor_summary.tsv"
        )
        descriptor_frame["evidence_scope"] = scope
        descriptor_frames.append(descriptor_frame)
        similarity_frame = read_table(
            directory / "nearest_G0_similarity.tsv"
        )
        similarity_frame["evidence_scope"] = scope
        similarity_frames.append(similarity_frame)
    sensitivity_descriptors = pd.concat(
        descriptor_frames, ignore_index=True
    )
    sensitivity_similarity = pd.concat(
        similarity_frames, ignore_index=True
    )
    sensitivity_descriptors["generation"] = pd.to_numeric(
        sensitivity_descriptors["generation"], errors="coerce"
    ).astype(int)
    sensitivity_descriptors["mean"] = pd.to_numeric(
        sensitivity_descriptors["mean"], errors="coerce"
    )
    sensitivity_descriptors["standard_deviation"] = pd.to_numeric(
        sensitivity_descriptors["standard_deviation"], errors="coerce"
    )
    sensitivity_similarity["generation"] = pd.to_numeric(
        sensitivity_similarity["generation"], errors="coerce"
    ).astype(int)
    sensitivity_similarity["nearest_G0_tanimoto"] = pd.to_numeric(
        sensitivity_similarity["nearest_G0_tanimoto"], errors="coerce"
    )
    sensitivity_descriptors = sensitivity_descriptors[
        sensitivity_descriptors["generation"].isin([0, 1])
    ]
    sensitivity_similarity = sensitivity_similarity[
        sensitivity_similarity["generation"] == 1
    ]
    descriptor_names = [
        "exact_mass",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "fraction_csp3",
    ]
    heat_rows = []
    for scope in scope_directories:
        group = sensitivity_descriptors[
            sensitivity_descriptors["evidence_scope"] == scope
        ]
        g0 = group[group["generation"] == 0].set_index("descriptor")
        g1 = group[group["generation"] == 1].set_index("descriptor")
        for descriptor in descriptor_names:
            heat_rows.append(
                {
                    "evidence_scope": scope,
                    "descriptor": descriptor,
                    "G1_mean_shift_in_G0_SD": (
                        float(g1.loc[descriptor, "mean"])
                        - float(g0.loc[descriptor, "mean"])
                    )
                    / max(
                        float(g0.loc[descriptor, "standard_deviation"]),
                        1e-12,
                    ),
                }
            )
    descriptor_shift = pd.DataFrame(heat_rows)
    descriptor_heat = descriptor_shift.pivot(
        index="evidence_scope",
        columns="descriptor",
        values="G1_mean_shift_in_G0_SD",
    ).loc[list(scope_directories), descriptor_names]
    figure, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    limit = max(
        1.0, float(np.nanmax(np.abs(descriptor_heat.to_numpy())))
    )
    image = axes[0].imshow(
        descriptor_heat.to_numpy(),
        aspect="auto",
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
    )
    axes[0].set_yticks(
        range(len(descriptor_heat)), descriptor_heat.index
    )
    axes[0].set_xticks(
        range(len(descriptor_names)),
        [_descriptor_label(value) for value in descriptor_names],
        rotation=40,
        ha="right",
    )
    axes[0].set_title("G1 physicochemical displacement")
    figure.colorbar(
        image,
        ax=axes[0],
        fraction=0.046,
        pad=0.04,
        label="Mean shift (G0 SD)",
    )
    _panel(axes[0], "A")
    similarity_values = [
        sensitivity_similarity.loc[
            sensitivity_similarity["evidence_scope"] == scope,
            "nearest_G0_tanimoto",
        ]
        .dropna()
        .to_numpy()
        for scope in scope_directories
    ]
    violin = axes[1].violinplot(
        similarity_values,
        positions=np.arange(len(scope_directories)),
        showmedians=True,
        widths=0.75,
    )
    for body, color in zip(
        violin["bodies"],
        [PALETTE["T1"], PALETTE["T2"], PALETTE["T3"]],
    ):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.72)
    axes[1].set_xticks(
        np.arange(len(scope_directories)), list(scope_directories)
    )
    axes[1].set_ylim(0, 1.02)
    axes[1].set_ylabel("Nearest-G0 Morgan Tanimoto")
    axes[1].set_title("G1 proximity to the known taxane space")
    _clean(axes[1])
    _panel(axes[1], "B")
    recorder.add(
        "FigS8",
        "A",
        descriptor_shift,
        source_path=";".join(
            str(path / "physicochemical_descriptor_summary.tsv")
            for path in scope_directories.values()
        ),
        description=(
            "G1 descriptor mean shifts relative to each evidence scope's "
            "identical G0 population."
        ),
    )
    recorder.add(
        "FigS8",
        "B",
        sensitivity_similarity,
        source_path=";".join(
            str(path / "nearest_G0_similarity.tsv")
            for path in scope_directories.values()
        ),
        description=(
            "Deterministically sampled nearest-G0 similarities for "
            "exclusive T1, T2, and T3 G1 spaces."
        ),
    )
    figure.tight_layout(pad=1.2, w_pad=2.0)
    outputs["FigS8"] = _save(
        figure, output_dir, "FigS8_evidence_layer_G1_sensitivity"
    )
    plt.close(figure)
    return outputs


def render_study_figures(
    analysis_dir: Path,
    provenance_dir: Path,
    sensitivity_dir: Path,
    external_benchmark_dir: Path,
    domain_benchmark_dir: Path,
    t2_analysis_dir: Path,
    t3_analysis_dir: Path,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir = ensure_dir(output_dir)
    recorder = SourceRecorder(output_dir)
    figure_paths = {
        "Fig1": _fig1(provenance_dir, analysis_dir, output_dir, recorder),
        "Fig2": _fig2(
            sensitivity_dir,
            external_benchmark_dir,
            domain_benchmark_dir,
            output_dir,
            recorder,
        ),
        "Fig3": _fig3(analysis_dir, output_dir, recorder),
        "Fig4": _fig4(analysis_dir, output_dir, recorder),
        "Fig5": _fig5(analysis_dir, output_dir, recorder),
    }
    supplementary_paths = _supplement_figures(
        analysis_dir,
        t2_analysis_dir,
        t3_analysis_dir,
        output_dir,
        recorder,
    )
    manifest_path = recorder.write(output_dir)
    summary_path = output_dir / "publication_figure_build_summary.json"
    summary = {
        "figures": {
            key: [str(path) for path in paths]
            for key, paths in figure_paths.items()
        },
        "supplementary_figures": {
            key: [str(path) for path in paths]
            for key, paths in supplementary_paths.items()
        },
        "source_data_manifest": str(manifest_path),
        "formats": ["PDF", "SVG", "PNG_400dpi"],
        "svg_text_remains_editable": True,
        "pdf_font_embedding": "TrueType_Type42",
    }
    write_json(summary, summary_path)
    return {
        "figures": figure_paths,
        "supplementary_figures": supplementary_paths,
        "source_manifest": manifest_path,
        "summary": summary_path,
    }
