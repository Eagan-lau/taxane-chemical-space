from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .io import ensure_dir, read_table, write_table


COLORS = {
    "G0": "#30343B",
    "G1": "#2A9D8F",
    "G2": "#E76F51",
    "G3": "#E9C46A",
    "blue": "#457B9D",
    "gray": "#8A8F98",
    "light_gray": "#D9DCE1",
    "red": "#C44536",
    "green": "#3A7D44",
}


def _setup_matplotlib():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.transparent": False,
        }
    )
    return plt


def _panel_label(axis, label: str) -> None:
    axis.text(
        -0.12,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
    )


def _clean_axis(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="#E6E8EB", linewidth=0.6, zorder=0)
    axis.set_axisbelow(True)


def _save_figure(figure, output_stem: Path) -> list[Path]:
    paths = []
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 400}),
    ):
        path = output_stem.with_suffix(suffix)
        figure.savefig(path, **kwargs)
        paths.append(path)
    return paths


def _read_json(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _generation_color(generation: int) -> str:
    return COLORS.get(f"G{generation}", COLORS["gray"])


def _figure_rule_grammar(
    grammar_summary: dict[str, Any],
    screen_summary: dict[str, Any],
    selection_summary: dict[str, Any],
    selected_grammar: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    plt = _setup_matplotlib()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
    axes = axes.ravel()

    stages = [
        ("T1 evidence\nrules", grammar_summary.get("rows_seen", 0)),
        ("QC-eligible\nrules", grammar_summary.get("eligible_rules", 0)),
        ("Semantic\ngroups", grammar_summary.get("semantic_groups", 0)),
        ("Executable\nrepresentatives", grammar_summary.get("selected_grammar_rules", 0)),
        ("G0-activated\nrules", screen_summary.get("activated_rules", 0)),
        ("Primary taxane\ngrammar", selection_summary.get("primary_rules", 0)),
    ]
    labels = [stage[0] for stage in stages]
    values = np.asarray([float(stage[1] or 0) for stage in stages])
    positions = np.arange(len(stages))
    axes[0].bar(
        positions,
        values,
        color=[
            COLORS["gray"],
            COLORS["blue"],
            COLORS["G1"],
            COLORS["G2"],
            COLORS["G3"],
            COLORS["green"],
        ],
        width=0.72,
    )
    axes[0].set_yscale("symlog", linthresh=10)
    axes[0].set_xticks(positions, labels, rotation=25, ha="right")
    axes[0].set_ylabel("Rule count (symlog)")
    axes[0].set_title("Evidence-to-grammar attrition")
    for position, value in zip(positions, values):
        axes[0].text(
            position,
            max(value, 1) * 1.15,
            f"{int(value):,}",
            ha="center",
            va="bottom",
            fontsize=6.5,
        )
    _clean_axis(axes[0])
    _panel_label(axes[0], "A")

    if not selected_grammar.empty and "template_sources" in selected_grammar:
        source_counts: dict[str, int] = {}
        for value in selected_grammar["template_sources"].fillna(""):
            for source in str(value).split(";"):
                source = source.strip()
                if source:
                    source_counts[source] = source_counts.get(source, 0) + 1
        source_frame = (
            pd.DataFrame(
                [{"source": key, "rule_count": value} for key, value in source_counts.items()]
            )
            .sort_values("rule_count", ascending=True)
            .tail(10)
        )
        axes[1].barh(
            source_frame["source"],
            source_frame["rule_count"],
            color=COLORS["blue"],
        )
        axes[1].set_xlabel("Selected rule records")
    axes[1].set_title("Database support of the executable grammar")
    _clean_axis(axes[1])
    _panel_label(axes[1], "B")

    if not selected_grammar.empty and "reaction_type" in selected_grammar:
        reaction_counts = (
            selected_grammar["reaction_type"]
            .replace("", "unassigned")
            .value_counts()
            .head(12)
            .sort_values()
        )
        axes[2].barh(
            reaction_counts.index,
            reaction_counts.values,
            color=COLORS["G1"],
        )
        axes[2].set_xlabel("Rule count")
    axes[2].set_title("Transformation ontology composition")
    _clean_axis(axes[2])
    _panel_label(axes[2], "C")

    if not selected_grammar.empty:
        delta_column = (
            "structural_element_delta"
            if "structural_element_delta" in selected_grammar
            else "reaction_delta_fingerprint"
        )
        delta_counts = (
            selected_grammar[delta_column]
            .replace("", "bond/stereo edit without heavy-element change")
            .value_counts()
            .head(12)
            .sort_values()
        )
        axes[3].barh(
            delta_counts.index,
            delta_counts.values,
            color=COLORS["G2"],
        )
        axes[3].set_xlabel("Rule count")
    axes[3].set_title("Most frequent structural element deltas")
    _clean_axis(axes[3])
    _panel_label(axes[3], "D")
    figure.tight_layout(pad=1.2)
    return _save_figure(figure, output_dir / "Fig1_reaction_grammar_construction")


def _figure_benchmark(
    benchmark_summary: dict[str, Any],
    positives: pd.DataFrame,
    groups: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    plt = _setup_matplotlib()
    figure, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))

    positive_rate = benchmark_summary.get("positive_connectivity_recovery_rate", 0) or 0
    decoy_rate = benchmark_summary.get("decoy_connectivity_match_rate", 0) or 0
    rates = [positive_rate, decoy_rate]
    intervals = [
        benchmark_summary.get("positive_recovery_wilson_95ci", [positive_rate, positive_rate]),
        benchmark_summary.get("decoy_match_wilson_95ci", [decoy_rate, decoy_rate]),
    ]
    errors = np.asarray(
        [
            [rate - interval[0] for rate, interval in zip(rates, intervals)],
            [interval[1] - rate for rate, interval in zip(rates, intervals)],
        ]
    )
    axes[0].bar(
        [0, 1],
        rates,
        color=[COLORS["G1"], COLORS["gray"]],
        yerr=errors,
        capsize=3,
        width=0.65,
    )
    axes[0].set_xticks([0, 1], ["Known pathway", "Matched decoys"])
    axes[0].set_ylabel("Connectivity recovery rate")
    axes[0].set_ylim(0, min(1, max(rates + [0.2]) * 1.5))
    axes[0].set_title("Leakage-controlled recovery")
    _clean_axis(axes[0])
    _panel_label(axes[0], "A")

    if not positives.empty and "recovered_connectivity" in positives:
        ordered = positives.sort_values(["benchmark_group_id", "benchmark_id"]).reset_index(drop=True)
        values = ordered["recovered_connectivity"].astype(bool).astype(int).to_numpy()[None, :]
        axes[1].imshow(
            values,
            aspect="auto",
            cmap="Greys",
            vmin=0,
            vmax=1,
            interpolation="nearest",
        )
        axes[1].set_yticks([])
        axes[1].set_xticks(
            np.arange(len(ordered)),
            ordered.get("step_label", ordered["benchmark_id"]),
            rotation=90,
            fontsize=5.5,
        )
        axes[1].set_xlabel("Curated pathway reactions")
    axes[1].set_title("Reaction-level recovery")
    _panel_label(axes[1], "B")

    if not groups.empty:
        display = groups.sort_values(
            ["reaction_count", "connectivity_recovered"], ascending=False
        ).head(15)
        x = np.arange(len(display))
        axes[2].bar(
            x,
            display["reaction_count"],
            color=COLORS["light_gray"],
            label="tested",
        )
        axes[2].bar(
            x,
            display["connectivity_recovered"],
            color=COLORS["blue"],
            label="recovered",
        )
        axes[2].set_xticks(x, [f"RCG{i+1}" for i in range(len(display))], rotation=90)
        axes[2].set_ylabel("Reaction count")
        axes[2].legend(frameon=False)
    axes[2].set_title("Reaction-center-group performance")
    _clean_axis(axes[2])
    _panel_label(axes[2], "C")
    figure.tight_layout(pad=1.0)
    return _save_figure(figure, output_dir / "Fig2_leakage_controlled_benchmark")


def _figure_generation(
    generation: pd.DataFrame,
    rejections: pd.DataFrame,
    convergence: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    plt = _setup_matplotlib()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    axes = axes.ravel()
    generation = generation.sort_values("generation")
    x = generation["generation"].astype(int).to_numpy()
    colors = [_generation_color(value) for value in x]

    axes[0].bar(
        x,
        generation["unique_nodes_first_observed"],
        color=colors,
        width=0.7,
    )
    axes[0].set_xticks(x, [f"G{value}" for value in x])
    axes[0].set_ylabel("Unique structures")
    axes[0].set_title("Generation-wise chemical-space expansion")
    _clean_axis(axes[0])
    _panel_label(axes[0], "A")

    generation_nonzero = generation[generation["generation"] > 0]
    gx = generation_nonzero["generation"].astype(int).to_numpy()
    axes[1].bar(
        gx - 0.22,
        generation_nonzero["raw_product_tuples"],
        width=0.22,
        color=COLORS["gray"],
        label="raw products",
    )
    axes[1].bar(
        gx,
        generation_nonzero["rejected_product_events"],
        width=0.22,
        color=COLORS["G2"],
        label="rejected",
    )
    axes[1].bar(
        gx + 0.22,
        generation_nonzero["derivation_events"],
        width=0.22,
        color=COLORS["G1"],
        label="accepted events",
    )
    axes[1].set_yscale("symlog", linthresh=10)
    axes[1].set_xticks(gx, [f"G{value}" for value in gx])
    axes[1].set_ylabel("Event count (symlog)")
    axes[1].legend(frameon=False)
    axes[1].set_title("Per-generation QC attrition")
    _clean_axis(axes[1])
    _panel_label(axes[1], "B")

    axes[2].bar(
        x,
        generation["cumulative_unique_nodes"],
        color=colors,
        width=0.7,
    )
    axes[2].set_xticks(x, [f"G{value}" for value in x])
    axes[2].set_ylabel("Cumulative unique structures")
    axes[2].set_title("Cumulative grammar-accessible space")
    _clean_axis(axes[2])
    _panel_label(axes[2], "C")

    generated = convergence[convergence["generation"] > 0].copy()
    if not generated.empty:
        parent_counts = pd.to_numeric(
            generated["unique_parent_count"], errors="coerce"
        ).fillna(0)
        bins = np.arange(0.5, min(parent_counts.max(), 20) + 1.5)
        axes[3].hist(
            np.clip(parent_counts, 1, 20),
            bins=bins,
            color=COLORS["blue"],
            edgecolor="white",
        )
        axes[3].set_xlabel("Unique immediate parents (values >20 clipped)")
        axes[3].set_ylabel("Generated structures")
        axes[3].set_yscale("log")
    axes[3].set_title("Convergent derivation multiplicity")
    _clean_axis(axes[3])
    _panel_label(axes[3], "D")
    figure.tight_layout(pad=1.2)
    return _save_figure(figure, output_dir / "Fig3_generation_expansion_and_QC")


def _figure_chemical_space(
    projection: pd.DataFrame,
    nearest: pd.DataFrame,
    descriptor_summary: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    plt = _setup_matplotlib()
    figure, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))

    for generation, group in projection.groupby("generation", sort=True):
        axes[0].scatter(
            group["molecular_fp_axis_1"],
            group["molecular_fp_axis_2"],
            s=4 if int(generation) > 0 else 8,
            alpha=0.35 if int(generation) > 0 else 0.8,
            linewidths=0,
            color=_generation_color(int(generation)),
            label=f"G{int(generation)}",
            rasterized=True,
        )
    axes[0].set_xlabel("Fingerprint axis 1")
    axes[0].set_ylabel("Fingerprint axis 2")
    axes[0].legend(frameon=False, markerscale=2)
    axes[0].set_title("Reaction-grammar chemical space")
    _clean_axis(axes[0])
    _panel_label(axes[0], "A")

    generations = sorted(int(value) for value in nearest["generation"].unique())
    values = [
        nearest.loc[
            nearest["generation"].astype(int) == generation,
            "nearest_G0_tanimoto",
        ].astype(float)
        for generation in generations
    ]
    if values:
        violin = axes[1].violinplot(
            values,
            positions=generations,
            showmedians=True,
            showextrema=False,
        )
        for body, generation in zip(violin["bodies"], generations):
            body.set_facecolor(_generation_color(generation))
            body.set_edgecolor("none")
            body.set_alpha(0.75)
        violin["cmedians"].set_color("#222222")
        axes[1].set_xticks(generations, [f"G{value}" for value in generations])
    axes[1].set_ylim(0, 1.02)
    axes[1].set_ylabel("Nearest-G0 Morgan Tanimoto")
    axes[1].set_title("Displacement from known taxanes")
    _clean_axis(axes[1])
    _panel_label(axes[1], "B")

    pivot = descriptor_summary.pivot(
        index="descriptor", columns="generation", values="median"
    )
    baseline = pivot.get(0)
    if baseline is not None:
        scale = pivot.std(axis=1).replace(0, 1)
        standardized = pivot.sub(baseline, axis=0).div(scale, axis=0)
        image = axes[2].imshow(
            standardized,
            aspect="auto",
            cmap="coolwarm",
            vmin=-2,
            vmax=2,
            interpolation="nearest",
        )
        axes[2].set_yticks(np.arange(len(standardized)), standardized.index)
        axes[2].set_xticks(
            np.arange(len(standardized.columns)),
            [f"G{int(value)}" for value in standardized.columns],
        )
        figure.colorbar(image, ax=axes[2], fraction=0.046, pad=0.03, label="Median shift from G0")
    axes[2].set_title("Physicochemical displacement")
    _panel_label(axes[2], "C")
    figure.tight_layout(pad=1.0)
    return _save_figure(figure, output_dir / "Fig4_chemical_space_displacement")


def _figure_grammar_and_routes(
    grammar_usage: pd.DataFrame,
    transitions: pd.DataFrame,
    edit_landscape: pd.DataFrame,
    bridges: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    plt = _setup_matplotlib()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.0))
    axes = axes.ravel()

    usage = grammar_usage.head(20).sort_values("derivation_event_count")
    axes[0].barh(
        [f"R{i+1}" for i in range(len(usage))],
        usage["derivation_event_count"],
        color=COLORS["blue"],
    )
    axes[0].set_xlabel("Derivation events")
    axes[0].set_title("Most frequently used grammar productions")
    _clean_axis(axes[0])
    _panel_label(axes[0], "A")

    if not grammar_usage.empty:
        fractions = (
            grammar_usage["derivation_event_count"].astype(float)
            / grammar_usage["derivation_event_count"].astype(float).sum()
        )
        axes[1].plot(
            np.arange(1, len(fractions) + 1),
            fractions.cumsum(),
            color=COLORS["G2"],
            linewidth=1.5,
        )
        axes[1].axhline(0.8, color=COLORS["gray"], linestyle="--", linewidth=0.8)
        axes[1].set_xscale("log")
        axes[1].set_ylim(0, 1.02)
        axes[1].set_xlabel("Ranked grammar productions (log)")
        axes[1].set_ylabel("Cumulative event fraction")
    axes[1].set_title("Grammar-use concentration")
    _clean_axis(axes[1])
    _panel_label(axes[1], "B")

    top_transitions = (
        transitions.groupby("functional_state_transition")["derivation_event_count"]
        .sum()
        .sort_values()
        .tail(15)
    )
    axes[2].barh(
        top_transitions.index,
        top_transitions.values,
        color=COLORS["G1"],
    )
    axes[2].set_xlabel("Derivation events")
    axes[2].set_title("Functional-state transitions")
    _clean_axis(axes[2])
    _panel_label(axes[2], "C")

    if not bridges.empty:
        candidate = bridges[bridges["latent_bridge_candidate"].astype(bool)].copy()
        candidate["bridge_score"] = (
            candidate["known_G0_ancestor_count"].astype(int)
            * candidate["known_G0_descendant_count"].astype(int).clip(lower=1)
        )
        candidate = candidate.sort_values("bridge_score", ascending=False).head(30)
        axes[3].scatter(
            candidate["known_G0_ancestor_count"],
            candidate["known_G0_descendant_count"],
            s=20 + 8 * candidate["bridge_score"],
            color=COLORS["G2"],
            alpha=0.7,
            edgecolor="white",
            linewidth=0.4,
        )
        axes[3].set_xlabel("Distinct known G0 ancestors")
        axes[3].set_ylabel("Distinct known G0 descendants")
    axes[3].set_title("Latent bridge candidates")
    _clean_axis(axes[3])
    _panel_label(axes[3], "D")
    figure.tight_layout(pad=1.2)
    return _save_figure(figure, output_dir / "Fig5_grammar_usage_and_route_convergence")


def render_publication_figures(
    analysis_dir: Path,
    output_dir: Path,
    *,
    grammar_summary_path: Path | None = None,
    screen_summary_path: Path | None = None,
    selection_summary_path: Path | None = None,
    selected_grammar_path: Path | None = None,
    benchmark_dir: Path | None = None,
) -> dict[str, list[Path]]:
    output_dir = ensure_dir(output_dir)
    grammar_summary = _read_json(grammar_summary_path)
    screen_summary = _read_json(screen_summary_path)
    selection_summary = _read_json(selection_summary_path)
    selected_grammar = (
        read_table(selected_grammar_path)
        if selected_grammar_path and selected_grammar_path.exists()
        else pd.DataFrame()
    )
    generation = read_table(analysis_dir / "generation_expansion_summary.tsv")
    descriptor_summary = read_table(
        analysis_dir / "physicochemical_descriptor_summary.tsv"
    )
    nearest = read_table(analysis_dir / "nearest_G0_similarity.tsv")
    projection = read_table(
        analysis_dir / "chemical_space_fingerprint_projection.tsv"
    )
    grammar_usage = read_table(analysis_dir / "reaction_grammar_usage.tsv")
    transitions = read_table(
        analysis_dir / "functional_state_transition_summary.tsv"
    )
    convergence = read_table(
        analysis_dir / "convergence_and_route_multiplicity.tsv"
    )
    bridges = read_table(analysis_dir / "latent_bridge_candidates.tsv")
    edit_landscape = read_table(analysis_dir / "reaction_edit_landscape.tsv")
    rejections = read_table(analysis_dir / "generation_rejection_summary.tsv")

    figure_paths = {
        "Fig1": _figure_rule_grammar(
            grammar_summary,
            screen_summary,
            selection_summary,
            selected_grammar,
            output_dir,
        ),
        "Fig3": _figure_generation(
            generation, rejections, convergence, output_dir
        ),
        "Fig4": _figure_chemical_space(
            projection, nearest, descriptor_summary, output_dir
        ),
        "Fig5": _figure_grammar_and_routes(
            grammar_usage, transitions, edit_landscape, bridges, output_dir
        ),
    }
    if benchmark_dir and (benchmark_dir / "benchmark_summary.json").exists():
        figure_paths["Fig2"] = _figure_benchmark(
            _read_json(benchmark_dir / "benchmark_summary.json"),
            read_table(benchmark_dir / "benchmark_positive_recovery.tsv"),
            read_table(benchmark_dir / "benchmark_reaction_center_groups.tsv"),
            output_dir,
        )
    chart_map_rows = []
    for figure_name, paths in figure_paths.items():
        chart_map_rows.append(
            {
                "figure": figure_name,
                "scientific_question": {
                    "Fig1": "How was the evidence library compressed into an executable grammar?",
                    "Fig2": "Does the grammar recover held-out pathway chemistry above matched decoys?",
                    "Fig3": "How does the grammar-accessible space expand and pass QC across G1-G3?",
                    "Fig4": "How far do generated structures move from known taxane chemical space?",
                    "Fig5": "Which grammar productions dominate and where do routes converge?",
                }.get(figure_name, ""),
                "source_tables": {
                    "Fig1": "grammar summaries; selected grammar",
                    "Fig2": "benchmark_positive_recovery.tsv; benchmark_reaction_center_groups.tsv",
                    "Fig3": "generation_expansion_summary.tsv; convergence_and_route_multiplicity.tsv",
                    "Fig4": "chemical_space_fingerprint_projection.tsv; nearest_G0_similarity.tsv; physicochemical_descriptor_summary.tsv",
                    "Fig5": "reaction_grammar_usage.tsv; functional_state_transition_summary.tsv; latent_bridge_candidates.tsv",
                }.get(figure_name, ""),
                "pdf": str(next(path for path in paths if path.suffix == ".pdf")),
                "svg": str(next(path for path in paths if path.suffix == ".svg")),
                "png": str(next(path for path in paths if path.suffix == ".png")),
            }
        )
    write_table(pd.DataFrame(chart_map_rows), output_dir / "figure_chart_map.tsv")
    return figure_paths
