from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .io import ensure_dir, read_table, write_json, write_table


PLACEHOLDER_PREFIX = "{{FINAL_"

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


def _percent(value: float, digits: int = 1) -> str:
    return f"{100.0 * float(value):.{digits}f}"


def _generation_row(frame: pd.DataFrame, generation: int) -> pd.Series:
    matches = frame[pd.to_numeric(frame["generation"]) == generation]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one G{generation} row, found {len(matches)}"
        )
    return matches.iloc[0]


def _concentration_row(frame: pd.DataFrame, label: str) -> pd.Series:
    matches = frame[frame["generation_scope"].astype(str) == label]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one concentration row for {label}, found {len(matches)}"
        )
    return matches.iloc[0]


def _trend_sentence(
    values: dict[int, float],
    *,
    metric: str,
    formatter,
) -> str:
    ordered = sorted(values)
    rendered = ", ".join(
        f"G{generation} {formatter(values[generation])}"
        for generation in ordered
    )
    return f"{metric} was {rendered}"


def _three_generation_decline_sentence(
    values: dict[int, float],
    *,
    metric: str,
) -> str:
    return (
        f"{metric} declined from {values[1]:.3f} at G1 to "
        f"{values[2]:.3f} at G2 and {values[3]:.3f} at G3"
    )


def _descriptor_sentence(descriptors: pd.DataFrame) -> str:
    descriptor_frame = descriptors.copy()
    for column in ("generation", "mean", "standard_deviation"):
        descriptor_frame[column] = pd.to_numeric(
            descriptor_frame[column], errors="coerce"
        )
    g0 = descriptor_frame[
        descriptor_frame["generation"] == 0
    ].set_index("descriptor")
    g3 = descriptor_frame[
        descriptor_frame["generation"] == 3
    ].set_index("descriptor")
    shared = g0.index.intersection(g3.index)
    shifts = []
    for descriptor in shared:
        standard_deviation = float(g0.loc[descriptor, "standard_deviation"])
        if not np.isfinite(standard_deviation) or standard_deviation == 0:
            continue
        shift = (
            float(g3.loc[descriptor, "mean"])
            - float(g0.loc[descriptor, "mean"])
        ) / standard_deviation
        shifts.append((abs(shift), descriptor, shift))
    shifts.sort(reverse=True)
    selected = shifts[:3]
    if not selected:
        return "no descriptor with a finite G0-standardized shift was available"
    pieces = [
        (
            f"{DESCRIPTOR_LABELS.get(descriptor, descriptor.replace('_', ' '))} "
            f"{'increased' if shift > 0 else 'decreased'} by "
            f"{abs(shift):.2f} G0 standard deviations"
        )
        for _magnitude, descriptor, shift in selected
    ]
    return "; ".join(pieces)


def _reaction_locality_sentence(edit_frame: pd.DataFrame) -> str:
    edit = edit_frame.copy()
    edit["observed_changed_source_atoms"] = pd.to_numeric(
        edit["observed_changed_source_atoms"], errors="coerce"
    )
    edit["derivation_event_count"] = pd.to_numeric(
        edit["derivation_event_count"], errors="coerce"
    ).fillna(0)
    total = float(edit["derivation_event_count"].sum())
    if total == 0:
        return "no accepted edit events were available"
    local = float(
        edit.loc[
            edit["observed_changed_source_atoms"] <= 1,
            "derivation_event_count",
        ].sum()
    )
    broad = float(
        edit.loc[
            edit["observed_changed_source_atoms"] > 3,
            "derivation_event_count",
        ].sum()
    )
    return (
        f"{_percent(local / total)}% changed at most one mapped source atom "
        f"and only {_percent(broad / total)}% changed more than three"
    )


def _bridge_context(
    bridges: pd.DataFrame, pair_frame: pd.DataFrame
) -> dict[str, Any]:
    candidate_mask = (
        bridges["latent_bridge_candidate"]
        .astype(str)
        .str.lower()
        .isin({"1", "true", "yes"})
    )
    candidates = bridges[candidate_mask].copy()
    values = pd.to_numeric(
        candidates.get(
            "distinct_G0_pair_bridge_count", pd.Series(dtype=float)
        ),
        errors="coerce",
    ).dropna()
    if len(values):
        breadth_sentence = (
            f"the median candidate supported {float(values.median()):.0f} "
            f"directed G0 pairs, the 90th percentile supported "
            f"{float(values.quantile(0.90)):.0f}, and the maximum was "
            f"{float(values.max()):.0f}"
        )
    else:
        breadth_sentence = "no conservative directed G0 bridge was detected"
    return {
        "latent_bridge_count": int(len(candidates)),
        "bridged_pair_count": int(len(pair_frame)),
        "bridge_breadth_sentence": breadth_sentence,
    }


def build_result_context(analysis_dir: Path) -> tuple[dict[str, Any], pd.DataFrame]:
    generation = read_table(
        analysis_dir / "generation_expansion_summary.tsv"
    )
    concentration = read_table(
        analysis_dir / "reaction_grammar_use_concentration.tsv"
    )
    grammar_usage = read_table(
        analysis_dir / "reaction_grammar_usage.tsv"
    )
    nearest = read_table(analysis_dir / "nearest_G0_similarity.tsv")
    descriptors = read_table(
        analysis_dir / "physicochemical_descriptor_summary.tsv"
    )
    convergence = read_table(
        analysis_dir / "convergence_and_route_multiplicity.tsv"
    )
    bridges = read_table(analysis_dir / "latent_bridge_candidates.tsv")
    pair_path = analysis_dir / "known_G0_pair_bridge_summary.tsv"
    if pair_path.exists() and pair_path.stat().st_size > 1:
        try:
            pairs = read_table(pair_path)
        except pd.errors.EmptyDataError:
            pairs = pd.DataFrame()
    else:
        pairs = pd.DataFrame()
    edits = read_table(analysis_dir / "reaction_edit_landscape.tsv")

    available_generations = set(
        pd.to_numeric(generation["generation"], errors="coerce")
        .dropna()
        .astype(int)
    )
    missing_generations = {0, 1, 2, 3} - available_generations
    if missing_generations:
        raise ValueError(
            "Final manuscript requires G0-G3 analysis; missing "
            f"{sorted(missing_generations)}"
        )

    g3 = _generation_row(generation, 3)
    all_concentration = _concentration_row(
        concentration, "all_generations"
    )
    generation_concentration = {
        generation_value: _concentration_row(
            concentration, f"G{generation_value}"
        )
        for generation_value in (1, 2, 3)
    }

    novelty_values = {
        generation_value: float(
            _generation_row(
                generation, generation_value
            )["unique_node_yield_per_raw_product"]
        )
        for generation_value in (1, 2, 3)
    }
    novelty_sentence = _three_generation_decline_sentence(
        novelty_values,
        metric="The unique-structure yield per raw product",
    )

    nearest["generation"] = pd.to_numeric(
        nearest["generation"], errors="coerce"
    )
    nearest["nearest_G0_tanimoto"] = pd.to_numeric(
        nearest["nearest_G0_tanimoto"], errors="coerce"
    )
    similarity_values = {
        generation_value: float(
            nearest.loc[
                nearest["generation"] == generation_value,
                "nearest_G0_tanimoto",
            ].median()
        )
        for generation_value in (1, 2, 3)
    }
    similarity_sentence = _three_generation_decline_sentence(
        similarity_values,
        metric="Median nearest-G0 Tanimoto similarity",
    )

    convergence["generation"] = pd.to_numeric(
        convergence["generation"], errors="coerce"
    )
    convergence["_is_convergent"] = (
        convergence["is_convergent"]
        .astype(str)
        .str.lower()
        .isin({"1", "true", "yes"})
    )
    convergence_fraction = (
        convergence.groupby("generation")["_is_convergent"].mean().to_dict()
    )
    bridge_context = _bridge_context(bridges, pairs)

    concentration_sentence = (
        "the top rule accounted for "
        + ", ".join(
            f"{_percent(float(row['top1_rule_event_fraction']))}% at G{generation_value}"
            for generation_value, row in generation_concentration.items()
        )
        + ", whereas the effective rule number was "
        + ", ".join(
            f"{float(row['effective_rule_number_exp_shannon']):.2f} at G{generation_value}"
            for generation_value, row in generation_concentration.items()
        )
    )

    context: dict[str, Any] = {
        "FINAL_G3_UNIQUE_STRUCTURES": (
            f"{int(g3['unique_nodes_first_observed']):,}"
        ),
        "FINAL_CUMULATIVE_STRUCTURES": (
            f"{int(g3['cumulative_unique_nodes']):,}"
        ),
        "FINAL_G3_ACCEPTED_EVENTS": f"{int(g3['derivation_events']):,}",
        "FINAL_TOP5_EVENT_FRACTION_PERCENT": _percent(
            float(all_concentration["top5_rule_event_fraction"])
        ),
        "FINAL_G0_G3_ACTIVE_RULE_COUNT": int(
            grammar_usage["grammar_rule_id"].nunique()
        ),
        "FINAL_G0_G3_ACTIVE_SEMANTIC_GROUP_COUNT": int(
            grammar_usage["semantic_group_id"].nunique()
        ),
        "FINAL_LATENT_BRIDGE_COUNT": bridge_context[
            "latent_bridge_count"
        ],
        "FINAL_BRIDGED_G0_PAIR_COUNT": bridge_context[
            "bridged_pair_count"
        ],
        "FINAL_NOVELTY_SENTENCE": novelty_sentence,
        "FINAL_SIMILARITY_SENTENCE": similarity_sentence,
        "FINAL_DESCRIPTOR_SENTENCE": _descriptor_sentence(descriptors),
        "FINAL_GRAMMAR_CONCENTRATION_SENTENCE": concentration_sentence,
        "FINAL_REACTION_LOCALITY_SENTENCE": _reaction_locality_sentence(
            edits
        ),
        "FINAL_G1_CONVERGENCE_PERCENT": _percent(
            float(convergence_fraction.get(1, 0.0))
        ),
        "FINAL_G2_CONVERGENCE_PERCENT": _percent(
            float(convergence_fraction.get(2, 0.0))
        ),
        "FINAL_G3_CONVERGENCE_PERCENT": _percent(
            float(convergence_fraction.get(3, 0.0))
        ),
        "FINAL_BRIDGE_BREADTH_SENTENCE": bridge_context[
            "bridge_breadth_sentence"
        ],
    }
    context_rows = pd.DataFrame(
        [
            {"placeholder": key, "value": value}
            for key, value in sorted(context.items())
        ]
    )
    return context, context_rows


def render_manuscript(
    template_path: Path,
    analysis_dir: Path,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    context, context_rows = build_result_context(analysis_dir)
    manuscript = template_path.read_text(encoding="utf-8")
    for key, value in context.items():
        manuscript = manuscript.replace(f"{{{{{key}}}}}", str(value))
    unresolved = sorted(
        {
            token.split("}}", 1)[0] + "}}"
            for token in manuscript.split("{{")
            if token.startswith("FINAL_")
        }
    )
    if unresolved or PLACEHOLDER_PREFIX in manuscript:
        raise ValueError(
            f"Unresolved final-result placeholders: {unresolved}"
        )

    manuscript_path = output_dir / "MANUSCRIPT_FINAL.md"
    manuscript_path.write_text(manuscript, encoding="utf-8")
    context_path = output_dir / "manuscript_result_context.tsv"
    write_table(context_rows, context_path)
    summary_path = output_dir / "manuscript_build_summary.json"
    summary = {
        "template": str(template_path),
        "analysis_dir": str(analysis_dir),
        "manuscript": str(manuscript_path),
        "result_context": str(context_path),
        "placeholder_count": int(len(context)),
        "unresolved_placeholder_count": 0,
    }
    write_json(summary, summary_path)
    return {
        "manuscript": manuscript_path,
        "context": context_path,
        "summary": summary_path,
    }
