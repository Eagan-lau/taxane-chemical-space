from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

from .io import ensure_dir, read_table, write_json, write_table


def _summary_record(
    label: str,
    summary_path: Path,
    g1_keys: set[str],
) -> dict[str, Any]:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    generation = next(
        (
            item
            for item in summary.get("generations", [])
            if int(item.get("generation", -1)) == 1
        ),
        {},
    )
    raw_products = int(generation.get("raw_product_tuples", 0))
    events = int(generation.get("accepted_derivation_events", 0))
    rejections = sum(
        int(value)
        for value in generation.get("rejection_reason_counts", {}).values()
    )
    return {
        "evidence_scope": label,
        "comparison_scope": "exclusive evidence layer; not cumulative",
        "compiled_rules": int(summary.get("compiled_grammar_rules", 0)),
        "activated_rules_G1": int(generation.get("activated_rules", 0)),
        "G0_seed_structures": int(
            summary.get("G0", {}).get("unique_full_stereo_structures", 0)
        ),
        "unique_G1_structures": int(len(g1_keys)),
        "accepted_G1_derivation_events": events,
        "matched_parent_rule_pairs": int(
            generation.get("matched_parent_rule_pairs", 0)
        ),
        "raw_product_tuples": raw_products,
        "rejected_product_events": rejections,
        "accepted_event_fraction_of_raw_products": (
            events / raw_products if raw_products else 0.0
        ),
        "known_G0_full_recovery_events": int(
            generation.get("known_G0_full_recovery_events", 0)
        ),
        "known_G0_connectivity_only_recovery_events": int(
            generation.get("known_G0_connectivity_only_recovery_events", 0)
        ),
        "grammar_compile_failures": int(
            summary.get("grammar_compile_failures", 0)
        ),
    }


def compare_g1_sensitivity_spaces(
    primary_nodes_path: Path,
    primary_summary_path: Path,
    t2_nodes_path: Path,
    t2_summary_path: Path,
    t3_nodes_path: Path,
    t3_summary_path: Path,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    inputs = {
        "T1_primary_plus_taxane_domain": (
            primary_nodes_path,
            primary_summary_path,
        ),
        "T2_exclusive": (t2_nodes_path, t2_summary_path),
        "T3_exclusive": (t3_nodes_path, t3_summary_path),
    }
    frames: dict[str, pd.DataFrame] = {}
    key_sets: dict[str, set[str]] = {}
    for label, (nodes_path, _summary_path) in inputs.items():
        nodes = read_table(nodes_path).fillna("")
        nodes["generation_first"] = pd.to_numeric(
            nodes["generation_first"], errors="coerce"
        ).fillna(-1)
        g1 = nodes[nodes["generation_first"] == 1].copy()
        frames[label] = g1
        key_sets[label] = set(g1["full_inchikey"].astype(str))

    tier_summary = pd.DataFrame(
        [
            _summary_record(label, summary_path, key_sets[label])
            for label, (_nodes_path, summary_path) in inputs.items()
        ]
    )
    primary_keys = key_sets["T1_primary_plus_taxane_domain"]
    tier_summary["G1_overlap_with_primary"] = tier_summary[
        "evidence_scope"
    ].map(lambda label: len(key_sets[label] & primary_keys))
    tier_summary["G1_unique_relative_to_primary"] = tier_summary[
        "evidence_scope"
    ].map(lambda label: len(key_sets[label] - primary_keys))
    tier_summary["fraction_of_layer_G1_also_in_primary"] = tier_summary[
        "evidence_scope"
    ].map(
        lambda label: (
            len(key_sets[label] & primary_keys) / len(key_sets[label])
            if key_sets[label]
            else 0.0
        )
    )

    overlap_rows = []
    for left, right in combinations(inputs, 2):
        intersection = key_sets[left] & key_sets[right]
        union = key_sets[left] | key_sets[right]
        overlap_rows.append(
            {
                "left_evidence_scope": left,
                "right_evidence_scope": right,
                "left_G1_structures": len(key_sets[left]),
                "right_G1_structures": len(key_sets[right]),
                "intersection_G1_structures": len(intersection),
                "union_G1_structures": len(union),
                "jaccard_similarity": (
                    len(intersection) / len(union) if union else 0.0
                ),
                "left_only_G1_structures": len(
                    key_sets[left] - key_sets[right]
                ),
                "right_only_G1_structures": len(
                    key_sets[right] - key_sets[left]
                ),
            }
        )
    overlap = pd.DataFrame(overlap_rows)

    representative: dict[str, dict[str, str]] = {}
    for label, frame in frames.items():
        for row in frame[
            ["full_inchikey", "connectivity_key", "smiles", "formula"]
        ].to_dict("records"):
            representative.setdefault(str(row["full_inchikey"]), row)
    membership_rows = []
    for key in sorted(set().union(*key_sets.values())):
        row = representative[key]
        memberships = [label for label in inputs if key in key_sets[label]]
        membership_rows.append(
            {
                "full_inchikey": key,
                "connectivity_key": row["connectivity_key"],
                "smiles": row["smiles"],
                "formula": row["formula"],
                **{
                    f"present_in_{label}": label in memberships
                    for label in inputs
                },
                "evidence_scope_membership_count": len(memberships),
                "evidence_scope_memberships": ";".join(memberships),
            }
        )
    membership = pd.DataFrame(membership_rows)

    paths = {
        "tier_summary": output_dir / "G1_evidence_layer_sensitivity_summary.tsv",
        "pairwise_overlap": output_dir / "G1_pairwise_structure_overlap.tsv",
        "structure_membership": output_dir / "G1_structure_layer_membership.tsv",
        "summary": output_dir / "G1_sensitivity_comparison_summary.json",
    }
    write_table(tier_summary, paths["tier_summary"])
    write_table(overlap, paths["pairwise_overlap"])
    write_table(membership, paths["structure_membership"])
    summary = {
        "mode": "full_stereochemistry_G1_exclusive_evidence_layer_comparison",
        "identity": "full_stereochemistry_aware_InChIKey",
        "layers_are_exclusive_not_cumulative": True,
        "G1_union_structures": int(len(membership)),
        "G1_shared_by_all_three_layers": int(
            (membership["evidence_scope_membership_count"] == 3).sum()
        ),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
