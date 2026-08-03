#!/usr/bin/env python3
"""Recompute route-layer outputs without repeating structure-level analyses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.analyze import _convergence_and_paths
from taxane_reaction_grammar_space.io import (
    ensure_dir,
    read_table,
    write_json,
    write_table,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = ensure_dir(args.output_dir)
    nodes = read_table(args.nodes).fillna("")
    events = read_table(args.events).fillna("")
    nodes["generation_first"] = pd.to_numeric(
        nodes["generation_first"], errors="raise"
    ).astype(int)
    for column in (
        "event_id",
        "generation",
        "target_is_new",
        "target_generation_first",
    ):
        events[column] = pd.to_numeric(
            events[column], errors="raise"
        ).astype(int)

    convergence, bridges, bridge_pairs = _convergence_and_paths(nodes, events)
    convergence_path = output_dir / "convergence_and_route_multiplicity.tsv"
    bridges_path = output_dir / "latent_bridge_candidates.tsv"
    bridge_pairs_path = output_dir / "known_G0_pair_bridge_summary.tsv"
    write_table(convergence, convergence_path)
    write_table(bridges, bridges_path)
    write_table(bridge_pairs, bridge_pairs_path)

    source_summary_path = (
        args.analysis_dir / "chemical_space_analysis_summary.json"
    )
    summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
    summary.update(
        {
            "convergent_generated_nodes": int(
                convergence.loc[
                    convergence["generation"] > 0, "is_convergent"
                ].astype(bool).sum()
            ),
            "latent_bridge_candidates": int(
                bridges["latent_bridge_candidate"].astype(bool).sum()
            ),
            "convergence_definition": (
                "at_least_two_distinct_parent_structures_in_the_targets_"
                "first_observed_generation"
            ),
            "path_count_layers": {
                "primary": "distinct_source_target_structural_edges",
                "semantic_audit": (
                    "distinct_source_target_semantic_group_edges"
                ),
                "raw_audit": "all_parent_rule_product_events",
            },
            "path_event_inclusion": (
                "all_events_recorded_in_the_targets_first_observed_"
                "generation_regardless_of_insertion_flag"
            ),
        }
    )
    summary["outputs"].update(
        {
            "convergence": str(convergence_path),
            "bridges": str(bridges_path),
            "bridge_pairs": str(bridge_pairs_path),
            "summary": str(
                output_dir / "chemical_space_analysis_summary.json"
            ),
        }
    )
    write_json(
        summary, output_dir / "chemical_space_analysis_summary.json"
    )
    write_json(
        {
            "status": "complete",
            "input_nodes": int(len(nodes)),
            "input_derivation_events": int(len(events)),
            "convergent_generated_nodes": int(
                summary["convergent_generated_nodes"]
            ),
            "latent_bridge_candidates": int(
                summary["latent_bridge_candidates"]
            ),
            "outputs": {
                "convergence": str(convergence_path),
                "bridges": str(bridges_path),
                "bridge_pairs": str(bridge_pairs_path),
            },
        },
        output_dir / "route_layer_recompute_summary.json",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
