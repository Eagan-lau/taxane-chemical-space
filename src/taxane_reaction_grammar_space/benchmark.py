from __future__ import annotations

import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .chemistry import formula_delta, require_rdkit, stable_hash
from .generate import _clear_atom_maps, _compile_rules, _fingerprint, _molecule_identity
from .io import ensure_dir, read_table, write_json, write_table
from .screen import _bit_vector_as_int


def _text(record: dict[str, Any], key: str) -> str:
    value = record.get(key, "")
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() == "nan" else value


def _is_taxol_derived(record: dict[str, Any]) -> bool:
    searchable = " ".join(
        _text(record, key)
        for key in (
            "template_sources",
            "benchmark_exclusion_flag",
            "curated_taxol_anchor",
            "curated_pathway_name",
            "curated_pathway_step_ids",
        )
    ).lower()
    return "taxol" in searchable or "curated_taxol" in searchable


def _example_pair_keys(record: dict[str, Any]) -> tuple[str, str] | None:
    Chem, *_rest = require_rdkit()
    substrate = Chem.MolFromSmiles(_text(record, "example_substrate_smiles"))
    product = Chem.MolFromSmiles(_text(record, "example_product_smiles"))
    if substrate is None or product is None:
        reaction = _text(record, "example_reaction_smiles")
        if ">>" in reaction:
            left, right = reaction.split(">>", 1)
            substrate = Chem.MolFromSmiles(left)
            product = Chem.MolFromSmiles(right)
    if substrate is None or product is None:
        return None
    try:
        return (
            _molecule_identity(substrate)["connectivity_key"],
            _molecule_identity(product)["connectivity_key"],
        )
    except ValueError:
        return None


def _benchmark_group(source, target) -> str:
    from rdkit import DataStructs

    source_identity = _molecule_identity(source)
    target_identity = _molecule_identity(target)
    source_fp = _fingerprint(source)
    target_fp = _fingerprint(target)
    xor_fp = source_fp ^ target_fp
    payload = {
        "element_delta": formula_delta(
            source_identity["formula"], target_identity["formula"]
        ),
        "fingerprint_xor": DataStructs.BitVectToBinaryText(xor_fp).hex(),
        "heavy_atom_delta": (
            target_identity["heavy_atom_count"] - source_identity["heavy_atom_count"]
        ),
        "ring_delta": target_identity["ring_count"] - source_identity["ring_count"],
    }
    return stable_hash(json.dumps(payload, sort_keys=True), length=20)


def _apply_rules_to_source(
    source,
    compiled_rules,
    example_pair_by_rule: dict[str, tuple[str, str] | None],
    *,
    benchmark_pair: tuple[str, str] | None,
    max_products_per_rule: int,
) -> tuple[dict[str, list[str]], dict[str, list[str]], int]:
    Chem, *_rest = require_rdkit()
    source_pattern_fp = _bit_vector_as_int(
        Chem.PatternFingerprint(source, fpSize=2048)
    )
    full_hits: dict[str, list[str]] = defaultdict(list)
    connectivity_hits: dict[str, list[str]] = defaultdict(list)
    leakage_excluded = 0
    for rule in compiled_rules:
        if benchmark_pair is not None:
            example_pair = example_pair_by_rule.get(rule.grammar_rule_id)
            if example_pair == benchmark_pair:
                leakage_excluded += 1
                continue
        if rule.query_fp & ~source_pattern_fp:
            continue
        if not source.HasSubstructMatch(rule.query):
            continue
        try:
            outcomes = rule.reaction.RunReactants(
                (source,), maxProducts=max_products_per_rule + 1
            )
        except Exception:
            continue
        if len(outcomes) > max_products_per_rule:
            continue
        for outcome in outcomes:
            if len(outcome) != 1:
                continue
            product = outcome[0]
            try:
                Chem.SanitizeMol(product)
                product = Chem.RemoveHs(product)
                if len(Chem.GetMolFrags(product)) != 1:
                    continue
                _clear_atom_maps(product)
                identity = _molecule_identity(product)
            except Exception:
                continue
            full_hits[identity["full_inchikey"]].append(rule.grammar_rule_id)
            connectivity_hits[identity["connectivity_key"]].append(
                rule.grammar_rule_id
            )
    return full_hits, connectivity_hits, leakage_excluded


def _wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return (float("nan"), float("nan"))
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half = (
        z
        * math.sqrt(
            proportion * (1 - proportion) / total + z * z / (4 * total * total)
        )
        / denominator
    )
    return (center - half, center + half)


def _fisher_exact(table: list[list[int]]) -> tuple[float, float]:
    try:
        from scipy.stats import fisher_exact

        result = fisher_exact(table)
        return float(result.statistic), float(result.pvalue)
    except Exception:
        return float("nan"), float("nan")


def _load_benchmark_reactions(path: Path) -> list[dict[str, Any]]:
    Chem, *_rest = require_rdkit()
    frame = read_table(path).fillna("")
    substrate_column = next(
        column
        for column in ("Substrate", "substrate", "source_smiles")
        if column in frame.columns
    )
    product_column = next(
        column
        for column in ("Product", "product", "target_smiles")
        if column in frame.columns
    )
    enzyme_column = next(
        (column for column in ("Enzyme", "enzyme", "step_id") if column in frame.columns),
        None,
    )
    reactions: list[dict[str, Any]] = []
    for index, row in enumerate(frame.to_dict("records"), start=1):
        source = Chem.MolFromSmiles(_text(row, substrate_column))
        target = Chem.MolFromSmiles(_text(row, product_column))
        if source is None or target is None:
            reactions.append(
                {
                    "benchmark_id": f"TAXOL_{index:03d}",
                    "valid": False,
                    "error": "invalid_substrate_or_product_smiles",
                }
            )
            continue
        source_identity = _molecule_identity(source)
        target_identity = _molecule_identity(target)
        reactions.append(
            {
                "benchmark_id": f"TAXOL_{index:03d}",
                "step_label": _text(row, enzyme_column) if enzyme_column else "",
                "valid": True,
                "source_mol": source,
                "target_mol": target,
                "source_smiles": source_identity["smiles"],
                "target_smiles": target_identity["smiles"],
                "source_full_inchikey": source_identity["full_inchikey"],
                "target_full_inchikey": target_identity["full_inchikey"],
                "source_connectivity_key": source_identity["connectivity_key"],
                "target_connectivity_key": target_identity["connectivity_key"],
                "source_formula": source_identity["formula"],
                "target_formula": target_identity["formula"],
                "mass_delta": round(
                    target_identity["exact_mass"] - source_identity["exact_mass"], 6
                ),
                "element_delta": json.dumps(
                    formula_delta(
                        source_identity["formula"], target_identity["formula"]
                    ),
                    sort_keys=True,
                ),
                "benchmark_group_id": _benchmark_group(source, target),
            }
        )
    return reactions


def _load_known_nodes(nodes_path: Path) -> list[dict[str, Any]]:
    Chem, *_rest = require_rdkit()
    frame = read_table(nodes_path).fillna("")
    smiles_column = next(
        column
        for column in ("standardized_smiles", "smiles", "smiles_raw")
        if column in frame.columns
    )
    id_column = next(
        column
        for column in ("molecule_id", "node_id", "source_id", "row_index")
        if column in frame.columns
    )
    unique: dict[str, dict[str, Any]] = {}
    for row in frame.to_dict("records"):
        mol = Chem.MolFromSmiles(_text(row, smiles_column))
        if mol is None:
            continue
        identity = _molecule_identity(mol)
        if identity["full_inchikey"] in unique:
            continue
        unique[identity["full_inchikey"]] = {
            "node_id": _text(row, id_column),
            "mol": mol,
            **identity,
            "fingerprint": _fingerprint(mol),
        }
    return list(unique.values())


def _choose_decoys(
    reactions: list[dict[str, Any]],
    known_nodes: list[dict[str, Any]],
    *,
    decoys_per_positive: int,
) -> pd.DataFrame:
    from rdkit import DataStructs

    rows: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str]] = set()
    positive_pairs = {
        (reaction["source_connectivity_key"], reaction["target_connectivity_key"])
        for reaction in reactions
        if reaction.get("valid")
    }
    for reaction in reactions:
        if not reaction.get("valid"):
            continue
        source_fp = _fingerprint(reaction["source_mol"])
        target_fp = _fingerprint(reaction["target_mol"])
        positive_similarity = DataStructs.TanimotoSimilarity(source_fp, target_fp)
        target_identity = _molecule_identity(reaction["target_mol"])
        candidates = []
        for node in known_nodes:
            pair = (reaction["source_connectivity_key"], node["connectivity_key"])
            if pair in positive_pairs or pair in seen_pairs:
                continue
            if node["connectivity_key"] == reaction["target_connectivity_key"]:
                continue
            source_similarity = DataStructs.TanimotoSimilarity(
                source_fp, node["fingerprint"]
            )
            target_similarity = DataStructs.TanimotoSimilarity(
                target_fp, node["fingerprint"]
            )
            formula_penalty = 0 if node["formula"] == target_identity["formula"] else 1
            mass_penalty = abs(node["exact_mass"] - target_identity["exact_mass"])
            similarity_penalty = abs(source_similarity - positive_similarity)
            score = (
                formula_penalty,
                round(mass_penalty, 6),
                round(similarity_penalty, 6),
                -round(target_similarity, 6),
                node["full_inchikey"],
            )
            candidates.append((score, node, source_similarity, target_similarity))
        candidates.sort(key=lambda item: item[0])
        for rank, (_score, node, source_similarity, target_similarity) in enumerate(
            candidates[:decoys_per_positive], start=1
        ):
            pair = (reaction["source_connectivity_key"], node["connectivity_key"])
            seen_pairs.add(pair)
            rows.append(
                {
                    "benchmark_id": reaction["benchmark_id"],
                    "benchmark_group_id": reaction["benchmark_group_id"],
                    "decoy_rank": rank,
                    "source_connectivity_key": reaction["source_connectivity_key"],
                    "decoy_node_id": node["node_id"],
                    "decoy_full_inchikey": node["full_inchikey"],
                    "decoy_connectivity_key": node["connectivity_key"],
                    "decoy_formula": node["formula"],
                    "formula_exact_to_true_target": (
                        node["formula"] == target_identity["formula"]
                    ),
                    "absolute_mass_difference_to_true_target": round(
                        abs(node["exact_mass"] - target_identity["exact_mass"]), 6
                    ),
                    "source_decoy_tanimoto": round(float(source_similarity), 6),
                    "true_target_decoy_tanimoto": round(float(target_similarity), 6),
                }
            )
    return pd.DataFrame(rows)


def run_leakage_controlled_benchmark(
    activated_grammar_path: Path,
    taxol_pathway_path: Path,
    nodes_path: Path,
    output_dir: Path,
    *,
    decoys_per_positive: int = 20,
    max_products_per_rule: int = 256,
    exclude_taxol_derived: bool = True,
) -> dict[str, Path]:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.error")
    output_dir = ensure_dir(output_dir)
    grammar = read_table(activated_grammar_path).fillna("")
    if exclude_taxol_derived:
        evaluated = grammar[
            ~grammar.apply(lambda row: _is_taxol_derived(row.to_dict()), axis=1)
        ].copy()
    else:
        evaluated = grammar.copy()
    compiled, compile_failures = _compile_rules(evaluated)
    record_by_rule = {
        (_text(record, "grammar_rule_id") or _text(record, "smarts_rule_id")): record
        for record in evaluated.to_dict("records")
    }
    example_pair_by_rule = {
        rule_id: _example_pair_keys(record)
        for rule_id, record in record_by_rule.items()
    }
    reactions = _load_benchmark_reactions(taxol_pathway_path)
    known_nodes = _load_known_nodes(nodes_path)
    decoys = _choose_decoys(
        reactions, known_nodes, decoys_per_positive=decoys_per_positive
    )
    decoys_by_benchmark = {
        key: group.copy()
        for key, group in decoys.groupby("benchmark_id", sort=False)
    }

    positive_rows: list[dict[str, Any]] = []
    decoy_rows: list[dict[str, Any]] = []
    hit_rows: list[dict[str, Any]] = []
    for reaction in reactions:
        if not reaction.get("valid"):
            positive_rows.append(reaction)
            continue
        pair = (
            reaction["source_connectivity_key"],
            reaction["target_connectivity_key"],
        )
        full_hits, connectivity_hits, leakage_excluded = _apply_rules_to_source(
            reaction["source_mol"],
            compiled,
            example_pair_by_rule,
            benchmark_pair=pair,
            max_products_per_rule=max_products_per_rule,
        )
        full_rule_ids = sorted(
            set(full_hits.get(reaction["target_full_inchikey"], []))
        )
        connectivity_rule_ids = sorted(
            set(connectivity_hits.get(reaction["target_connectivity_key"], []))
        )
        recovered_full = bool(full_rule_ids)
        recovered_connectivity = bool(connectivity_rule_ids)
        positive_row = {
            key: value
            for key, value in reaction.items()
            if key not in {"source_mol", "target_mol"}
        }
        positive_row.update(
            {
                "recovered_full_stereo": recovered_full,
                "recovered_connectivity": recovered_connectivity,
                "full_match_rule_count": len(full_rule_ids),
                "connectivity_match_rule_count": len(connectivity_rule_ids),
                "full_match_rule_ids": ";".join(full_rule_ids),
                "connectivity_match_rule_ids": ";".join(connectivity_rule_ids),
                "exact_pair_leakage_rules_excluded": leakage_excluded,
            }
        )
        positive_rows.append(positive_row)
        for rule_id in connectivity_rule_ids:
            hit_rows.append(
                {
                    "benchmark_id": reaction["benchmark_id"],
                    "target_type": "positive",
                    "target_connectivity_key": reaction["target_connectivity_key"],
                    "grammar_rule_id": rule_id,
                    "full_stereo_match": rule_id in full_rule_ids,
                }
            )

        decoy_group = decoys_by_benchmark.get(reaction["benchmark_id"])
        if decoy_group is None:
            continue
        for decoy in decoy_group.to_dict("records"):
            full_rule_ids = sorted(
                set(full_hits.get(decoy["decoy_full_inchikey"], []))
            )
            connectivity_rule_ids = sorted(
                set(connectivity_hits.get(decoy["decoy_connectivity_key"], []))
            )
            decoy.update(
                {
                    "matched_full_stereo": bool(full_rule_ids),
                    "matched_connectivity": bool(connectivity_rule_ids),
                    "full_match_rule_count": len(full_rule_ids),
                    "connectivity_match_rule_count": len(connectivity_rule_ids),
                    "full_match_rule_ids": ";".join(full_rule_ids),
                    "connectivity_match_rule_ids": ";".join(connectivity_rule_ids),
                }
            )
            decoy_rows.append(decoy)
            for rule_id in connectivity_rule_ids:
                hit_rows.append(
                    {
                        "benchmark_id": reaction["benchmark_id"],
                        "target_type": "decoy",
                        "target_connectivity_key": decoy["decoy_connectivity_key"],
                        "grammar_rule_id": rule_id,
                        "full_stereo_match": rule_id in full_rule_ids,
                    }
                )

    positives = pd.DataFrame(positive_rows)
    decoy_results = pd.DataFrame(decoy_rows)
    hits = pd.DataFrame(hit_rows)
    valid_positives = positives[positives.get("valid", False).astype(bool)].copy()
    positive_successes = int(valid_positives["recovered_connectivity"].astype(bool).sum())
    positive_total = int(len(valid_positives))
    decoy_successes = int(decoy_results["matched_connectivity"].astype(bool).sum())
    decoy_total = int(len(decoy_results))
    odds_ratio, fisher_p = _fisher_exact(
        [
            [positive_successes, positive_total - positive_successes],
            [decoy_successes, decoy_total - decoy_successes],
        ]
    )
    positive_interval = _wilson_interval(positive_successes, positive_total)
    decoy_interval = _wilson_interval(decoy_successes, decoy_total)
    group_summary = (
        valid_positives.groupby("benchmark_group_id", dropna=False)
        .agg(
            reaction_count=("benchmark_id", "size"),
            connectivity_recovered=("recovered_connectivity", "sum"),
            full_stereo_recovered=("recovered_full_stereo", "sum"),
            step_labels=("step_label", lambda values: ";".join(sorted(set(values)))),
        )
        .reset_index()
    )
    paths = {
        "positives": output_dir / "benchmark_positive_recovery.tsv",
        "decoys": output_dir / "benchmark_matched_decoys.tsv",
        "hits": output_dir / "benchmark_rule_hits.tsv",
        "groups": output_dir / "benchmark_reaction_center_groups.tsv",
        "compile_failures": output_dir / "benchmark_grammar_compile_failures.tsv",
        "summary": output_dir / "benchmark_summary.json",
    }
    write_table(positives, paths["positives"])
    write_table(decoy_results, paths["decoys"])
    write_table(hits, paths["hits"])
    write_table(group_summary, paths["groups"])
    write_table(pd.DataFrame(compile_failures), paths["compile_failures"])
    summary = {
        "mode": (
            "external_evidence_exact_pair_leakage_controlled_recovery"
            if exclude_taxol_derived
            else "domain_informed_pathway_replay_calibration"
        ),
        "activated_grammar_input": str(activated_grammar_path),
        "taxol_pathway_input": str(taxol_pathway_path),
        "known_nodes_input": str(nodes_path),
        "grammar_rules_before_taxol_exclusion": int(len(grammar)),
        "exclude_taxol_derived": bool(exclude_taxol_derived),
        "evaluated_grammar_rules": int(len(evaluated)),
        "external_grammar_rules": (
            int(len(evaluated)) if exclude_taxol_derived else None
        ),
        "compiled_external_rules": (
            int(len(compiled)) if exclude_taxol_derived else None
        ),
        "compiled_evaluated_rules": int(len(compiled)),
        "compile_failures": int(len(compile_failures)),
        "positive_reactions": positive_total,
        "positive_connectivity_recovered": positive_successes,
        "positive_connectivity_recovery_rate": (
            positive_successes / positive_total if positive_total else None
        ),
        "positive_recovery_wilson_95ci": list(positive_interval),
        "positive_full_stereo_recovered": int(
            valid_positives["recovered_full_stereo"].astype(bool).sum()
        ),
        "decoy_pairs": decoy_total,
        "decoy_connectivity_matched": decoy_successes,
        "decoy_connectivity_match_rate": (
            decoy_successes / decoy_total if decoy_total else None
        ),
        "decoy_match_wilson_95ci": list(decoy_interval),
        "fisher_exact_odds_ratio": odds_ratio,
        "fisher_exact_p_value": fisher_p,
        "reaction_center_group_count": int(
            valid_positives["benchmark_group_id"].nunique()
        ),
        "leakage_controls": (
            [
                "all TaxolKnownPathway-derived rules excluded globally",
                "any external rule carrying the exact benchmark connectivity pair excluded per reaction",
                "recovery summarized by reaction-center fingerprint group",
                "connectivity and full-stereochemistry recovery reported separately",
            ]
            if exclude_taxol_derived
            else [
                "domain-derived rules deliberately retained as an internal replay calibration",
                "no claim of independent pathway prediction is made for this mode",
                "connectivity and full-stereochemistry replay reported separately",
            ]
        ),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
