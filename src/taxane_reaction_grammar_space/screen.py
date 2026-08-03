from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import pandas as pd

from .chemistry import require_rdkit, stable_hash
from .io import ensure_dir, read_table, write_json, write_table


def _detect_column(columns: list[str], candidates: list[str], label: str) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"Could not detect {label}; tried {candidates}")


def _bit_vector_as_int(bit_vector) -> int:
    return int(bit_vector.ToBitString(), 2)


def load_seed_molecules(
    nodes_path: Path,
) -> tuple[pd.DataFrame, list[tuple[str, Any, int]]]:
    Chem, *_rest = require_rdkit()
    nodes = read_table(nodes_path).fillna("")
    smiles_column = _detect_column(
        list(nodes.columns),
        ["standardized_smiles", "smiles", "smiles_raw"],
        "seed SMILES column",
    )
    id_column = _detect_column(
        list(nodes.columns),
        ["molecule_id", "node_id", "source_id", "row_index"],
        "seed identifier column",
    )
    valid_records: list[dict[str, str]] = []
    molecules: list[tuple[str, Any, int]] = []
    for row in nodes.to_dict("records"):
        seed_id = str(row.get(id_column, "")).strip()
        smiles = str(row.get(smiles_column, "")).strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        fingerprint = Chem.PatternFingerprint(mol, fpSize=2048)
        valid_records.append(
            {
                "seed_id": seed_id,
                "seed_smiles": Chem.MolToSmiles(
                    mol, canonical=True, isomericSmiles=True
                ),
            }
        )
        molecules.append((seed_id, mol, _bit_vector_as_int(fingerprint)))
    if not molecules:
        raise RuntimeError("No valid seed molecules were loaded")
    return pd.DataFrame(valid_records), molecules


def screen_grammar_against_seeds(
    grammar_path: Path,
    nodes_path: Path,
    output_dir: Path,
    *,
    validate_prefilter_rules: int = 100,
    max_site_matches_per_seed: int = 512,
    matched_seed_id_sample_size: int = 25,
) -> dict[str, Path]:
    start = time.time()
    Chem, AllChem, *_rest = require_rdkit()
    output_dir = ensure_dir(output_dir)
    grammar = read_table(grammar_path).fillna("")
    release_tier = "T1"
    if "exclusive_release_tier" in grammar.columns:
        values = {
            str(value).strip().upper()
            for value in grammar["exclusive_release_tier"]
            if str(value).strip()
        }
        if len(values) == 1:
            candidate = next(iter(values)).split("_", 1)[0]
            if candidate in {"T1", "T2", "T3"}:
                release_tier = candidate
    seeds, molecules = load_seed_molecules(nodes_path)
    activated: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    prefilter_validation_failures: list[dict[str, str]] = []

    for index, row in enumerate(grammar.to_dict("records")):
        smarts = str(row.get("reaction_smarts", "")).strip()
        rule_id = str(row.get("grammar_rule_id", row.get("smarts_rule_id", "")))
        status = "ok"
        error = ""
        matched_seed_ids: list[str] = []
        site_match_count = 0
        prefilter_candidates = 0
        try:
            reaction = AllChem.ReactionFromSmarts(smarts)
            if reaction is None or reaction.GetNumReactantTemplates() != 1:
                raise ValueError("expected one reactant template")
            query = reaction.GetReactantTemplate(0)
            query_fp = _bit_vector_as_int(Chem.PatternFingerprint(query, fpSize=2048))
            for seed_id, mol, mol_fp in molecules:
                passes_prefilter = query_fp & ~mol_fp == 0
                direct_match = False
                if passes_prefilter:
                    prefilter_candidates += 1
                    direct_match = mol.HasSubstructMatch(query)
                if index < validate_prefilter_rules:
                    validation_match = mol.HasSubstructMatch(query)
                    if validation_match and not passes_prefilter:
                        prefilter_validation_failures.append(
                            {"grammar_rule_id": rule_id, "seed_id": seed_id}
                        )
                    direct_match = validation_match
                if direct_match:
                    matched_seed_ids.append(seed_id)
                    matches = mol.GetSubstructMatches(
                        query,
                        uniquify=True,
                        maxMatches=max_site_matches_per_seed,
                    )
                    site_match_count += len(matches)
        except Exception as exc:
            status = "compile_or_match_failed"
            error = str(exc)

        audit_record = {
            "grammar_rule_id": rule_id,
            "smarts_rule_id": str(row.get("smarts_rule_id", "")),
            "semantic_group_id": str(row.get("semantic_group_id", "")),
            "reaction_type": str(row.get("reaction_type", "")),
            "screen_status": status,
            "screen_error": error,
            "g0_match_count": len(matched_seed_ids),
            "g0_site_match_count": site_match_count,
            "g0_prefilter_candidate_count": prefilter_candidates,
            "g0_match_seed_ids_sample": ";".join(
                matched_seed_ids[:matched_seed_id_sample_size]
            ),
            "g0_match_seed_set_hash": stable_hash(
                "\x1f".join(matched_seed_ids), length=24
            ),
            "g0_activated": bool(matched_seed_ids),
        }
        audit.append(audit_record)
        if matched_seed_ids and status == "ok":
            row.update(audit_record)
            activated.append(row)
        if (index + 1) % 1000 == 0:
            print(
                f"[screen-grammar] rules={index + 1} activated={len(activated)}",
                flush=True,
            )

    if prefilter_validation_failures:
        raise RuntimeError(
            "Pattern-fingerprint prefilter produced false negatives; "
            f"see {len(prefilter_validation_failures)} validation failures"
        )

    activated_frame = pd.DataFrame(activated)
    audit_frame = pd.DataFrame(audit)
    paths = {
        "activated": output_dir
        / f"generative_grammar.{release_tier}_G0_activated.tsv",
        "audit": output_dir
        / f"generative_grammar.{release_tier}_G0_screen_audit.tsv",
        "seeds": output_dir / "G0_seed_inventory.tsv",
        "summary": output_dir
        / f"generative_grammar.{release_tier}_G0_screen_summary.json",
    }
    write_table(activated_frame, paths["activated"])
    write_table(audit_frame, paths["audit"])
    write_table(seeds, paths["seeds"])
    summary = {
        "mode": "G0_substructure_compatibility_screen",
        "release_tier": release_tier,
        "grammar_input": str(grammar_path),
        "nodes_input": str(nodes_path),
        "grammar_rules": int(len(grammar)),
        "valid_G0_seeds": int(len(seeds)),
        "activated_rules": int(len(activated_frame)),
        "activated_semantic_groups": int(
            activated_frame.get(
                "semantic_group_id", pd.Series(dtype=str)
            ).nunique()
        ),
        "total_G0_rule_matches": int(
            audit_frame.get("g0_match_count", pd.Series(dtype=int)).sum()
        ),
        "total_G0_site_matches": int(
            audit_frame.get("g0_site_match_count", pd.Series(dtype=int)).sum()
        ),
        "compile_or_match_failures": int(
            (audit_frame["screen_status"] != "ok").sum()
        ),
        "prefilter_validation_rules": min(
            validate_prefilter_rules, int(len(grammar))
        ),
        "prefilter_validation_false_negatives": 0,
        "elapsed_seconds": round(time.time() - start, 3),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
