from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import pandas as pd

from .abstraction import (
    estimate_rxnmapper_token_count,
    extract_rdchiral_template,
    map_reaction_smiles_batch,
)
from .chem import reaction_smarts_valid, replay_reaction_smarts_on_pair
from .release import export_reaction_smarts_library, validate_smarts_release
from .rules import build_rule_library
from .schema import RULE_COLUMNS, SMARTS_LIBRARY_COLUMNS, TEMPLATE_COLUMNS
from .templates import deduplicate_templates, qc_templates
from .utils import clean_text, ensure_dir, join_values, sha256_text, write_json


ANCHOR_GENERALIZATION_REPORT_COLUMNS = [
    "record_id",
    "source_database",
    "source_reaction_id",
    "evidence_layer",
    "canonical_reaction_smiles",
    "main_substrate_smiles",
    "main_product_smiles",
    "mapped_reaction_smiles",
    "rxnmapper_confidence",
    "extraction_status",
    "reaction_smarts",
    "template_hash",
    "template_qc_status",
    "smarts_replay_pass",
    "smarts_replay_note",
    "template_qc_note",
]


def _mode_suffix(release_mode: str) -> str:
    return "" if clean_text(release_mode).lower() == "strict" else f".{clean_text(release_mode).lower()}"


def _smarts_release_stem(release_mode: str) -> str:
    mode = clean_text(release_mode).lower()
    return "anchor_derived" if mode == "strict" else f"anchor_derived_{mode}"


def _empty_template_df() -> pd.DataFrame:
    return pd.DataFrame(columns=TEMPLATE_COLUMNS)


def _empty_rule_df() -> pd.DataFrame:
    return pd.DataFrame(columns=RULE_COLUMNS)


def _empty_smarts_df() -> pd.DataFrame:
    return pd.DataFrame(columns=SMARTS_LIBRARY_COLUMNS)


def _text_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([""] * len(df), index=df.index, dtype=str)
    return df[col].astype(str)


def _multi_value_contains_any(value: Any, enabled: set[str]) -> bool:
    text = clean_text(value)
    if not text:
        return False
    parts = [p.strip() for p in text.replace("|", ";").replace(",", ";").split(";") if p.strip()]
    return any(p in enabled for p in parts)


def _anchor_candidate_mask(templates: pd.DataFrame, enabled_layers: set[str]) -> pd.Series:
    if templates.empty:
        return pd.Series([], dtype=bool)
    scope_ok = _text_series(templates, "template_scope").eq("exact_anchor")
    anchor_ok = _text_series(templates, "anchor_edge_use").str.lower().isin(["true", "1", "yes"])
    rxn_ok = _text_series(templates, "canonical_reaction_smiles").str.contains(">>", regex=False)
    layer_ok = _text_series(templates, "evidence_layer").map(lambda x: _multi_value_contains_any(x, enabled_layers))
    return scope_ok & anchor_ok & rxn_ok & layer_ok


def _report_row(
    row: dict[str, Any],
    mapped_rxn: str,
    confidence: float,
    status: str,
    smarts: str,
    template_hash: str,
    qc_status: str,
    replay_pass: str,
    replay_note: str,
    qc_note: str,
) -> dict[str, Any]:
    return {
        "record_id": clean_text(row.get("record_id")),
        "source_database": clean_text(row.get("source_database")),
        "source_reaction_id": clean_text(row.get("source_reaction_id")),
        "evidence_layer": clean_text(row.get("evidence_layer")),
        "canonical_reaction_smiles": clean_text(row.get("canonical_reaction_smiles")),
        "main_substrate_smiles": clean_text(row.get("main_substrate_smiles")),
        "main_product_smiles": clean_text(row.get("main_product_smiles")),
        "mapped_reaction_smiles": clean_text(mapped_rxn),
        "rxnmapper_confidence": f"{float(confidence):.4f}" if confidence else "",
        "extraction_status": clean_text(status),
        "reaction_smarts": clean_text(smarts),
        "template_hash": clean_text(template_hash),
        "template_qc_status": clean_text(qc_status),
        "smarts_replay_pass": clean_text(replay_pass),
        "smarts_replay_note": clean_text(replay_note),
        "template_qc_note": clean_text(qc_note),
    }


def _promote_anchor_to_template(
    row: dict[str, Any],
    smarts: str,
    confidence: float,
    template_index: int,
    *,
    release_mode: str = "strict",
    acceptance_status: str = "ok",
    replay_ok: bool = True,
    replay_note: str = "replay_product_match",
    mapper_confidence_ok: bool = True,
    qc_note: str = "ok",
) -> dict[str, Any]:
    mode = clean_text(release_mode).lower() or "strict"
    status = clean_text(acceptance_status) or "ok"
    out = dict(row)
    out["reaction_smarts"] = clean_text(smarts)
    prefix = "ANCGEN_TPL" if mode == "strict" else f"ANCGEN_{mode.upper()}_TPL"
    out["template_id"] = f"{prefix}_{template_index:09d}"
    out["template_hash"] = sha256_text(smarts)[:20]
    out["template_origin"] = "exact_anchor_derived_smarts" if mode == "strict" else f"exact_anchor_derived_smarts_{mode}"
    out["template_scope"] = "generalized_template"
    out["template_generalization"] = "anchor_generalization_rxnmapper_rdchiral" if mode == "strict" else f"anchor_generalization_rxnmapper_rdchiral_{mode}"
    out["source_reaction_smarts"] = clean_text(smarts)
    out["predictive_rule_use"] = "true"
    out["anchor_edge_use"] = "false"
    out["template_extraction_status"] = status
    out["template_qc_status"] = "ok"
    out["template_qc_note"] = join_values(
        [
            "anchor_generalization_replay_validated_against_original_pair" if replay_ok else "anchor_generalization_replay_failed_kept_permissive",
            "rxnmapper_confidence_pass" if mapper_confidence_ok else "rxnmapper_low_confidence_kept_permissive",
            replay_note,
            qc_note,
        ]
    )
    out["abstracted_from_exact_reaction"] = "true"
    out["derived_from_exact_anchor"] = "true"
    out["rxnmapper_confidence"] = f"{float(confidence):.4f}"
    out["rdchiral_extraction_status"] = "ok"
    out["abstracted_smarts_applies_to_original_pair"] = "true" if replay_ok else "false"
    out["exact_abstraction_qc_status"] = "pass" if replay_ok and mapper_confidence_ok else status
    out["benchmark_exclusion_flag"] = (
        "anchor_derived_exclude_from_primary_benchmark"
        if mode == "strict"
        else f"anchor_derived_{mode}_{status}_exclude_from_core"
    )
    out["source_evidence_text"] = join_values(
        [
            out.get("source_evidence_text", ""),
            "exact_anchor_generalized_by_rxnmapper_rdchiral",
            "smarts_replay_validated_against_original_pair" if replay_ok else "smarts_replay_failed_retained_as_permissive_rule",
            "rxnmapper_confidence_pass" if mapper_confidence_ok else "rxnmapper_low_confidence_retained_as_permissive_rule",
            f"anchor_generalization_mode={mode}",
            f"anchor_generalization_status={status}",
            f"rxnmapper_confidence={float(confidence):.4f}",
        ]
    )
    return out


def _write_empty_outputs(release_dir: Path, report: pd.DataFrame, summary: dict[str, Any], *, release_mode: str = "strict") -> dict[str, Any]:
    suffix = _mode_suffix(release_mode)
    stem = _smarts_release_stem(release_mode)
    report.to_csv(release_dir / f"anchor_generalization{suffix}.report.tsv", sep="\t", index=False)
    _empty_template_df().to_csv(release_dir / f"anchor_generalization{suffix}.templates.raw.tsv", sep="\t", index=False)
    _empty_template_df().to_csv(release_dir / f"anchor_generalization{suffix}.templates.deduplicated.tsv", sep="\t", index=False)
    _empty_rule_df().to_csv(release_dir / f"anchor_generalization{suffix}.rules.audit.tsv", sep="\t", index=False)
    _empty_smarts_df().to_csv(release_dir / f"reaction_smarts_library.{stem}.tsv", sep="\t", index=False)
    _empty_smarts_df().to_csv(release_dir / f"reaction_smarts_rules.{stem}.tsv", sep="\t", index=False)
    issues, validation = validate_smarts_release(_empty_smarts_df())
    issues.to_csv(release_dir / f"anchor_generalization{suffix}.validation_issues.tsv", sep="\t", index=False)
    write_json(release_dir / f"anchor_generalization{suffix}.validation.json", validation)
    summary = {
        **summary,
        "raw_anchor_derived_templates": 0,
        "deduplicated_anchor_derived_templates": 0,
        "anchor_derived_audit_rules": 0,
        "anchor_derived_released_smarts": 0,
        "validation": validation,
    }
    write_json(release_dir / f"anchor_generalization{suffix}.summary.json", summary)
    return summary


def generalize_exact_anchor_templates(
    templates: pd.DataFrame,
    family_evidence: pd.DataFrame,
    release_dir: str | Path,
    *,
    enabled_layers: list[str] | None = None,
    batch_size: int = 16,
    min_mapper_confidence: float = 0.50,
    max_records: int | None = None,
    max_reaction_chars: int = 2000,
    max_mapper_tokens: int = 500,
    mapper_timeout_seconds: int = 120,
    release_mode: str = "strict",
) -> dict[str, Any]:
    """Generalize exact anchor pairs into a separate anchor-derived SMARTS release.

    This flow intentionally does not mutate the main template/rule tables and does
    not add anchor-derived SMARTS to the strict core release. It is an exploratory
    supplement with explicit replay validation and failure reporting.
    """
    release = ensure_dir(release_dir)
    mode = clean_text(release_mode).lower() or "strict"
    if mode not in {"strict", "permissive"}:
        raise ValueError(f"Unknown anchor generalization release_mode: {release_mode!r}")
    permissive = mode == "permissive"
    suffix = _mode_suffix(mode)
    stem = _smarts_release_stem(mode)
    layers = set(enabled_layers or ["T1_Bio_Core", "T2_Bio_Extended", "T3_Chem_like"])
    df = templates.copy().fillna("")
    cand_mask = _anchor_candidate_mask(df, layers)
    candidate_idx = list(df.index[cand_mask])
    if max_records is not None:
        candidate_idx = candidate_idx[: int(max_records)]

    start = time.time()
    base_summary: dict[str, Any] = {
        "enabled": True,
        "candidate_anchor_records": int(cand_mask.sum()),
        "selected_anchor_records": int(len(candidate_idx)),
        "attempted_records": 0,
        "successful_templates": 0,
        "failed_records": 0,
        "enabled_layers": sorted(layers),
        "min_mapper_confidence": float(min_mapper_confidence),
        "max_reaction_chars": int(max_reaction_chars) if max_reaction_chars else None,
        "max_mapper_tokens": int(max_mapper_tokens) if max_mapper_tokens else None,
        "mapper_timeout_seconds": int(mapper_timeout_seconds) if mapper_timeout_seconds else None,
        "release_mode": mode,
        "accept_replay_failed": bool(permissive),
        "accept_low_mapper_confidence": bool(permissive),
        "release_policy": "separate_anchor_derived_release_not_mixed_into_core" if mode == "strict" else "permissive_all_extractable_anchor_smarts_annotated_not_for_core",
    }
    if not candidate_idx:
        report = pd.DataFrame(columns=ANCHOR_GENERALIZATION_REPORT_COLUMNS)
        return _write_empty_outputs(release, report, {**base_summary, "elapsed_seconds": 0.0}, release_mode=mode)

    try:
        from rxnmapper import RXNMapper
    except Exception as exc:  # pragma: no cover - depends on optional environment
        reports = [
            _report_row(df.loc[idx].to_dict(), "", 0.0, "rxnmapper_unavailable", "", "", "not_attempted", "false", "", str(exc))
            for idx in candidate_idx
        ]
        report = pd.DataFrame(reports, columns=ANCHOR_GENERALIZATION_REPORT_COLUMNS)
        summary = {
            **base_summary,
            "failed_records": int(len(candidate_idx)),
            "optional_dependency_status": "rxnmapper_unavailable",
            "optional_dependency_error": str(exc),
            "elapsed_seconds": round(time.time() - start, 3),
        }
        return _write_empty_outputs(release, report, summary, release_mode=mode)

    mapper = RXNMapper()
    reports: list[dict[str, Any]] = []
    promoted: list[dict[str, Any]] = []
    attempted = successful = failed = 0
    replay_failed_accepted = 0
    low_confidence_accepted = 0

    for offset in range(0, len(candidate_idx), max(1, batch_size)):
        batch_idx = candidate_idx[offset: offset + max(1, batch_size)]
        if offset == 0 or offset % max(1, batch_size * 100) == 0:
            print(f"[anchor_generalization] processed_candidates={offset}/{len(candidate_idx)}", flush=True)
        batch_rows = [df.loc[i].to_dict() for i in batch_idx]
        batch_rxns: list[str] = []
        valid_positions: list[int] = []

        for pos, row in enumerate(batch_rows):
            rxn = clean_text(row.get("canonical_reaction_smiles"))
            token_count = estimate_rxnmapper_token_count(rxn)
            too_long_chars = bool(max_reaction_chars and len(rxn) > max_reaction_chars)
            too_many_tokens = bool(max_mapper_tokens and token_count > max_mapper_tokens)
            if not rxn or ">>" not in rxn or too_long_chars or too_many_tokens:
                failed += 1
                note = join_values(
                    [
                        "missing_or_invalid_reaction_smiles" if (not rxn or ">>" not in rxn) else "",
                        f"reaction_chars={len(rxn)}",
                        f"max_reaction_chars={max_reaction_chars}",
                        f"estimated_rxnmapper_tokens={token_count}",
                        f"max_mapper_tokens={max_mapper_tokens}",
                        "skipped_before_rxnmapper_to_avoid_transformer_length_failure" if too_long_chars or too_many_tokens else "",
                    ]
                )
                reports.append(_report_row(row, "", 0.0, "invalid_or_too_long_exact_anchor", "", "", "not_attempted", "false", "", note))
                continue
            batch_rxns.append(rxn)
            valid_positions.append(pos)
            attempted += 1

        mapped_results = map_reaction_smiles_batch(mapper, batch_rxns, timeout_seconds=mapper_timeout_seconds) if batch_rxns else []
        if len(mapped_results) < len(valid_positions):
            mapped_results = list(mapped_results) + [
                {"mapped_rxn": "", "confidence": 0.0, "mapper_error": "mapper_returned_fewer_results_than_inputs"}
                for _ in range(len(valid_positions) - len(mapped_results))
            ]
        for pos, mapped_row in zip(valid_positions, mapped_results):
            row = batch_rows[pos]
            mapped_rxn = clean_text(mapped_row.get("mapped_rxn"))
            confidence = float(mapped_row.get("confidence") or 0.0)
            mapper_confidence_ok = confidence >= min_mapper_confidence
            if not mapped_rxn or (not mapper_confidence_ok and not permissive):
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, "rxnmapper_failed_or_low_confidence", "", "", "not_attempted", "false", "", clean_text(mapped_row.get("mapper_error"))))
                continue

            extracted = extract_rdchiral_template(clean_text(row.get("record_id")) or clean_text(row.get("source_reaction_id")), mapped_rxn)
            smarts = clean_text(extracted.get("reaction_smarts"))
            if extracted.get("extraction_status") != "ok" or not smarts:
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, clean_text(extracted.get("extraction_status")), "", "", "not_attempted", "false", "", clean_text(extracted.get("template_qc_note"))))
                continue

            ok, qc_note = reaction_smarts_valid(smarts)
            if not ok:
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, "reaction_smarts_qc_failed", smarts, sha256_text(smarts)[:20], "invalid", "false", "parse_qc_failed", qc_note))
                continue

            replay_ok, replay_note = replay_reaction_smarts_on_pair(
                smarts,
                clean_text(row.get("main_substrate_smiles")),
                clean_text(row.get("main_product_smiles")),
            )
            if not replay_ok and not permissive:
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, "reaction_smarts_replay_failed", smarts, sha256_text(smarts)[:20], "invalid", "false", replay_note, qc_note))
                continue

            successful += 1
            if permissive and not replay_ok:
                replay_failed_accepted += 1
            if permissive and not mapper_confidence_ok:
                low_confidence_accepted += 1
            template_hash = sha256_text(smarts)[:20]
            if replay_ok and mapper_confidence_ok:
                status = "ok"
            elif replay_ok and not mapper_confidence_ok:
                status = "ok_low_mapper_confidence"
            elif (not replay_ok) and mapper_confidence_ok:
                status = "reaction_smarts_replay_failed"
            else:
                status = "reaction_smarts_replay_failed_low_mapper_confidence"
            promoted.append(
                _promote_anchor_to_template(
                    row,
                    smarts,
                    confidence,
                    len(promoted) + 1,
                    release_mode=mode,
                    acceptance_status=status,
                    replay_ok=replay_ok,
                    replay_note=replay_note,
                    mapper_confidence_ok=mapper_confidence_ok,
                    qc_note=qc_note,
                )
            )
            reports.append(_report_row(row, mapped_rxn, confidence, status, smarts, template_hash, "ok", "true" if replay_ok else "false", replay_note, qc_note))

    report = pd.DataFrame(reports, columns=ANCHOR_GENERALIZATION_REPORT_COLUMNS)
    report.to_csv(release / f"anchor_generalization{suffix}.report.tsv", sep="\t", index=False)

    raw_templates = pd.DataFrame(promoted).fillna("") if promoted else _empty_template_df()
    for col in TEMPLATE_COLUMNS:
        if col not in raw_templates.columns:
            raw_templates[col] = ""
    raw_templates = raw_templates[TEMPLATE_COLUMNS]
    raw_templates.to_csv(release / f"anchor_generalization{suffix}.templates.raw.tsv", sep="\t", index=False)

    qc = qc_templates(raw_templates) if len(raw_templates) else _empty_template_df()
    dedup = deduplicate_templates(qc) if len(qc) else _empty_template_df()
    dedup.to_csv(release / f"anchor_generalization{suffix}.templates.deduplicated.tsv", sep="\t", index=False)

    rules = build_rule_library(dedup, family_evidence) if len(dedup) else _empty_rule_df()
    rules.to_csv(release / f"anchor_generalization{suffix}.rules.audit.tsv", sep="\t", index=False)

    smarts = export_reaction_smarts_library(rules, tier="all", require_qc_ok=not permissive) if len(rules) else _empty_smarts_df()
    if len(smarts):
        smarts = smarts.copy()
        smarts["smarts_library_tier"] = stem
        prefix = "SMRT_ANCHOR" if mode == "strict" else f"SMRT_ANCHOR_{mode.upper()}"
        smarts["smarts_rule_id"] = [f"{prefix}_{i+1:09d}" for i in range(len(smarts))]
    smarts.to_csv(release / f"reaction_smarts_library.{stem}.tsv", sep="\t", index=False)
    smarts.to_csv(release / f"reaction_smarts_rules.{stem}.tsv", sep="\t", index=False)

    issues, validation = validate_smarts_release(smarts)
    issues.to_csv(release / f"anchor_generalization{suffix}.validation_issues.tsv", sep="\t", index=False)
    validation = {
        **validation,
        "anchor_derived_released_smarts": int(len(smarts)),
        "strict_validation_passed": bool(len(issues) == 0),
        "permissive_release_valid": bool(
            not any(smarts.get("template_scope", pd.Series(dtype=str)).astype(str) != "generalized_template")
            and not any(smarts.get("reaction_smarts", pd.Series(dtype=str)).astype(str).str.strip() == "")
        ) if len(smarts) else True,
        "validation_passed": bool(len(issues) == 0) if mode == "strict" else True,
        "validation_issue_records": int(len(issues)),
        "release_mode": mode,
    }
    if permissive:
        validation["note"] = (
            "Permissive anchor-derived release keeps every extractable RXNMapper/RDChiral SMARTS. "
            "Replay failures and low mapper confidence are retained as explicit QC annotations and are not core-ready."
        )
    write_json(release / f"anchor_generalization{suffix}.validation.json", validation)

    summary = {
        **base_summary,
        "attempted_records": int(attempted),
        "successful_templates": int(successful),
        "failed_records": int(failed),
        "replay_failed_accepted_templates": int(replay_failed_accepted),
        "low_confidence_accepted_templates": int(low_confidence_accepted),
        "raw_anchor_derived_templates": int(len(raw_templates)),
        "deduplicated_anchor_derived_templates": int(len(dedup)),
        "anchor_derived_audit_rules": int(len(rules)),
        "anchor_derived_released_smarts": int(len(smarts)),
        "validation": validation,
        "elapsed_seconds": round(time.time() - start, 2),
        "report_path": str(release / f"anchor_generalization{suffix}.report.tsv"),
        "release_path": str(release / f"reaction_smarts_library.{stem}.tsv"),
    }
    write_json(release / f"anchor_generalization{suffix}.summary.json", summary)
    return summary
