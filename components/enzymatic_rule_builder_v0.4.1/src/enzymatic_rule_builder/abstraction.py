from __future__ import annotations

import time
import re
import signal
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd

from .chem import reaction_smarts_valid, replay_reaction_smarts_on_pair
from .utils import clean_text, ensure_dir, join_values, write_json


ABSTRACTION_REPORT_COLUMNS = [
    "record_id", "source_database", "source_reaction_id", "evidence_layer",
    "canonical_reaction_smiles", "main_substrate_smiles", "main_product_smiles",
    "mapped_reaction_smiles", "rxnmapper_confidence", "extraction_status",
    "reaction_smarts", "template_qc_status", "smarts_replay_pass", "smarts_replay_note",
    "template_qc_note",
]

RXNMAPPER_TOKEN_RE = re.compile(
    r"\[[^\]]+\]|Br|Cl|Si|Se|Li|Na|Mg|Al|Ca|Fe|Zn|@@?|>>|[A-Za-z]|\d+|."
)


def estimate_rxnmapper_token_count(reaction_smiles: str) -> int:
    """Conservative token-count estimate for RXNMapper transformer inputs."""
    rxn = clean_text(reaction_smiles)
    return len(RXNMAPPER_TOKEN_RE.findall(rxn)) if rxn else 0


class MapperTimeoutError(TimeoutError):
    pass


@contextmanager
def mapper_time_limit(seconds: int | None):
    """Unix alarm-based timeout around RXNMapper calls."""
    if not seconds or seconds <= 0:
        yield
        return

    def _handle_timeout(signum, frame):  # type: ignore[no-untyped-def]
        raise MapperTimeoutError(f"rxnmapper_timeout_seconds={seconds}")

    old_handler = signal.signal(signal.SIGALRM, _handle_timeout)
    old_alarm = signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(old_alarm)
        signal.signal(signal.SIGALRM, old_handler)


def extract_rdchiral_template(source_id: str, mapped_reaction_smiles: str) -> dict[str, str]:
    """Extract a forward generalized reaction SMARTS from an atom-mapped reaction."""
    try:
        from rdchiral.template_extractor import extract_from_reaction
    except Exception as exc:  # pragma: no cover - depends on optional environment
        return {"extraction_status": "rdchiral_unavailable", "template_qc_note": str(exc)}

    rxn = clean_text(mapped_reaction_smiles)
    if ">>" not in rxn:
        return {"extraction_status": "mapped_reaction_missing_arrow", "template_qc_note": ""}

    reactants, products = rxn.split(">>", 1)
    try:
        extracted = extract_from_reaction({"_id": source_id, "reactants": reactants, "products": products})
    except Exception as exc:
        return {"extraction_status": "rdchiral_extract_failed", "template_qc_note": str(exc)}

    reactant_template = clean_text(extracted.get("reactants") if extracted else "")
    product_template = clean_text(extracted.get("products") if extracted else "")
    if not reactant_template or not product_template:
        return {"extraction_status": "rdchiral_no_template", "template_qc_note": ""}

    # RDChiral reports retrosynthetic SMARTS. For downstream molecular network construction
    # we need the forward substrate-to-product direction.
    return {
        "extraction_status": "ok",
        "reaction_smarts": reactant_template + ">>" + product_template,
        "rdchiral_retro_smarts": clean_text(extracted.get("reaction_smarts")),
        "template_qc_note": "",
    }


def map_reaction_smiles_batch(mapper: Any, reaction_smiles: list[str], *, timeout_seconds: int | None = 120) -> list[dict[str, Any]]:
    """Map a batch with RXNMapper, falling back to per-reaction mapping on batch failure."""
    try:
        with mapper_time_limit(timeout_seconds):
            mapped = mapper.get_attention_guided_atom_maps(reaction_smiles)
        if isinstance(mapped, list):
            return mapped
        return [
            {"mapped_rxn": "", "confidence": 0.0, "mapper_error": f"unexpected_mapper_result_type:{type(mapped).__name__}"}
            for _ in reaction_smiles
        ]
    except Exception as batch_exc:
        out: list[dict[str, Any]] = []
        for rxn in reaction_smiles:
            try:
                with mapper_time_limit(timeout_seconds):
                    result = mapper.get_attention_guided_atom_maps([rxn])
                out.append(result[0] if isinstance(result, list) and result else {"mapped_rxn": "", "confidence": 0.0, "mapper_error": "single_empty_result"})
            except Exception as exc:
                out.append({"mapped_rxn": "", "confidence": 0.0, "mapper_error": f"{type(batch_exc).__name__}: {batch_exc}; single: {exc}"})
        return out


def _candidate_mask(main_df: pd.DataFrame, enabled_layers: set[str]) -> pd.Series:
    if main_df.empty:
        return pd.Series([], dtype=bool)
    has_rxn = main_df.get("canonical_reaction_smiles", "").astype(str).str.contains(">>", regex=False)
    no_smarts = main_df.get("reaction_smarts", "").astype(str).str.strip().eq("")
    layer_ok = main_df.get("evidence_layer", "").astype(str).isin(enabled_layers)
    return has_rxn & no_smarts & layer_ok


def abstract_exact_reactions_to_smarts(
    main_df: pd.DataFrame,
    output_dir: str | Path,
    *,
    enabled_layers: list[str] | None = None,
    batch_size: int = 16,
    min_mapper_confidence: float = 0.50,
    max_records: int | None = None,
    max_reaction_chars: int = 2000,
    max_mapper_tokens: int = 500,
    mapper_timeout_seconds: int = 120,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Convert selected exact reactions to generalized reaction SMARTS.

    The function deliberately operates only on selected evidence layers. By default it
    includes T0 curated exact anchors because curated pathway reactions are high-quality
    data sources. T0-derived SMARTS are flagged for benchmark exclusion so independent
    external-recall tests do not leak known-pathway information.
    """
    layers = set(enabled_layers or ["T1_Bio_Core", "T2_Bio_Extended"])
    out_dir = ensure_dir(output_dir)
    report_path = out_dir / "exact_reaction_abstraction.tsv"
    summary_path = out_dir / "exact_reaction_abstraction.summary.json"

    df = main_df.copy().fillna("")
    cand_mask = _candidate_mask(df, layers)
    candidate_idx = list(df.index[cand_mask])
    if max_records is not None:
        candidate_idx = candidate_idx[: int(max_records)]

    start = time.time()
    if not candidate_idx:
        report = pd.DataFrame(columns=ABSTRACTION_REPORT_COLUMNS)
        summary = {
            "candidate_records": 0,
            "attempted_records": 0,
            "successful_templates": 0,
            "failed_records": 0,
            "enabled_layers": sorted(layers),
            "elapsed_seconds": 0.0,
        }
        report.to_csv(report_path, sep="\t", index=False)
        write_json(summary_path, summary)
        return df, report, summary

    try:
        from rxnmapper import RXNMapper
    except Exception as exc:  # pragma: no cover - depends on optional environment
        reports = []
        for idx in candidate_idx:
            row = df.loc[idx].to_dict()
            reports.append(_report_row(row, "", 0.0, "rxnmapper_unavailable", "", "not_attempted", "false", "rxnmapper_unavailable", str(exc)))
        report = pd.DataFrame(reports, columns=ABSTRACTION_REPORT_COLUMNS)
        report.to_csv(report_path, sep="\t", index=False)
        summary = {
            "candidate_records": int(len(candidate_idx)),
            "attempted_records": 0,
            "successful_templates": 0,
            "failed_records": int(len(candidate_idx)),
            "enabled_layers": sorted(layers),
            "elapsed_seconds": round(time.time() - start, 3),
            "optional_dependency_status": "rxnmapper_unavailable",
            "optional_dependency_error": str(exc),
        }
        write_json(summary_path, summary)
        return df, report, summary

    mapper = RXNMapper()
    reports: list[dict[str, Any]] = []
    attempted = successful = failed = 0

    for offset in range(0, len(candidate_idx), max(1, batch_size)):
        batch_idx = candidate_idx[offset: offset + max(1, batch_size)]
        if offset == 0 or offset % max(1, batch_size * 100) == 0:
            print(f"[exact_abstraction] processed_candidates={offset}/{len(candidate_idx)}", flush=True)
        batch_rows = [df.loc[i].to_dict() for i in batch_idx]
        batch_rxns = []
        valid_positions = []
        for pos, row in enumerate(batch_rows):
            rxn = clean_text(row.get("canonical_reaction_smiles"))
            token_count = estimate_rxnmapper_token_count(rxn)
            too_long_chars = bool(max_reaction_chars and len(rxn) > max_reaction_chars)
            too_many_tokens = bool(max_mapper_tokens and token_count > max_mapper_tokens)
            if not rxn or ">>" not in rxn or too_long_chars or too_many_tokens:
                failed += 1
                note = join_values([
                    "missing_or_invalid_reaction_smiles" if (not rxn or ">>" not in rxn) else "",
                    f"reaction_chars={len(rxn)}",
                    f"max_reaction_chars={max_reaction_chars}",
                    f"estimated_rxnmapper_tokens={token_count}",
                    f"max_mapper_tokens={max_mapper_tokens}",
                    "skipped_before_rxnmapper_to_avoid_transformer_length_failure" if too_long_chars or too_many_tokens else "",
                ])
                reports.append(_report_row(row, "", 0.0, "invalid_or_too_long_exact_reaction_smiles", "", "not_attempted", "false", "not_attempted", note))
                continue
            batch_rxns.append(rxn)
            valid_positions.append(pos)
            attempted += 1

        mapped_results = map_reaction_smiles_batch(mapper, batch_rxns, timeout_seconds=mapper_timeout_seconds) if batch_rxns else []
        for pos, mapped_row in zip(valid_positions, mapped_results):
            df_index = batch_idx[pos]
            row = batch_rows[pos]
            mapped_rxn = clean_text(mapped_row.get("mapped_rxn"))
            confidence = float(mapped_row.get("confidence") or 0.0)
            if not mapped_rxn or confidence < min_mapper_confidence:
                failed += 1
                note = clean_text(mapped_row.get("mapper_error"))
                reports.append(_report_row(row, mapped_rxn, confidence, "rxnmapper_failed_or_low_confidence", "", "not_attempted", "false", "not_attempted", note))
                continue

            extracted = extract_rdchiral_template(clean_text(row.get("record_id")) or clean_text(row.get("source_reaction_id")), mapped_rxn)
            smarts = clean_text(extracted.get("reaction_smarts"))
            if extracted.get("extraction_status") != "ok" or not smarts:
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, clean_text(extracted.get("extraction_status")), "", "not_attempted", "false", "not_attempted", clean_text(extracted.get("template_qc_note"))))
                continue

            ok, qc_note = reaction_smarts_valid(smarts)
            if not ok:
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, "reaction_smarts_qc_failed", smarts, "invalid", "false", "parse_qc_failed", qc_note))
                continue

            replay_ok, replay_note = replay_reaction_smarts_on_pair(
                smarts,
                clean_text(row.get("main_substrate_smiles")),
                clean_text(row.get("main_product_smiles")),
            )
            if not replay_ok:
                failed += 1
                reports.append(_report_row(row, mapped_rxn, confidence, "reaction_smarts_replay_failed", smarts, "invalid", "false", replay_note, qc_note))
                continue

            successful += 1
            df.at[df_index, "reaction_smarts"] = smarts
            df.at[df_index, "abstracted_from_exact_reaction"] = "true"
            curated_anchor = clean_text(row.get("curated_taxol_anchor", "")).lower() in {"true", "1", "yes"}
            df.at[df_index, "derived_from_exact_anchor"] = "true" if curated_anchor else "false"
            df.at[df_index, "rxnmapper_confidence"] = f"{confidence:.4f}"
            df.at[df_index, "rdchiral_extraction_status"] = "ok"
            df.at[df_index, "abstracted_smarts_applies_to_original_pair"] = "true"
            df.at[df_index, "exact_abstraction_qc_status"] = "pass"
            df.at[df_index, "benchmark_exclusion_flag"] = "curated_taxol_anchor_derived_exclude_from_external_recall" if curated_anchor else "none"
            df.at[df_index, "source_evidence_text"] = join_values([
                row.get("source_evidence_text", ""),
                "exact_reaction_generalized_by_rxnmapper_rdchiral",
                "smarts_replay_validated_against_original_pair",
                f"rxnmapper_confidence={confidence:.4f}",
            ])
            reports.append(_report_row(row, mapped_rxn, confidence, "ok", smarts, "ok", "true", replay_note, qc_note))

    report = pd.DataFrame(reports, columns=ABSTRACTION_REPORT_COLUMNS)
    report.to_csv(report_path, sep="\t", index=False)
    summary = {
        "candidate_records": int(cand_mask.sum()),
        "attempted_records": int(attempted),
        "successful_templates": int(successful),
        "failed_records": int(failed),
        "enabled_layers": sorted(layers),
        "min_mapper_confidence": float(min_mapper_confidence),
        "max_reaction_chars": int(max_reaction_chars) if max_reaction_chars else None,
        "max_mapper_tokens": int(max_mapper_tokens) if max_mapper_tokens else None,
        "mapper_timeout_seconds": int(mapper_timeout_seconds) if mapper_timeout_seconds else None,
        "elapsed_seconds": round(time.time() - start, 2),
        "report_path": str(report_path),
    }
    write_json(summary_path, summary)
    return df, report, summary


def _report_row(
    row: dict[str, Any],
    mapped_rxn: str,
    confidence: float,
    status: str,
    smarts: str,
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
        "template_qc_status": clean_text(qc_status),
        "smarts_replay_pass": clean_text(replay_pass),
        "smarts_replay_note": clean_text(replay_note),
        "template_qc_note": clean_text(qc_note),
    }
