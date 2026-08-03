from __future__ import annotations

from typing import Any

import pandas as pd

from .ec import ec_prefix_match, normalize_ec_number, normalize_ec_prefix
from .utils import clean_text, join_values, split_multi_value


def _ec_prefixes(ec: str) -> list[str]:
    ec = normalize_ec_number(ec)
    if not ec:
        return []
    parts = ec.split(".")
    vals = []
    for i in range(len(parts), 0, -1):
        prefix = ".".join(parts[:i])
        # Convert trailing wildcard form such as 2.3.1.- to its constrained prefix 2.3.1
        vals.append(prefix.rstrip(".-"))
    return [v for v in vals if v]


def _row_confidence(row: pd.Series) -> float:
    try:
        return float(row.get("confidence", 1.0) or 1.0)
    except Exception:
        return 1.0


def _ids(rule_like: dict[str, Any], *keys: str) -> set[str]:
    vals: list[str] = []
    for key in keys:
        vals.extend(split_multi_value(rule_like.get(key, "")))
    return set(vals)


def annotate_families(rule_like: dict[str, Any], family_evidence: pd.DataFrame | None) -> tuple[str, str, str, str]:
    """Assign enzyme families only from an external evidence table.

    Supported evidence modes include source reaction ID, Rhea ID, KEGG ID, MetaNetX ID, protein ID, EC exact,
    and EC prefix. There is intentionally no reaction-type-to-family hard-coded mapping here.
    """
    if family_evidence is None or family_evidence.empty:
        return "", "", "none", ""
    fe = family_evidence.fillna("").copy()
    hits: list[tuple[float, str, str, str]] = []
    source_ids = _ids(rule_like, "source_reaction_id", "source_reaction_ids")
    rhea_ids = _ids(rule_like, "rhea_ids")
    kegg_ids = _ids(rule_like, "kegg_ids")
    metanetx_ids = _ids(rule_like, "metanetx_ids")
    protein_ids = _ids(rule_like, "protein_ids")
    source_dbs = _ids(rule_like, "source_database", "template_sources")
    ecs_raw = split_multi_value(rule_like.get("ec_numbers", "")) + split_multi_value(rule_like.get("source_ec_numbers", "")) + split_multi_value(rule_like.get("candidate_ec_numbers", "")) + split_multi_value(rule_like.get("template_ec_candidates", ""))
    ecs = [normalize_ec_number(ec) for ec in ecs_raw if normalize_ec_number(ec)]
    ec_prefixes = set(p for ec in ecs for p in _ec_prefixes(ec))

    for _, row in fe.iterrows():
        primary = clean_text(row.get("primary_family", "")) or clean_text(row.get("family", ""))
        secondary = clean_text(row.get("secondary_families", ""))
        if not primary and not secondary:
            continue
        matched = False
        reason: list[str] = []
        mt = clean_text(row.get("match_type", "")).lower()
        mv = clean_text(row.get("match_value", ""))
        row_db = clean_text(row.get("source_database", ""))

        def mark(label: str) -> None:
            nonlocal matched
            matched = True
            reason.append(label)

        if mt == "source_reaction_id" and mv in source_ids:
            mark(f"source_reaction_id:{mv}")
        if mt == "source_database_reaction_id" and mv in source_ids and (not row_db or row_db in source_dbs):
            mark(f"source_database_reaction_id:{row_db}:{mv}")
        if mt == "rhea_id" and mv in rhea_ids:
            mark(f"rhea_id:{mv}")
        if mt == "kegg_id" and mv in kegg_ids:
            mark(f"kegg_id:{mv}")
        if mt in {"metanetx_id", "mnx_id"} and mv in metanetx_ids:
            mark(f"metanetx_id:{mv}")
        if mt == "protein_id" and mv in protein_ids:
            mark(f"protein_id:{mv}")
        mv_ec = normalize_ec_number(mv)
        if mt in {"ec", "ec_exact"} and mv_ec in ecs:
            mark(f"ec_exact:{mv_ec}")
        mv_prefix = normalize_ec_prefix(mv)
        if mt == "ec_prefix" and mv_prefix and any(ec_prefix_match(mv_prefix, ec) for ec in ecs):
            mark(f"ec_prefix:{mv_prefix}")

        row_sid = clean_text(row.get("source_reaction_id", ""))
        row_rhea = clean_text(row.get("rhea_id", ""))
        row_kegg = clean_text(row.get("kegg_id", ""))
        row_mnx = clean_text(row.get("metanetx_id", ""))
        row_ec = clean_text(row.get("ec_number", ""))
        row_protein = clean_text(row.get("protein_id", ""))
        if row_sid and row_sid in source_ids and (not row_db or not source_dbs or row_db in source_dbs):
            mark(f"source_reaction_id:{row_sid}")
        if row_rhea and row_rhea in rhea_ids:
            mark(f"rhea_id:{row_rhea}")
        if row_kegg and row_kegg in kegg_ids:
            mark(f"kegg_id:{row_kegg}")
        if row_mnx and row_mnx in metanetx_ids:
            mark(f"metanetx_id:{row_mnx}")
        if row_protein and row_protein in protein_ids:
            mark(f"protein_id:{row_protein}")
        row_ec_norm = normalize_ec_number(row_ec)
        if row_ec_norm and (row_ec_norm in ecs or any(ec_prefix_match(row_ec_norm, ec) for ec in ecs)):
            mark(f"ec:{row_ec_norm}")
        if matched:
            conf = _row_confidence(row)
            ev = clean_text(row.get("evidence_source", "")) or clean_text(row.get("provenance", "")) or "family_evidence_table"
            hits.append((conf, primary, secondary, ev + "[" + ",".join(sorted(set(reason))) + "]"))
    if not hits:
        return "", "", "none", ""
    hits.sort(key=lambda x: x[0], reverse=True)
    primaries = join_values(h[1] for h in hits if h[1])
    secondaries = join_values(h[2] for h in hits if h[2])
    evidence = join_values(h[3] for h in hits)
    return primaries, secondaries, "external_evidence", evidence
