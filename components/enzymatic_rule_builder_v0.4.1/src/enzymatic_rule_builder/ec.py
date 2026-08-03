from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Any

from .schema import EC_CLASS_LABELS
from .source_layers import best_layer
from .utils import clean_text, join_values, split_multi_value, unique_stable

# Supports full and partial EC numbers, e.g. 1.14.14.17, 2.3.1.-, 1.14.-.-.
EC_RE = re.compile(r"(?<![A-Za-z0-9_.-])(?:\d{1,2}|-)(?:\.(?:\d{1,3}|-)){3}(?![A-Za-z0-9_.-])")

# Evidence-type confidence represents source/evidence strength, not EC granularity.
EC_EVIDENCE_WEIGHTS = {
    "source_direct": 1.00,
    "template_direct": 0.95,
    "database_cross_reference": 0.72,
    "reaction_type_prior": 0.22,
}

EC_LAYER_MULTIPLIERS = {
    "T1_Bio_Core": 0.95,
    "T2_Bio_Extended": 0.78,
    "T3_Chem_like": 0.35,
    "Unknown": 0.50,
}

# EC granularity is deliberately not a confidence penalty. A partial EC such as
# 2.3.1.- from a curated source can be a high-confidence broad annotation. The
# granularity/specificity is reported separately in top_ec_granularity and
# top_ec_specificity.
EC_STATUS_MULTIPLIERS = {
    "full4": 1.00,
    "partial": 1.00,
    "missing": 0.00,
}

EC_SOURCE_FIELDS = [
    ("ec_numbers", "source_direct"),
    ("source_ec_numbers", "source_direct"),
    ("template_ec_candidates", "template_direct"),
    ("database_ec_candidates", "database_cross_reference"),
    ("ec_prior_candidates", "reaction_type_prior"),
]

EVIDENCE_PRIORITY = {
    "source_direct": 4,
    "template_direct": 3,
    "database_cross_reference": 2,
    "reaction_type_prior": 1,
}

EC_CLASS_BY_NUMBER = {
    "1": "EC1_oxidoreductase",
    "2": "EC2_transferase",
    "3": "EC3_hydrolase",
    "4": "EC4_lyase",
    "5": "EC5_isomerase",
    "6": "EC6_ligase",
    "7": "EC7_translocase",
}


def normalize_ec_numbers(value: Any) -> str:
    text = clean_text(value).replace("EC:", "").replace("EC ", "").replace("ec:", "").replace("ec ", "")
    text = text.replace("n", "-").replace("N", "-")
    found = EC_RE.findall(text)
    # Preserve values already split but not caught due unusual separators.
    for part in split_multi_value(text):
        if EC_RE.fullmatch(part):
            found.append(part)
    return join_values(found)


def normalize_ec_number(value: Any) -> str:
    vals = split_multi_value(normalize_ec_numbers(value))
    return vals[0] if vals else ""


def normalize_ec_prefix(value: Any) -> str:
    text = clean_text(value).replace("EC:", "").replace("EC ", "").replace("ec:", "").replace("ec ", "")
    text = text.replace("n", "-").replace("N", "-").strip().strip(";,.| ")
    norm = normalize_ec_number(text)
    if norm:
        return norm
    if re.fullmatch(r"(?:\d{1,2}|-)(?:\.(?:\d{1,3}|-)){0,3}", text):
        return text
    return ""


def ec_status(ec: str) -> str:
    ec = normalize_ec_number(ec)
    if not ec:
        return "missing"
    if "-" in ec:
        return "partial"
    return "full4"


def ec_parts(ec: str) -> list[str]:
    norm = normalize_ec_number(ec)
    if norm:
        return norm.split(".")
    text = clean_text(ec).replace("n", "-").replace("N", "-")
    # Family-evidence prefixes may be supplied as 1, 1.14, 2.3.1, or 2.3.1.-.
    if re.fullmatch(r"(?:\d{1,2}|-)(?:\.(?:\d{1,3}|-)){0,3}", text):
        return text.split(".")
    return []


def ec_specificity(ec: str) -> int:
    """Number of non-wildcard EC levels."""
    return sum(1 for p in ec_parts(ec) if p and p != "-")


def ec_granularity(ec: str) -> str:
    """Human-readable EC annotation granularity, independent of confidence."""
    spec = ec_specificity(ec)
    status = ec_status(ec)
    if status == "missing" or spec == 0:
        return "missing"
    if status == "full4":
        return "full4"
    return {
        1: "class_level",
        2: "subclass_level",
        3: "subsubclass_level",
        4: "partial4",
    }.get(spec, "partial")


def ec_prefix_match(prefix: str, ec: str) -> bool:
    """Hierarchical EC prefix match.

    This intentionally does not use string startswith, so 1.1 does not match 1.14.14.176.
    Wildcards in the prefix are treated as unconstrained suffix levels.
    """
    p_raw = ec_parts(prefix)
    e = ec_parts(ec)
    if not p_raw or not e:
        return False
    # A partial EC such as 2.3.1.- is treated as the constrained prefix 2.3.1.
    p = [x for x in p_raw if x != "-"]
    if len(p) > len(e):
        return False
    return all(pi == ei for pi, ei in zip(p, e))


def broad_ec_classes(*ec_fields: Any) -> str:
    labels = []
    for field in ec_fields:
        for ec in split_multi_value(normalize_ec_numbers(field)):
            if not ec:
                continue
            head = ec.split(".", 1)[0]
            if head in EC_CLASS_LABELS:
                labels.append(EC_CLASS_LABELS[head])
    return join_values(labels)


def _ec_class_numbers(ec_fields: list[str]) -> list[str]:
    vals: list[str] = []
    for field in ec_fields:
        for ec in split_multi_value(normalize_ec_numbers(field)):
            if not ec:
                continue
            head = ec.split(".", 1)[0]
            if head in EC_CLASS_BY_NUMBER and head not in vals:
                vals.append(head)
    return vals


def _ec_class_labels(numbers: list[str]) -> str:
    return join_values(EC_CLASS_BY_NUMBER.get(n, f"EC{n}") for n in numbers if n)


def expected_ec_classes_for_reaction_type(row_like: dict[str, Any]) -> tuple[list[str], str]:
    """Conservative reaction-type to EC-class expectation used only as QC.

    This does not predict an EC and must not overwrite source EC evidence.  It is
    a consistency screen: if a high-confidence EC label is incompatible with a
    directional reaction type, the rule should not be treated as a strict EC
    annotation until manually reviewed or independently supported.
    """
    text = " ".join(
        clean_text(row_like.get(key, ""))
        for key in [
            "reaction_type", "reaction_type_source", "reaction_subtype",
            "reaction_subtype_source", "main_functional_group_changes",
            "transferred_group", "transferred_group_class",
        ]
    ).lower()
    if not text or text.strip() in {"unassigned", "none"}:
        return [], "reaction_type_unassigned"

    # Direction-specific rules first.  For example, adding an acetyl group is
    # transferase-like, whereas removing one by ester hydrolysis is hydrolase-like.
    if "acetylation_or_deacetylation_like_acyl_transfer" in text:
        return ["2"], "directional_acyl_transfer_expected_transferase"
    if "deacetylation_or_acetyl_ester_hydrolysis" in text:
        return ["3"], "directional_deacylation_expected_hydrolase"
    if "deacetylation" in text or "deacylation" in text or "hydrolysis" in text:
        return ["3"], "hydrolytic_loss_expected_hydrolase"

    if any(token in text for token in ["hydroxylation", "oxygenation", "oxidation", "dehydrogenation", "reduction", "hydrogenation"]):
        return ["1"], "redox_or_oxygenation_expected_oxidoreductase"

    if any(token in text for token in ["cinnamoylation", "benzoylation", "acyl_transfer", "acetylation"]):
        return ["2"], "group_transfer_expected_transferase"

    if "glycosylation_or_sugar_loss" in text:
        return ["2", "3"], "glycosyl_transfer_or_glycoside_hydrolysis_expected_transferase_or_hydrolase"
    if "glycosylation" in text or "sugar_transfer" in text:
        return ["2"], "glycosyl_transfer_expected_transferase"
    if "sugar_loss" in text or "deglycosylation" in text:
        return ["3"], "glycoside_loss_expected_hydrolase"

    if "phosphorylation_or_dephosphorylation" in text:
        return ["2", "3"], "phosphoryl_transfer_or_hydrolysis_expected_transferase_or_hydrolase"
    if "dephosphorylation" in text:
        return ["3"], "phosphate_loss_expected_hydrolase"
    if "phosphorylation" in text:
        return ["2"], "phosphoryl_transfer_expected_transferase"

    if "methylation_or_demethylation" in text:
        return ["1", "2"], "methyl_transfer_or_demethylation_expected_transferase_or_oxidoreductase"
    if "methylation" in text:
        return ["2"], "methyl_transfer_expected_transferase"
    if "demethylation" in text:
        return ["1"], "demethylation_often_expected_oxidoreductase"

    if "dehydration" in text or "deoxygenation" in text:
        return ["1", "4"], "dehydration_or_deoxygenation_expected_oxidoreductase_or_lyase"
    if "cyclization" in text or "ring_rearrangement" in text:
        return ["4", "5"], "cyclization_or_rearrangement_expected_lyase_or_isomerase"
    if "large_side_chain_transfer_or_loss" in text:
        return ["2", "3"], "large_side_chain_transfer_or_loss_expected_transferase_or_hydrolase"
    return [], "no_conservative_ec_class_expectation"


def assess_ec_reaction_type_consistency(
    row_like: dict[str, Any],
    supported_ecs: list[str],
    top_ec: str = "",
) -> dict[str, str]:
    expected, mode = expected_ec_classes_for_reaction_type(row_like)
    observed = _ec_class_numbers(supported_ecs)
    top_classes = _ec_class_numbers([top_ec]) if top_ec else []
    expected_labels = _ec_class_labels(expected)
    observed_labels = _ec_class_labels(observed)
    top_label = _ec_class_labels(top_classes)

    if not supported_ecs:
        status = "missing_ec_not_assessed"
        note = "No supported EC evidence is available for reaction-type consistency QC."
    elif not expected:
        status = "not_assessed"
        note = f"No conservative EC-class expectation for this reaction type ({mode})."
    elif any(cls in expected for cls in top_classes) and all(cls in expected for cls in observed):
        status = "consistent"
        note = "Observed EC class is compatible with the directional reaction-type expectation."
    elif any(cls in expected for cls in top_classes):
        status = "mixed_supported_classes"
        note = "Top EC is compatible, but additional supported EC classes are outside the expected set."
    elif any(cls in expected for cls in observed):
        status = "mixed_top_inconsistent"
        note = "At least one candidate EC class is compatible, but the top-ranked EC class is not."
    else:
        status = "inconsistent"
        note = "Supported EC class is not compatible with the conservative directional reaction-type expectation."
    return {
        "ec_reaction_type_consistency": status,
        "ec_reaction_type_expected_classes": expected_labels,
        "ec_reaction_type_observed_classes": observed_labels,
        "ec_reaction_type_top_class": top_label,
        "ec_reaction_type_consistency_note": note,
        "ec_reaction_type_consistency_mode": mode,
    }


def ec_directionality_qc(row_like: dict[str, Any]) -> dict[str, str]:
    handling = clean_text(row_like.get("direction_handling", "")).lower()
    variant = clean_text(row_like.get("direction_variant", "")).lower()
    source_direction = clean_text(row_like.get("direction", row_like.get("source_direction", ""))).lower()
    if "unknown_direction" in handling or source_direction == "unknown":
        warning = "Source direction is unknown; EC annotation is retained only for the emitted left-to-right rule and must not be reused for the reverse edge."
    elif "split_reversible" in handling or source_direction == "reversible":
        warning = "Reversible source was split into directional rules; each direction must carry its own EC evidence and the reverse edge must not reuse this row blindly."
    elif variant == "reverse" or "reversed_from_source" in handling:
        warning = "Source orientation was corrected; EC annotation applies to the emitted substrate-to-product rule only."
    else:
        warning = ""
    return {
        "ec_directionality_scope": "direction_specific_rule",
        "ec_directionality_warning": warning,
        "reverse_ec_inheritance_policy": "do_not_reuse_for_reverse_edge;use_explicit_reverse_directional_rule_if_present",
    }


def all_candidate_ecs(*ec_fields: Any) -> str:
    return join_values(ec for field in ec_fields for ec in split_multi_value(normalize_ec_numbers(field)))


def _ec_overlaps(a: str, b: str) -> bool:
    """Whether two EC annotations are redundant parent/child labels."""
    return bool(a and b and (a == b or ec_prefix_match(a, b) or ec_prefix_match(b, a)))


def hierarchical_deduplicate_ec_records(records: list[dict[str, Any]], *, top_n: int | None = None) -> list[dict[str, Any]]:
    """De-redundify ranked EC records by parent/child overlap.

    Unlike earlier versions, this does not automatically prefer full four-level EC
    numbers. It keeps the highest-confidence representative in each overlapping
    EC group. Specificity is used only as a tie-breaker, because EC confidence is
    a property of source quality, while EC specificity is a separate annotation
    granularity field.
    """
    out: list[dict[str, Any]] = []
    for rec in records:
        ec = clean_text(rec.get("ec", ""))
        if not ec:
            continue
        if any(_ec_overlaps(ec, clean_text(existing.get("ec", ""))) for existing in out):
            continue
        out.append(rec)
        if top_n is not None and len(out) >= top_n:
            break
    return out


@dataclass
class EcEvidence:
    ec: str
    evidence_type: str
    layer: str = "Unknown"
    sources: tuple[str, ...] = field(default_factory=tuple)

    @property
    def event_weight(self) -> float:
        base = EC_EVIDENCE_WEIGHTS.get(self.evidence_type, 0.0)
        layer_mult = EC_LAYER_MULTIPLIERS.get(self.layer, EC_LAYER_MULTIPLIERS["Unknown"])
        status_mult = EC_STATUS_MULTIPLIERS.get(ec_status(self.ec), 0.0)
        return max(0.0, min(0.99, base * layer_mult * status_mult))


def _combined_confidence(weights: list[float]) -> float:
    if not weights:
        return 0.0
    miss = 1.0
    for w in weights:
        miss *= (1.0 - max(0.0, min(0.99, w)))
    return round(max(0.0, min(0.99, 1.0 - miss)), 4)


def _source_tokens(value: Any) -> tuple[str, ...]:
    return tuple(split_multi_value(value))


def _collapse_duplicate_events(events: list[EcEvidence]) -> list[EcEvidence]:
    """Avoid over-counting the same EC copied across multiple EC columns.

    If the same EC is reported by the same source/layer more than once, keep the
    strongest evidence type. Multiple independent sources are still allowed to
    increase confidence.
    """
    best: dict[tuple[str, tuple[str, ...], str], EcEvidence] = {}
    for ev in events:
        key = (ev.ec, tuple(sorted(ev.sources)), ev.layer)
        prev = best.get(key)
        if prev is None:
            best[key] = ev
            continue
        if EVIDENCE_PRIORITY.get(ev.evidence_type, 0) > EVIDENCE_PRIORITY.get(prev.evidence_type, 0):
            best[key] = ev
    return list(best.values())


def _gather_ec_evidence(row_like: dict[str, Any], layer_best: str | None = None) -> list[EcEvidence]:
    layer = clean_text(layer_best) or best_layer(row_like.get("evidence_layer", row_like.get("evidence_layers_all", "")))
    source_tokens = _source_tokens(row_like.get("source_database", "")) or _source_tokens(row_like.get("template_sources", ""))
    if not source_tokens:
        source_tokens = (clean_text(row_like.get("source_database", "Unknown")) or "Unknown",)
    events: list[EcEvidence] = []
    for field_name, evidence_type in EC_SOURCE_FIELDS:
        for ec in split_multi_value(normalize_ec_numbers(row_like.get(field_name, ""))):
            if ec:
                events.append(EcEvidence(ec=ec, evidence_type=evidence_type, layer=layer, sources=source_tokens))
    return _collapse_duplicate_events(events)


def _score_ecs(events: list[EcEvidence], *, include_priors: bool) -> list[dict[str, Any]]:
    by_ec: dict[str, list[EcEvidence]] = {}
    for ev in events:
        if not include_priors and ev.evidence_type == "reaction_type_prior":
            continue
        by_ec.setdefault(ev.ec, []).append(ev)
    scored = []
    for ec, evs in by_ec.items():
        etypes = unique_stable(ev.evidence_type for ev in evs)
        sources = unique_stable(s for ev in evs for s in ev.sources)
        layers = unique_stable(ev.layer for ev in evs)
        conf = _combined_confidence([ev.event_weight for ev in evs])
        strongest = max((ev.event_weight for ev in evs), default=0.0)
        scored.append({
            "ec": ec,
            "confidence": conf,
            "strongest_event_weight": round(strongest, 4),
            "status": ec_status(ec),
            "granularity": ec_granularity(ec),
            "specificity": ec_specificity(ec),
            "evidence_count": len(evs),
            "evidence_types": ";".join(etypes),
            "sources": ";".join(sources),
            "layers": ";".join(layers),
        })
    return scored


def _rank_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        records,
        key=lambda r: (
            float(r.get("confidence", 0.0)),
            int(r.get("evidence_count", 0)),
            int(r.get("specificity", 0)),  # specificity is a tie-breaker only
            1 if r.get("status") == "full4" else 0,
            clean_text(r.get("ec", "")),
        ),
        reverse=True,
    )


def _assignment_mode(records: list[dict[str, Any]]) -> str:
    if not records:
        return "missing"
    evidence_types = set(t for r in records for t in split_multi_value(r.get("evidence_types", "")))
    if "source_direct" in evidence_types:
        return "source_direct"
    if "template_direct" in evidence_types:
        return "template_direct"
    if "database_cross_reference" in evidence_types:
        return "database_cross_reference"
    if "reaction_type_prior" in evidence_types:
        return "prior_only"
    return "supported"


def summarize_ec_evidence(row_like: dict[str, Any], layer_best: str | None = None, *, top_n: int = 3) -> dict[str, Any]:
    """Aggregate, de-redundify, rank and summarize EC evidence for one rule.

    Confidence is evidence/source-quality confidence. EC specificity is reported
    separately as granularity and does not reduce confidence. Priors are retained
    separately and are used only as a fallback when no supported EC exists.
    """
    events = _gather_ec_evidence(row_like, layer_best)
    supported_ranked_all = _rank_records(_score_ecs(events, include_priors=False))
    prior_ranked_all = _rank_records([r for r in _score_ecs(events, include_priors=True) if "reaction_type_prior" in split_multi_value(r.get("evidence_types", ""))])

    supported_records = hierarchical_deduplicate_ec_records(supported_ranked_all)
    prior_records = hierarchical_deduplicate_ec_records(prior_ranked_all)
    ranked = supported_records if supported_records else prior_records
    assignment_mode = _assignment_mode(supported_records if supported_records else prior_records)

    all_supported_ecs = [r["ec"] for r in supported_records]
    all_supported_events_ecs = unique_stable(ev.ec for ev in events if ev.evidence_type != "reaction_type_prior")
    full_ecs = unique_stable(ec for ec in all_supported_events_ecs if ec_status(ec) == "full4")
    supported_partial_ecs = unique_stable(ec for ec in all_supported_events_ecs if ec_status(ec) == "partial")
    prior_ecs = unique_stable(ev.ec for ev in events if ev.evidence_type == "reaction_type_prior")

    broad = broad_ec_classes(";".join(all_supported_ecs or [r["ec"] for r in prior_records]))
    broad_vals = split_multi_value(broad)

    top = ranked[0] if ranked else {}
    top3 = ranked[:top_n]
    top3_ecs = [r["ec"] for r in top3]
    top3_confs = [f"{float(r['confidence']):.4f}" for r in top3]
    top3_modes = [r["evidence_types"] for r in top3]
    top3_sources = [r["sources"] for r in top3]
    top3_specificities = [str(r.get("specificity", "")) for r in top3]
    top3_granularities = [r.get("granularity", "") for r in top3]

    candidate_count = len(all_supported_ecs)
    broad_count = len(broad_vals)
    close_top2 = False
    if len(top3) > 1:
        close_top2 = abs(float(top3[0]["confidence"]) - float(top3[1]["confidence"])) < 0.10
    # Parent/child EC labels are already de-redundified and are not treated as conflicts.
    conflict_flag = bool(broad_count > 1 or (candidate_count > 1 and close_top2))
    if broad_count > 1:
        conflict_level = "high"
    elif candidate_count > 3 or close_top2:
        conflict_level = "medium"
    elif candidate_count > 1:
        conflict_level = "low"
    else:
        conflict_level = "none"

    top_conf = float(top.get("confidence", 0.0) or 0.0)
    consistency = assess_ec_reaction_type_consistency(row_like, all_supported_ecs, top.get("ec", "") if top else "")
    # Strict EC annotation now means high-confidence source/template-supported EC annotation.
    # It does not require full four-level EC precision.
    strict = bool(
        top
        and assignment_mode in {"source_direct", "template_direct"}
        and top_conf >= 0.75
        and conflict_level not in {"high"}
        and broad_count <= 1
        and consistency["ec_reaction_type_consistency"] not in {"inconsistent", "mixed_top_inconsistent"}
    )

    ec_status_summary = "missing"
    if full_ecs and supported_partial_ecs:
        ec_status_summary = "mixed"
    elif full_ecs:
        ec_status_summary = "full4"
    elif supported_partial_ecs:
        ec_status_summary = "partial"
    elif prior_ecs:
        ec_status_summary = "prior_only"

    evidence_json = json.dumps({
        "ranked_supported_ecs": supported_records,
        "ranked_prior_ecs": prior_records,
        "all_supported_ec_events": [ev.__dict__ | {"event_weight": ev.event_weight} for ev in events if ev.evidence_type != "reaction_type_prior"],
        "conflict_level": conflict_level,
        "assignment_mode": assignment_mode,
        "reaction_type_consistency": consistency,
        "confidence_definition": "source/evidence quality; EC granularity is reported separately and does not penalize confidence",
    }, ensure_ascii=False, sort_keys=True)

    top_specificity = top.get("specificity", "") if top else ""
    top_granularity = top.get("granularity", "") if top else ""

    directionality = ec_directionality_qc(row_like)
    return {
        "source_ec_numbers": all_candidate_ecs(row_like.get("ec_numbers", ""), row_like.get("source_ec_numbers", "")),
        "candidate_ec_numbers": join_values(all_supported_ecs),
        "full_ec_numbers": join_values(full_ecs),
        "partial_ec_numbers": join_values(supported_partial_ecs),
        "supported_partial_ec_numbers": join_values(supported_partial_ecs),
        "prior_ec_numbers": join_values(prior_ecs),
        "ec_prior_candidates": join_values(prior_ecs),
        "ranked_ec_numbers": join_values(r["ec"] for r in ranked),
        "top_ec_number": top.get("ec", ""),
        "top_ec_confidence": round(top_conf, 4) if top else "",
        "top_ec_assignment_mode": assignment_mode if top else "missing",
        "top_ec_evidence_types": top.get("evidence_types", ""),
        "top_ec_sources": top.get("sources", ""),
        "top_ec_specificity": top_specificity,
        "top_ec_granularity": top_granularity,
        "top3_ec_numbers": join_values(top3_ecs),
        "top3_ec_confidences": ";".join(top3_confs),
        "top3_ec_assignment_modes": join_values(top3_modes),
        "top3_ec_sources": "|".join(top3_sources),
        "top3_ec_specificities": ";".join(top3_specificities),
        "top3_ec_granularities": join_values(top3_granularities),
        "ec_status": ec_status_summary,
        "ec_annotation_scope": top_granularity,
        "ec_candidate_count": candidate_count,
        "broad_ec_classes": broad,
        "broad_ec_class_count": broad_count,
        "ec_conflict_flag": conflict_flag,
        "ec_conflict_level": conflict_level,
        **consistency,
        **directionality,
        "strict_ec_annotation_use": strict,
        "ec_evidence_summary_json": evidence_json,
    }
