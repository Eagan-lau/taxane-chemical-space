from __future__ import annotations

import json
import re
from typing import Any

from .utils import clean_text, join_values, split_multi_value

KNOWN_COMPOSITE_TOKENS = [
    "net reaction",
    "overall reaction",
    "overall transformation",
    "merged reaction",
    "composite reaction",
    "multi-step",
    "multistep",
    "two-step",
    "pathway-level",
    "pathway level",
    "cascade reaction",
    "coupled reaction",
    "multiple enzymes",
    "multi enzyme",
    "multi-enzyme",
]

CHANGE_KEYWORDS = [
    ("glycosyl", "glycosylation_or_deglycosylation"),
    ("sugar", "glycosylation_or_deglycosylation"),
    ("acetyl", "acetylation_or_deacetylation"),
    ("deacetyl", "acetylation_or_deacetylation"),
    ("benzoyl", "benzoylation_or_debenzoylation"),
    ("debenzoyl", "benzoylation_or_debenzoylation"),
    ("acyl", "acyl_transfer_or_acyl_loss"),
    ("methyl", "methylation_or_demethylation"),
    ("hydroxyl", "hydroxylation_or_dehydroxylation"),
    ("oxygenation", "oxygenation_or_deoxygenation"),
    ("deoxygen", "oxygenation_or_deoxygenation"),
    ("oxidation", "oxidation_or_reduction"),
    ("reduction", "oxidation_or_reduction"),
    ("dehydrogenation", "oxidation_or_reduction"),
    ("hydrogenation", "oxidation_or_reduction"),
    ("ring", "ring_rearrangement_or_cyclization"),
    ("cyclization", "ring_rearrangement_or_cyclization"),
    ("rearrangement", "ring_rearrangement_or_cyclization"),
    ("side_chain", "side_chain_transfer_or_loss"),
    ("side-chain", "side_chain_transfer_or_loss"),
]


def _load_delta(row: dict[str, Any]) -> dict[str, Any]:
    text = clean_text(row.get("reaction_delta_json", ""))
    if not text:
        return {}
    try:
        data = json.loads(text)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _known_composite_by_source(row: dict[str, Any]) -> bool:
    text = " ".join(
        clean_text(row.get(k, "")).lower()
        for k in [
            "source_evidence_text",
            "reaction_type_source",
            "reaction_subtype_source",
            "reaction_equation",
            "template_qc_note",
            "notes",
        ]
    )
    return any(tok in text for tok in KNOWN_COMPOSITE_TOKENS)


def _changes_from_text(*texts: Any) -> list[str]:
    joined = " ".join(clean_text(t).lower() for t in texts if clean_text(t))
    # Split explicit compound labels like hydroxylation+acetylation or oxidation;acylation.
    explicit_parts = re.split(r"\s*(?:\+|;|/|,|\band\b)\s*", joined)
    changes: list[str] = []
    for key, label in CHANGE_KEYWORDS:
        if key in joined:
            changes.append(label)
    # If explicit text contains multiple recognizable transformation words, keep all unique labels.
    # Fall through to a stable de-duplication below.
    seen = set()
    out = []
    for c in changes:
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def _changes_from_delta(row: dict[str, Any], existing: list[str]) -> list[str]:
    delta = _load_delta(row)
    if not delta:
        return existing
    atom_delta = delta.get("atom_delta", {}) or {}
    heavy = int(delta.get("heavy_atom_delta", 0) or 0)
    mass = float(delta.get("exact_mass_delta", 0.0) or 0.0)
    ring = int(delta.get("ring_count_delta", 0) or 0)
    c = int(atom_delta.get("C", 0) or 0)
    h = int(atom_delta.get("H", 0) or 0)
    n = int(atom_delta.get("N", 0) or 0)
    o = int(atom_delta.get("O", 0) or 0)
    changes = list(existing)
    if c == 0 and n == 0 and o > 0 and heavy <= 2:
        changes.append("oxygenation_or_hydroxylation")
    if c == 0 and n == 0 and o < 0 and heavy >= -2:
        changes.append("deoxygenation_or_dehydration")
    if c == 0 and o == 0 and h <= -2 and abs(mass) < 4:
        changes.append("oxidation_or_dehydrogenation")
    if c == 0 and o == 0 and h >= 2 and abs(mass) < 4:
        changes.append("reduction_or_hydrogenation")
    if abs(c) == 2 and abs(o) == 1 and 35 <= abs(mass) <= 50:
        changes.append("acetylation_or_deacetylation")
    if abs(c) >= 5 and abs(o) >= 4 and 120 <= abs(mass) <= 210:
        changes.append("glycosylation_or_deglycosylation")
    if abs(c) >= 6 and abs(o) >= 1 and 80 <= abs(mass) <= 160:
        changes.append("aromatic_acyl_transfer_or_loss")
    if ring != 0:
        changes.append("ring_rearrangement_or_cyclization")
    if abs(heavy) >= 10 and not changes:
        changes.append("large_group_transfer_or_loss")
    seen = set()
    out = []
    for c0 in changes:
        if c0 and c0 not in seen:
            seen.add(c0)
            out.append(c0)
    return out


def annotate_biochemical_step_granularity(row: dict[str, Any], reaction_type: str = "", reaction_subtype: str = "") -> dict[str, str]:
    """Annotate the biochemical granularity of a single SMARTS rule.

    This annotation deliberately does not claim ground-truth enzyme mechanism. It
    distinguishes a single SMARTS application from the biological interpretation
    of that transform. The output is intended to be carried into downstream edge
    tables as explanatory metadata.
    """
    scope = clean_text(row.get("template_scope", ""))
    smarts = clean_text(row.get("reaction_smarts", ""))
    delta_available = bool(clean_text(row.get("reaction_delta_json", "")) or clean_text(row.get("reaction_delta_fingerprint", "")) or clean_text(reaction_type))
    if scope in {"exact_anchor", "no_template"} or (scope and scope != "generalized_template") or (not smarts and not delta_available):
        return {
            "rule_application_unit": "not_applicable_exact_anchor",
            "biochemical_step_granularity": "uncertain",
            "biochemical_step_granularity_confidence": "0.00",
            "granularity_assignment_mode": "not_a_predictive_smarts_rule",
            "granularity_evidence_summary": "exact anchors are not SMARTS rules",
            "composite_rule_flag": "false",
            "reaction_center_count": "",
            "independent_reaction_center_count": "",
            "functional_group_change_count": "",
            "main_functional_group_changes": "",
        }

    if _known_composite_by_source(row):
        changes = _changes_from_text(reaction_type, reaction_subtype, row.get("reaction_delta_fingerprint", ""))
        changes = _changes_from_delta(row, changes)
        count = max(1, len(changes)) if changes else ""
        return {
            "rule_application_unit": "single_smarts_application",
            "biochemical_step_granularity": "known_composite_step",
            "biochemical_step_granularity_confidence": "0.95",
            "granularity_assignment_mode": "source_metadata_known_composite",
            "granularity_evidence_summary": "source metadata indicates net/overall/multi-step/composite transformation",
            "composite_rule_flag": "true",
            "reaction_center_count": str(count),
            "independent_reaction_center_count": str(count),
            "functional_group_change_count": str(count),
            "main_functional_group_changes": join_values(changes),
        }

    changes = _changes_from_text(reaction_type, reaction_subtype, row.get("reaction_delta_fingerprint", ""))
    changes = _changes_from_delta(row, changes)
    # If source explicitly encodes multiple reaction-type tokens in the field, treat it as composite suspicion.
    source_text = clean_text(row.get("reaction_type_source", "") or reaction_type)
    explicit_multi = bool(re.search(r"\b(?:and|plus)\b|[+;]", source_text.lower()))
    change_count = max(1, len(changes)) if changes else 0

    if change_count >= 2 or explicit_multi:
        label = "possible_composite_step"
        conf = 0.78 if change_count >= 2 else 0.68
        mode = "multi_functional_group_delta" if change_count >= 2 else "multi_token_reaction_type"
        composite = "true"
    elif change_count == 1:
        label = "likely_single_enzyme_step"
        conf = 0.86
        mode = "single_functional_group_delta"
        composite = "false"
    else:
        # A computable SMARTS but no representative delta can still be used for network construction.
        label = "uncertain"
        conf = 0.40
        mode = "insufficient_structural_delta"
        composite = "false"

    center_count = str(change_count if change_count else "")
    return {
        "rule_application_unit": "single_smarts_application",
        "biochemical_step_granularity": label,
        "biochemical_step_granularity_confidence": f"{conf:.2f}",
        "granularity_assignment_mode": mode,
        "granularity_evidence_summary": "source metadata + reaction-type text + functional-group delta heuristic",
        "composite_rule_flag": composite,
        "reaction_center_count": center_count,
        "independent_reaction_center_count": center_count,
        "functional_group_change_count": str(change_count if change_count else ""),
        "main_functional_group_changes": join_values(changes),
    }
