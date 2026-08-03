from __future__ import annotations

from typing import Any

from .utils import clean_text

DEFAULT_WEIGHTS = {
    # Family annotation is deliberately excluded from rule-confidence scoring.
    # Rule validity is based on source quality, SMARTS QC/replay, reaction/EC
    # evidence, direction and scope. Family evidence is reported separately as
    # downstream genome-mining support.
    "components": {"source": 0.35, "qc": 0.30, "family": 0.00, "reaction_type": 0.20, "direction": 0.05, "scope": 0.10},
    "source_score": {"T1_Bio_Core": 0.95, "T2_Bio_Extended": 0.65, "T3_Chem_like": 0.25, "Unknown": 0.1},
    "qc_score": {"ok": 1.0, "exact_pair_ok": 0.70, "rdkit_unavailable": 0.35, "no_template": 0.0, "invalid": 0.0},
    "family_score": {"external_evidence": 1.0, "ec_supported_family_unassigned": 0.55, "none": 0.0},
    "reaction_type_score": {"source": 1.0, "role_aware_projection": 0.85, "structural_delta_specific": 0.80, "structural_delta": 0.60, "ec_broad": 0.55, "none": 0.0},
    "direction_score": {
        "direction_qc_ok": 1.0,
        "direction_qc_reversed_to_substrate_product": 1.0,
        "direction_qc_reversible_split": 0.85,
        "direction_qc_unknown_exploratory_only": 0.10,
        "source_forward": 1.0,
        "source_reverse_corrected": 1.0,
        "source_reversible": 0.8,
        "source_reverse": 0.6,
        "unknown": 0.10,
    },
    "scope_score": {"generalized_template": 1.0, "exact_anchor": 0.4, "no_template": 0.0},
}


def score_lookup(table_name: str, key: Any, weights: dict | None = None, default: float = 0.0) -> float:
    cfg = weights or DEFAULT_WEIGHTS
    table = cfg.get(table_name, DEFAULT_WEIGHTS.get(table_name, {}))
    return float(table.get(clean_text(key), default))


def combine_scores(parts: dict[str, float], weights: dict | None = None) -> float:
    cfg = weights or DEFAULT_WEIGHTS
    comp = cfg.get("components", DEFAULT_WEIGHTS["components"])
    total = 0.0
    denom = 0.0
    for name, weight in comp.items():
        w = float(weight)
        total += w * float(parts.get(name, 0.0))
        denom += w
    return round(total / denom, 4) if denom else 0.0
