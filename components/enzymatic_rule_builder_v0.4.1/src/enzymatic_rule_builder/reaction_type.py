from __future__ import annotations

import json
from typing import Any

from .direction import direction_mode_from_record, normalize_direction_label
from .ec import broad_ec_classes
from .utils import clean_text


def normalize_direction_label(direction: Any = "", is_reversible: Any = "") -> str:
    """Normalize source direction to forward/reverse/reversible/unknown."""
    mode = direction_mode_from_record({"direction": direction, "is_reversible": is_reversible})
    if mode == "source_reversible":
        return "reversible"
    if mode == "source_reverse":
        return "reverse"
    if mode == "source_forward":
        return "forward"
    return "unknown"


def _load_delta(row: dict[str, Any]) -> dict[str, Any]:
    text = clean_text(row.get("reaction_delta_json", ""))
    if not text:
        return {}
    try:
        data = json.loads(text)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def infer_structural_reaction_type(row: dict[str, Any]) -> tuple[str, str, str]:
    transferred = clean_text(row.get("transferred_group_class", "") or row.get("transferred_group", "")).lower()
    acceptor = clean_text(row.get("acceptor_atom_class", "")).lower()
    if "cinnamoyl" in transferred:
        subtype = transferred or clean_text(row.get("reaction_delta_fingerprint", ""))
        if acceptor.startswith("amine"):
            return "N_cinnamoylation_or_hydroxycinnamoyl_acyl_transfer", subtype, "role_aware_projection"
        if acceptor.startswith("alcohol"):
            return "O_cinnamoylation_or_cinnamate_ester_transfer", subtype, "role_aware_projection"
        return "cinnamoylation_or_aromatic_acyl_transfer", subtype, "role_aware_projection"
    if "benzoyl" in transferred or "aromatic_acyl" in transferred:
        subtype = transferred or clean_text(row.get("reaction_delta_fingerprint", ""))
        if acceptor.startswith("amine"):
            return "N_benzoylation_or_aromatic_acyl_transfer", subtype, "role_aware_projection"
        if acceptor.startswith("alcohol"):
            return "benzoylation_or_aromatic_acyl_transfer", subtype, "role_aware_projection"
        return "benzoylation_or_aromatic_acyl_transfer", subtype, "role_aware_projection"
    if "glycosyl" in transferred or "sugar" in transferred:
        subtype = transferred or clean_text(row.get("reaction_delta_fingerprint", ""))
        return "glycosylation_or_sugar_loss", subtype, "role_aware_projection"
    if "methyl" in transferred:
        subtype = transferred or clean_text(row.get("reaction_delta_fingerprint", ""))
        if acceptor.startswith("amine"):
            return "N_methylation_or_demethylation", subtype, "role_aware_projection"
        if acceptor.startswith("alcohol"):
            return "O_methylation_or_demethylation", subtype, "role_aware_projection"
        return "methylation_or_demethylation", subtype, "role_aware_projection"
    if "phosphoryl" in transferred or "phosphate" in transferred:
        subtype = transferred or clean_text(row.get("reaction_delta_fingerprint", ""))
        return "phosphorylation_or_dephosphorylation", subtype, "role_aware_projection"

    delta = _load_delta(row)
    if not delta:
        fp = clean_text(row.get("reaction_delta_fingerprint", ""))
        return ("unassigned_structural_delta", fp, "structural_delta") if fp else ("", "", "none")

    atom_delta = delta.get("atom_delta", {}) or {}
    heavy = int(delta.get("heavy_atom_delta", 0) or 0)
    mass = float(delta.get("exact_mass_delta", 0.0) or 0.0)
    ring = int(delta.get("ring_count_delta", 0) or 0)
    c = int(atom_delta.get("C", 0) or 0)
    h = int(atom_delta.get("H", 0) or 0)
    n = int(atom_delta.get("N", 0) or 0)
    o = int(atom_delta.get("O", 0) or 0)
    subtype = clean_text(row.get("reaction_delta_fingerprint", ""))

    if c == 0 and n == 0 and o > 0 and heavy <= 2:
        return "hydroxylation_or_oxygenation", subtype, "structural_delta_specific"
    if c == 0 and n == 0 and o < 0 and heavy >= -2:
        return "deoxygenation_or_dehydration", subtype, "structural_delta_specific"
    if c == 0 and o == 0 and h <= -2 and abs(mass) < 4:
        return "oxidation_or_dehydrogenation", subtype, "structural_delta_specific"
    if c == 0 and o == 0 and h >= 2 and abs(mass) < 4:
        return "reduction_or_hydrogenation", subtype, "structural_delta_specific"
    if c == 2 and o == 1 and 1 <= heavy <= 4 and 40 <= mass <= 45:
        return "acetylation_or_deacetylation_like_acyl_transfer", subtype, "structural_delta_specific"
    if c == -2 and o == -1 and -4 <= heavy <= -1 and -45 <= mass <= -40:
        return "deacetylation_or_acetyl_ester_hydrolysis", subtype, "structural_delta_specific"
    if c >= 6 and o >= 1 and n == 0 and 90 <= abs(mass) <= 130:
        return "benzoylation_or_aromatic_acyl_transfer", subtype, "structural_delta_specific"
    if c >= 8 and o >= 1 and n == 0 and 120 <= abs(mass) <= 155:
        return "cinnamoylation_or_aromatic_acyl_transfer", subtype, "structural_delta_specific"
    if c >= 5 and o >= 4 and 120 <= abs(mass) <= 190:
        return "glycosylation_or_sugar_loss", subtype, "structural_delta_specific"
    if heavy >= 8 and (n > 0 or o > 0):
        return "large_side_chain_transfer_or_loss", subtype, "structural_delta_specific"
    if ring != 0:
        return "cyclization_or_ring_rearrangement", subtype, "structural_delta_specific"
    if subtype:
        return "unassigned_structural_delta", subtype, "structural_delta"
    return "", "", "none"


def classify_reaction_type(row: dict[str, Any]) -> tuple[str, str, str]:
    source_type = clean_text(row.get("reaction_type_source", ""))
    source_subtype = clean_text(row.get("reaction_subtype_source", ""))
    if source_type or source_subtype:
        return source_type or source_subtype, source_subtype, "source"
    rt, subtype, mode = infer_structural_reaction_type(row)
    if rt:
        return rt, subtype, mode
    broad = broad_ec_classes(row.get("ec_numbers", ""), row.get("template_ec_candidates", ""), row.get("database_ec_candidates", ""))
    if broad:
        return "unassigned_ec_supported", broad, "ec_broad"
    return "unassigned", "", "none"


def direction_mode(row: dict[str, Any]) -> str:
    handling = clean_text(row.get("direction_handling", "")).lower()
    if handling.startswith("reversed_from_source"):
        return "source_reverse_corrected"
    if handling.startswith("split_reversible"):
        return "source_reversible"
    return direction_mode_from_record(row)
