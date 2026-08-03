from __future__ import annotations

from typing import Any
import re

from .utils import clean_text, split_multi_value
from .schema import LAYER_ORDER


LEGACY_LAYER_ALIASES = {
    "T0_Taxol_Curated": "T1_Bio_Core",
}


def normalize_evidence_layer(layer: Any) -> str:
    text = clean_text(layer)
    if not text:
        return "Unknown"
    return LEGACY_LAYER_ALIASES.get(text, text if text in LAYER_ORDER else "Unknown")


def normalize_evidence_layers(layers: Any) -> str:
    vals = []
    for val in split_multi_value(layers):
        norm = normalize_evidence_layer(val)
        if norm and norm not in vals:
            vals.append(norm)
    return ";".join(vals) if vals else "Unknown"


def is_chemical_like_source(*texts: Any) -> bool:
    low = " ".join(clean_text(t).lower() for t in texts if clean_text(t))
    if not low:
        return False
    # Use explicit tokens/phrases rather than broad substring matching. In particular,
    # do not let "biosynthesis" or "biosynthetic" match chemical synthesis.
    normalized = re.sub(r"[_-]+", " ", low)
    chemical_patterns = [
        r"(^|[^a-z0-9])uspto([^a-z0-9]|$)",
        r"(^|[^a-z0-9])np like([^a-z0-9]|$)",
        r"(^|[^a-z0-9])chem like([^a-z0-9]|$)",
        r"(^|[^a-z0-9])organic chemistry([^a-z0-9]|$)",
        r"(^|[^a-z0-9])organic synthesis([^a-z0-9]|$)",
        r"(^|[^a-z0-9])chemical synthesis([^a-z0-9]|$)",
        r"(^|[^a-z0-9])synthetic reaction([^a-z0-9]|$)",
        r"(^|[^a-z0-9])synthetic template([^a-z0-9]|$)",
        r"(^|[^a-z0-9])npl([^a-z0-9]|$)",
    ]
    return any(re.search(pattern, normalized) for pattern in chemical_patterns)


def infer_evidence_layer(source_database: Any = "", source_text: Any = "", fallback: str = "Unknown") -> str:
    low = " ".join([clean_text(source_database).lower(), clean_text(source_text).lower()])
    explicit = clean_text(fallback)
    if explicit and explicit != "Unknown":
        return explicit
    if "taxol" in low or "knownpathway" in low or "curated_taxol" in low:
        return "T1_Bio_Core"
    if "bionavi" in low and "biochem" in low:
        return "T2_Bio_Extended"
    if is_chemical_like_source(low):
        return "T3_Chem_like"
    if any(token in low for token in ["rhea", "retrorules", "retro_rules", "kegg", "rclass", "metanetx", "mnx"]):
        return "T1_Bio_Core"
    if "bionavi" in low:
        return "T2_Bio_Extended"
    return "Unknown"


def infer_evidence_layer_from_record(record: dict[str, Any], fallback: str = "Unknown") -> str:
    """Infer evidence layer with awareness of computability and annotation-only sources.

    This is stricter than the legacy name-only classifier. KEGG/MetaNetX cross
    references are valuable evidence, but they should not by themselves upgrade an
    annotation-only row to a direct computable biochemical template.
    """
    explicit = clean_text(fallback or record.get("evidence_layer", ""))
    if explicit and explicit != "Unknown":
        return normalize_evidence_layer(explicit)

    source_database = clean_text(record.get("source_database", ""))
    text = " ".join(
        clean_text(record.get(key, ""))
        for key in [
            "source_database", "parser_name", "source_evidence_text", "source_reaction_id",
            "reaction_type_source", "reaction_subtype_source",
        ]
    ).lower()
    has_smarts = bool(clean_text(record.get("reaction_smarts", "")))
    has_exact = bool(clean_text(record.get("reaction_smiles", "")) or (clean_text(record.get("substrate_smiles", "")) and clean_text(record.get("product_smiles", ""))))
    has_ec = bool(clean_text(record.get("ec_numbers", "")) or clean_text(record.get("template_ec_candidates", "")) or clean_text(record.get("database_ec_candidates", "")))
    has_rhea = bool(clean_text(record.get("rhea_ids", ""))) or "rhea" in text
    has_kegg = bool(clean_text(record.get("kegg_ids", ""))) or "kegg" in text or "rclass" in text
    has_mnx = bool(clean_text(record.get("metanetx_ids", ""))) or "metanetx" in text or "mnx" in text
    has_retrorules = "retrorules" in text or "retro_rules" in text
    has_bionavi_biochem = "bionavi" in text and "biochem" in text

    if "taxol" in text or "knownpathway" in text or "curated_taxol" in text:
        return "T1_Bio_Core"
    if is_chemical_like_source(source_database, text):
        return "T3_Chem_like"

    direct_bio_template = has_smarts and (has_retrorules or has_rhea or has_kegg or has_mnx)
    supported_exact_bio = has_exact and (has_rhea or has_bionavi_biochem)
    if direct_bio_template and (has_ec or has_rhea or has_mnx or has_kegg):
        return "T1_Bio_Core"
    if supported_exact_bio and (has_ec or has_rhea):
        return "T1_Bio_Core"
    if has_retrorules and has_smarts:
        return "T2_Bio_Extended"
    if has_bionavi_biochem or ("bionavi" in text and not is_chemical_like_source(text)):
        return "T2_Bio_Extended"
    if (has_kegg or has_mnx) and has_ec:
        return "T2_Bio_Extended"
    return infer_evidence_layer(source_database, text, "Unknown")


def best_layer(layers: Any) -> str:
    vals = [normalize_evidence_layer(x) for x in split_multi_value(layers)]
    vals = [x for x in vals if x and x != "Unknown"]
    if not vals:
        return "Unknown"
    return sorted(vals, key=lambda x: LAYER_ORDER.get(x, LAYER_ORDER["Unknown"]))[0]


def is_t3_only(layers: Any) -> bool:
    vals = set(normalize_evidence_layer(x) for x in split_multi_value(layers))
    return vals == {"T3_Chem_like"}
