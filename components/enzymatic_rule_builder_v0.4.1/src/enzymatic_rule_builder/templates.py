from __future__ import annotations

import pandas as pd

from .chem import build_reaction_smiles, canonical_reaction_smiles, reaction_delta, reaction_smarts_valid, reverse_reaction_transform
from .direction import direction_qc_from_handling, source_direction_mode
from .schema import TEMPLATE_COLUMNS
from .utils import clean_text, join_values, sha256_text


def _smarts_variants(smarts: str, direction: str, is_reversible: str) -> list[dict[str, str]]:
    """Return directional SMARTS variants for a single source template.

    Each returned SMARTS is intended to be applied left-to-right by downstream
    graph construction. `reverse` is treated as a direction-correction request;
    `reversible` is split into two separate directional rules.
    """
    s = clean_text(smarts)
    if not s:
        return []
    mode = source_direction_mode(direction, is_reversible)
    if mode == "source_reverse":
        rev = reverse_reaction_transform(s)
        return [{
            "reaction_smarts": rev or s,
            "direction_handling": "reversed_from_source",
            "direction_variant": "forward_after_source_reverse_correction",
            "reverse_template_hash": "",
        }]
    if mode == "source_reversible":
        rev = reverse_reaction_transform(s)
        if rev and rev != s:
            forward_hash = sha256_text(s)[:20]
            reverse_hash = sha256_text(rev)[:20]
            return [
                {
                    "reaction_smarts": s,
                    "direction_handling": "split_reversible_forward",
                    "direction_variant": "forward",
                    "reverse_template_hash": reverse_hash,
                },
                {
                    "reaction_smarts": rev,
                    "direction_handling": "split_reversible_reverse",
                    "direction_variant": "reverse",
                    "reverse_template_hash": forward_hash,
                },
            ]
        return [{
            "reaction_smarts": s,
            "direction_handling": "split_reversible_single_symmetric_or_unparsed",
            "direction_variant": "forward",
            "reverse_template_hash": "",
        }]
    if mode == "source_forward":
        return [{
            "reaction_smarts": s,
            "direction_handling": "kept_forward",
            "direction_variant": "forward",
            "reverse_template_hash": "",
        }]
    return [{
        "reaction_smarts": s,
        "direction_handling": "unknown_direction_kept_left_to_right",
        "direction_variant": "left_to_right",
        "reverse_template_hash": "",
    }]


def _set_direction_qc(out: dict) -> None:
    handling = out.get("direction_handling", "")
    variant = out.get("direction_variant", "")
    qc, note = direction_qc_from_handling(handling, variant)
    out["normalized_direction"] = "substrate_to_product"
    out["direction_qc_status"] = qc
    out["direction_qc_note"] = note


def _reverse_example_pair(out: dict) -> None:
    sub = clean_text(out.get("main_substrate_smiles", ""))
    prod = clean_text(out.get("main_product_smiles", ""))
    if not sub or not prod:
        return
    out["main_substrate_smiles"], out["main_product_smiles"] = prod, sub
    out["canonical_substrate_smiles"], out["canonical_product_smiles"] = prod, sub
    out["canonical_reaction_smiles"] = build_reaction_smiles([prod], [sub])
    fp, js = reaction_delta(prod, sub)
    out["reaction_delta_fingerprint"] = fp
    out["reaction_delta_json"] = js


def build_templates(main_pairs: pd.DataFrame, *, allow_exact_pairs_as_predictive_rules: bool = False) -> pd.DataFrame:
    rows = []
    for i, row in main_pairs.fillna("").iterrows():
        base = row.to_dict()
        smarts = clean_text(base.get("reaction_smarts", ""))
        can_rxn = clean_text(base.get("canonical_reaction_smiles", "")) or canonical_reaction_smiles(base.get("reaction_smiles", ""))
        source_direction = clean_text(base.get("direction", ""))
        source_reversible = clean_text(base.get("is_reversible", ""))
        variants = _smarts_variants(smarts, source_direction, source_reversible) if smarts else []
        if variants:
            group_id = clean_text(base.get("reversible_group_id", ""))
            if len(variants) > 1 and not group_id:
                group_id = sha256_text(f"DIRGROUP::{smarts}::{base.get('record_id','')}::{base.get('source_reaction_id','')}")[:20]
            for spec in variants:
                variant_smarts = clean_text(spec["reaction_smarts"])
                out = dict(base)
                out["template_origin"] = "exact_reaction_derived_smarts" if str(out.get("abstracted_from_exact_reaction", "")).lower() in {"true", "1", "yes"} else "source_reaction_smarts"
                out["template_scope"] = "generalized_template"
                out["template_generalization"] = "generalized_directional"
                out["source_reaction_smarts"] = smarts
                out["reaction_smarts"] = variant_smarts
                out["direction_handling"] = spec["direction_handling"]
                out["direction_variant"] = spec["direction_variant"]
                out["reverse_template_hash"] = spec["reverse_template_hash"]
                out["reversible_group_id"] = group_id if len(variants) > 1 else ""
                out["predictive_rule_use"] = "true"
                out["anchor_edge_use"] = "false"
                out["template_extraction_status"] = "imported_source_template_directionalized"
                if spec["direction_handling"] == "split_reversible_reverse":
                    _reverse_example_pair(out)
                _set_direction_qc(out)
                out["template_id"] = f"TPL_{len(rows)+1:09d}"
                out["template_hash"] = sha256_text(variant_smarts)[:20]
                rows.append(out)
        elif can_rxn:
            out = dict(base)
            out["template_origin"] = "exact_substrate_product_pair"
            out["template_scope"] = "exact_anchor"
            out["template_generalization"] = "exact_pair_not_generalized"
            out["source_reaction_smarts"] = ""
            out["reverse_template_hash"] = ""
            out["predictive_rule_use"] = "true" if allow_exact_pairs_as_predictive_rules else "false"
            out["anchor_edge_use"] = "true"
            out["template_extraction_status"] = "exact_pair_retained_not_generalized"
            out["direction_handling"] = out.get("direction_handling", "") or "exact_pair_left_to_right"
            out["direction_variant"] = out.get("direction_variant", "") or "forward"
            _set_direction_qc(out)
            out["reaction_smarts"] = ""
            out["template_id"] = f"TPL_{len(rows)+1:09d}"
            out["template_hash"] = sha256_text(f"EXACT::{can_rxn}")[:20]
            rows.append(out)
        else:
            out = dict(base)
            out["template_origin"] = "none"
            out["template_scope"] = "no_template"
            out["template_generalization"] = "none"
            out["source_reaction_smarts"] = ""
            out["reverse_template_hash"] = ""
            out["predictive_rule_use"] = "false"
            out["anchor_edge_use"] = "false"
            out["direction_handling"] = out.get("direction_handling", "") or "missing"
            out["direction_variant"] = out.get("direction_variant", "")
            _set_direction_qc(out)
            out["template_extraction_status"] = "missing_template_and_pair"
            out["direction_handling"] = out.get("direction_handling", "") or "missing_template_no_direction_handling"
            out["reaction_smarts"] = ""
            out["template_id"] = f"TPL_{len(rows)+1:09d}"
            out["template_hash"] = sha256_text(f"MISSING::{i}")[:20]
            rows.append(out)
    df = pd.DataFrame(rows).fillna("")
    for col in TEMPLATE_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[TEMPLATE_COLUMNS]


def qc_templates(templates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in templates.fillna("").iterrows():
        out = row.to_dict()
        scope = clean_text(out.get("template_scope"))
        smarts = clean_text(out.get("reaction_smarts"))
        can_rxn = clean_text(out.get("canonical_reaction_smiles"))
        if scope == "generalized_template":
            ok, note = reaction_smarts_valid(smarts)
            out["template_qc_status"] = "ok" if ok else ("rdkit_unavailable" if note == "rdkit_unavailable" else "invalid")
            out["template_qc_note"] = note
        elif scope == "exact_anchor" and can_rxn:
            out["template_qc_status"] = "exact_pair_ok"
            out["template_qc_note"] = "exact substrate-product pair retained as anchor; not a generalized predictive template"
        else:
            out["template_qc_status"] = "no_template"
            out["template_qc_note"] = "missing generalized template and exact pair"
        rows.append(out)
    df = pd.DataFrame(rows).fillna("")
    for col in TEMPLATE_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[TEMPLATE_COLUMNS]


def deduplicate_templates(qc_df: pd.DataFrame) -> pd.DataFrame:
    if qc_df.empty:
        return qc_df.copy()
    df = qc_df.fillna("").copy()
    merge_cols = [
        "source_database", "evidence_layer", "source_reaction_id", "record_id", "ec_numbers", "template_ec_candidates",
        "database_ec_candidates", "ec_prior_candidates", "rhea_ids", "kegg_ids", "metanetx_ids", "reaction_type_source",
        "reaction_subtype_source", "cofactor_or_donor_class", "enzyme_name", "protein_ids", "source_file", "source_evidence_text",
        "donor_class", "acceptor_atom_class", "transferred_group_class", "main_pair_projection_method", "main_pair_projection_note",
        "abstracted_from_exact_reaction", "derived_from_exact_anchor", "rxnmapper_confidence",
        "rdchiral_extraction_status", "abstracted_smarts_applies_to_original_pair",
        "exact_abstraction_qc_status", "benchmark_exclusion_flag", "direction_handling", "direction_variant", "reversible_group_id",
        "source_reaction_smarts", "reverse_template_hash",
        "consensus_group_id", "consensus_generation_mode", "consensus_evidence_rows",
        "consensus_source_database_count", "consensus_evidence_layer_support", "consensus_qc_status",
        "consensus_representative_rule_ids", "consensus_supporting_reaction_types",
    ]

    # Most large releases are dominated by singleton template hashes. Keep those
    # rows directly and only perform expensive provenance joins for true
    # duplicates; this preserves the old semantics while avoiding one Python
    # group loop per unique rule.
    counts = df["template_hash"].value_counts(sort=False, dropna=False)
    out = df.drop_duplicates("template_hash", keep="first").set_index("template_hash", drop=False)
    out["template_count"] = counts.reindex(out.index).fillna(1).astype(int).astype(str)
    out["source_record_count"] = "1"

    duplicate_hashes = counts[counts > 1].index
    if len(duplicate_hashes):
        dup = df[df["template_hash"].isin(duplicate_hashes)].copy()
        grouped = dup.groupby("template_hash", sort=False, dropna=False)

        for col in merge_cols:
            if col in dup.columns:
                merged = grouped[col].agg(lambda s: join_values(s.astype(str).tolist()))
                out.loc[merged.index, col] = merged

        if "template_qc_status" in dup.columns:
            status_values = dup["template_qc_status"].astype(str)
            status = pd.Series("no_template", index=duplicate_hashes, dtype=object)
            for label in ["invalid", "rdkit_unavailable", "exact_pair_ok", "ok"]:
                mask = status_values.eq(label).groupby(dup["template_hash"]).any()
                status.loc[mask[mask].index] = label
            out.loc[status.index, "template_qc_status"] = status

        if "template_qc_note" in dup.columns:
            notes = grouped["template_qc_note"].agg(lambda s: join_values(s.astype(str).tolist()))
            out.loc[notes.index, "template_qc_note"] = notes

        if "record_id" in dup.columns:
            source_counts = grouped["record_id"].agg(lambda s: len(set(s.astype(str).tolist())))
            out.loc[source_counts.index, "source_record_count"] = source_counts.astype(int).astype(str)

        truthy_values = ["true", "1", "yes"]
        if "predictive_rule_use" in dup.columns:
            predictive = dup["predictive_rule_use"].astype(str).str.lower().isin(truthy_values).groupby(dup["template_hash"]).any()
            out.loc[predictive.index, "predictive_rule_use"] = predictive.map(lambda x: "true" if x else "false")
        if "anchor_edge_use" in dup.columns:
            anchor = dup["anchor_edge_use"].astype(str).str.lower().isin(truthy_values).groupby(dup["template_hash"]).any()
            out.loc[anchor.index, "anchor_edge_use"] = anchor.map(lambda x: "true" if x else "false")

    out = out.reset_index(drop=True).fillna("")
    for col in TEMPLATE_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    return out[TEMPLATE_COLUMNS]
