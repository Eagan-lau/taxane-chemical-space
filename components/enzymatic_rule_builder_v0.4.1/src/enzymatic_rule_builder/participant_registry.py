from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .chem import canonical_smiles, heavy_atom_count, parse_reaction_smiles
from .schema import PARTICIPANT_REGISTRY_COLUMNS
from .sources import load_datasets
from .utils import clean_text, ensure_dir, join_values, sha256_text, write_yaml


def _atom_token_count(smiles: str, token: str) -> int:
    return clean_text(smiles).count(token)


def _classify_recurring_participant(smiles: str, reactant_count: int, product_count: int) -> dict[str, str]:
    """Assign an auditable, structure-derived external-participant role.

    This is not an enzyme-family mapping.  It only marks recurring non-core
    reaction participants so main-pair extraction does not confuse cofactors or
    donors with the scaffold being transformed.
    """
    s = clean_text(smiles)
    h = heavy_atom_count(s)
    n_p = _atom_token_count(s, "P")
    n_s = _atom_token_count(s, "S")
    role = "external_participant_candidate"
    registry_class = "external_participant_candidate"
    transferred = ""
    leaving = ""
    mode = "recurring_participant_frequency"

    if h <= 3:
        role = "small_molecule_reagent_or_byproduct"
        registry_class = "small_molecule"
        mode = "small_recurring_participant"
    elif n_s >= 1 and n_p >= 1 and h >= 25:
        role = "acyl_CoA_or_CoA_carrier"
        registry_class = "CoA_carrier_or_acyl_CoA_donor"
        transferred = "acyl"
        leaving = "CoA"
        mode = "structure_contains_sulfur_phosphate_CoA_like_carrier"
    elif n_p >= 2 and h >= 25 and n_s == 0:
        role = "nucleotide_sugar_or_phosphate_cofactor"
        registry_class = "nucleotide_sugar_or_phosphate_donor"
        transferred = "glycosyl_or_phosphoryl"
        leaving = "nucleotide_diphosphate_or_ADP"
        mode = "structure_contains_multiph phosphate_nucleotide_like_carrier"
    elif n_s >= 1 and n_p == 0 and 15 <= h <= 45:
        role = "SAM_or_sulfur_methyl_donor_candidate"
        registry_class = "SAM_or_sulfur_methyl_donor"
        transferred = "methyl"
        leaving = "SAH_or_sulfur_carrier"
        mode = "structure_contains_sulfur_adenosyl_like_candidate"
    elif reactant_count > product_count * 2:
        role = "left_side_external_donor_candidate"
        registry_class = "left_side_external_donor"
        mode = "reactant_enriched_recurring_participant"
    elif product_count > reactant_count * 2:
        role = "right_side_external_product_candidate"
        registry_class = "right_side_external_product"
        mode = "product_enriched_recurring_participant"
    return {
        "role_class": role,
        "registry_class": registry_class,
        "transferred_group": transferred,
        "leaving_group_class": leaving,
        "role_assignment_mode": mode,
    }


def build_participant_registry_from_sources(
    source_df: pd.DataFrame,
    *,
    min_occurrence: int = 3,
    max_heavy_atoms: int = 80,
) -> pd.DataFrame:
    """Build a data-derived external-participant candidate registry.

    This intentionally does not hard-code cofactor names. It ranks recurring
    participants from the supplied source reactions. The resulting table/YAML is a
    provenance-tracked registry that can be reviewed and then supplied to
    `build-rules` as cofactor/donor evidence.
    """
    stats: dict[str, dict[str, Any]] = {}
    for _, row in source_df.fillna("").iterrows():
        # Source-curated substrate/product rows encode the node-level main transform,
        # not the full biochemical participant set; do not infer external participants
        # from those exact pairs.
        if clean_text(row.get("substrate_smiles", "")) and clean_text(row.get("product_smiles", "")):
            continue
        rxn = clean_text(row.get("reaction_smiles", ""))
        reactants, products = parse_reaction_smiles(rxn)
        # A single-reactant/single-product exact reaction supplies no evidence for
        # external cofactors/donors. Only multi-component reactions can contribute
        # recurring participant candidates.
        if len(reactants) + len(products) <= 2:
            continue
        for side, comps in [("reactant", reactants), ("product", products)]:
            for comp in comps:
                can = canonical_smiles(comp)
                if not can:
                    continue
                h = heavy_atom_count(can)
                rec = stats.setdefault(can, {
                    "participant_smiles": can,
                    "participant_hash": sha256_text(can)[:20],
                    "heavy_atom_count": h,
                    "occurrence_count": 0,
                    "reactant_count": 0,
                    "product_count": 0,
                    "source_databases": set(),
                    "source_reaction_ids": set(),
                })
                rec["occurrence_count"] += 1
                rec[f"{side}_count"] += 1
                if clean_text(row.get("source_database", "")):
                    rec["source_databases"].add(clean_text(row.get("source_database")))
                if clean_text(row.get("source_reaction_id", "")):
                    rec["source_reaction_ids"].add(clean_text(row.get("source_reaction_id")))
    rows = []
    for can, rec in stats.items():
        if rec["occurrence_count"] < int(min_occurrence):
            continue
        if rec["heavy_atom_count"] > int(max_heavy_atoms):
            continue
        source_dbs = sorted(rec["source_databases"])
        source_ids = sorted(rec["source_reaction_ids"])
        role_info = _classify_recurring_participant(can, int(rec["reactant_count"]), int(rec["product_count"]))
        # Frequency across multiple reactions/sources is a candidate signal, not final biological role truth.
        if len(source_dbs) >= 2 or rec["occurrence_count"] >= max(int(min_occurrence) * 2, 6):
            conf = 0.70
            mode = role_info["role_assignment_mode"] + ";recurring_participant_multi_source_or_high_frequency"
        else:
            conf = 0.50
            mode = role_info["role_assignment_mode"] + ";recurring_participant_single_source"
        rows.append({
            "participant_smiles": rec["participant_smiles"],
            "participant_hash": rec["participant_hash"],
            "heavy_atom_count": int(rec["heavy_atom_count"]),
            "occurrence_count": int(rec["occurrence_count"]),
            "reactant_count": int(rec["reactant_count"]),
            "product_count": int(rec["product_count"]),
            "source_database_count": int(len(source_dbs)),
            "source_databases": join_values(source_dbs),
            "source_reaction_ids": join_values(source_ids[:50]),
            "role_class": role_info["role_class"],
            "role_assignment_mode": mode,
            "role_confidence": f"{conf:.2f}",
            "registry_class": role_info["registry_class"] if role_info["registry_class"] != "external_participant_candidate" else f"participant_{rec['participant_hash']}",
            "transferred_group": role_info["transferred_group"],
            "leaving_group_class": role_info["leaving_group_class"],
            "provenance": "database_derived_recurring_participant_registry",
        })
    df = pd.DataFrame(rows, columns=PARTICIPANT_REGISTRY_COLUMNS).fillna("")
    if not df.empty:
        df = df.sort_values(["occurrence_count", "source_database_count", "heavy_atom_count"], ascending=[False, False, True]).reset_index(drop=True)
    return df


def registry_to_yaml_payload(df: pd.DataFrame) -> dict[str, Any]:
    participants = []
    for _, row in df.fillna("").iterrows():
        participants.append({
            "class": clean_text(row.get("registry_class", "external_participant_candidate")),
            "role": clean_text(row.get("role_class", "external_participant_candidate")),
            "smiles": [clean_text(row.get("participant_smiles", ""))],
            "confidence": clean_text(row.get("role_confidence", "")),
            "transferred_group": clean_text(row.get("transferred_group", "")),
            "leaving_group_class": clean_text(row.get("leaving_group_class", "")),
            "provenance": clean_text(row.get("provenance", "")),
            "source_databases": clean_text(row.get("source_databases", "")),
            "occurrence_count": int(row.get("occurrence_count", 0) or 0),
        })
    return {"participants": participants}


def build_participant_registry_from_manifest(
    manifest: str | Path,
    output_dir: str | Path,
    *,
    min_occurrence: int = 3,
    max_heavy_atoms: int = 80,
) -> dict[str, Any]:
    out = ensure_dir(output_dir)
    source_df, used = load_datasets(manifest)
    registry = build_participant_registry_from_sources(source_df, min_occurrence=min_occurrence, max_heavy_atoms=max_heavy_atoms)
    tsv = out / "participant_role_registry.tsv"
    yaml_path = out / "participant_role_registry.yaml"
    registry.to_csv(tsv, sep="\t", index=False)
    write_yaml(yaml_path, registry_to_yaml_payload(registry))
    return {
        "source_records": int(len(source_df)),
        "input_files": used,
        "participant_registry_records": int(len(registry)),
        "participant_role_registry_tsv": str(tsv.resolve()),
        "participant_role_registry_yaml": str(yaml_path.resolve()),
        "min_occurrence": int(min_occurrence),
        "max_heavy_atoms": int(max_heavy_atoms),
    }
