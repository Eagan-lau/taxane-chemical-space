from __future__ import annotations

import pandas as pd

from .schema import ANCHOR_EDGE_COLUMNS


def build_anchor_edges(templates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if templates.empty:
        return pd.DataFrame(columns=ANCHOR_EDGE_COLUMNS)
    anchor_df = templates[templates.get("anchor_edge_use", "").astype(str).str.lower().isin(["true", "1", "yes"])]
    for i, row in anchor_df.fillna("").iterrows():
        rows.append({
            "anchor_edge_id": f"ANCHOR_{len(rows)+1:09d}",
            "source_database": row.get("source_database", ""),
            "source_reaction_id": row.get("source_reaction_id", ""),
            "enzyme_name": row.get("enzyme_name", ""),
            "substrate_smiles": row.get("main_substrate_smiles", ""),
            "product_smiles": row.get("main_product_smiles", ""),
            "canonical_reaction_smiles": row.get("canonical_reaction_smiles", ""),
            "ec_numbers": row.get("ec_numbers", ""),
            "evidence_layer": row.get("evidence_layer", ""),
            "direction": row.get("direction", ""),
            "source_reaction_reversibility": row.get("is_reversible", ""),
            "source_file": row.get("source_file", ""),
            "curated_taxol_anchor": row.get("curated_taxol_anchor", ""),
            "curated_pathway_name": row.get("curated_pathway_name", ""),
            "curated_pathway_step_id": row.get("curated_pathway_step_id", ""),
        })
    out = pd.DataFrame(rows)
    for col in ANCHOR_EDGE_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    return out[ANCHOR_EDGE_COLUMNS]
