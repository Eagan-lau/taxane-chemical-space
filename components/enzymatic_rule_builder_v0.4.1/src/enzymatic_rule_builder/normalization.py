from __future__ import annotations

from pathlib import Path

import pandas as pd

from .chem import reaction_delta
from .cofactors import CofactorRegistry
from .schema import MAIN_PAIR_COLUMNS
from .utils import sha256_text


def add_main_pairs(source_df: pd.DataFrame, cofactor_yaml: str | Path | None = None) -> pd.DataFrame:
    registry = CofactorRegistry.from_yaml(cofactor_yaml)
    rows = []
    for _, row in source_df.fillna("").iterrows():
        out = row.to_dict()
        out.update(registry.strip_reaction(
            out.get("reaction_smiles", ""),
            out.get("substrate_smiles", ""),
            out.get("product_smiles", ""),
            out.get("cofactor_or_donor_class", ""),
            out.get("direction", ""),
            out.get("is_reversible", ""),
        ))
        if str(out.get("direction_handling", "")).startswith("source_reversible") and not out.get("reversible_group_id"):
            out["reversible_group_id"] = sha256_text(
                "REVERSIBLE::" + str(out.get("record_id", "")) + "::" + str(out.get("reaction_smarts", "")) + "::" + str(out.get("canonical_reaction_smiles", ""))
            )[:20]
        fp, js = reaction_delta(out.get("main_substrate_smiles", ""), out.get("main_product_smiles", ""))
        out["reaction_delta_fingerprint"] = fp
        out["reaction_delta_json"] = js
        rows.append(out)
    df = pd.DataFrame(rows).fillna("")
    for col in MAIN_PAIR_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[MAIN_PAIR_COLUMNS]
