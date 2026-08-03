from __future__ import annotations

import random

import pandas as pd

from .chem import canonical_reaction_smiles, rule_applies_to_pair
from .utils import clean_text


def _truthy_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([False] * len(df), index=df.index)
    return df[col].astype(str).str.lower().isin(["true", "1", "yes"])


def known_pathway_recall(rules: pd.DataFrame, source_reactions: pd.DataFrame) -> pd.DataFrame:
    known = source_reactions[_truthy_series(source_reactions, "curated_taxol_anchor")].copy()
    # Leakage control: curated exact anchors are excluded from the candidate rules.
    candidates = rules[
        (~_truthy_series(rules, "curated_taxol_anchor"))
        & (rules.get("template_scope", "").astype(str).eq("generalized_template"))
        & (rules.get("predictive_rule_use", "").astype(str).str.lower().isin(["true", "1", "yes"]))
    ].copy()
    rows = []
    for _, k in known.fillna("").iterrows():
        sub = clean_text(k.get("substrate_smiles", ""))
        prod = clean_text(k.get("product_smiles", ""))
        exact = canonical_reaction_smiles(f"{sub}>>{prod}") if sub and prod else ""
        hit_ids = []
        best_id = ""
        best_score = 0.0
        for _, r in candidates.fillna("").iterrows():
            smarts = clean_text(r.get("reaction_smarts", ""))
            applies = bool(smarts and sub and prod and rule_applies_to_pair(smarts, sub, prod))
            if applies:
                rid = clean_text(r.get("rule_id", ""))
                hit_ids.append(rid)
                try:
                    score = float(r.get("final_rule_confidence", 0) or 0)
                except Exception:
                    score = 0.0
                if score > best_score:
                    best_score = score
                    best_id = rid
        rows.append({
            "known_reaction_id": k.get("source_reaction_id", ""),
            "enzyme_name": k.get("enzyme_name", ""),
            "known_ec": k.get("ec_numbers", ""),
            "canonical_known_reaction_smiles": exact,
            "recovered_by_external_generalized_rule": bool(hit_ids),
            "best_rule_id": best_id,
            "best_rule_confidence": best_score,
            "all_hit_rule_ids": ";".join(hit_ids),
            "leakage_control": "curated_taxol_exact_anchors_excluded;only_external_generalized_templates_tested",
        })
    return pd.DataFrame(rows)


def decoy_pair_benchmark(rules: pd.DataFrame, source_reactions: pd.DataFrame, max_decoys: int = 500, seed: int = 13) -> pd.DataFrame:
    known = source_reactions[_truthy_series(source_reactions, "curated_taxol_anchor")].copy()
    cols = ["decoy_id", "substrate_reaction_id", "product_reaction_id", "matched_rule_count", "max_rule_confidence", "matched_rule_ids"]
    if known.empty:
        return pd.DataFrame(columns=cols)
    candidates = rules[
        (~_truthy_series(rules, "curated_taxol_anchor"))
        & (rules.get("template_scope", "").astype(str).eq("generalized_template"))
        & (rules.get("predictive_rule_use", "").astype(str).str.lower().isin(["true", "1", "yes"]))
    ].copy()
    records = known.fillna("").to_dict(orient="records")
    rng = random.Random(seed)
    n = min(max_decoys, max(0, len(records) * (len(records) - 1)))
    rows = []
    seen = set()
    attempts = 0
    while len(rows) < n and attempts < n * 20 + 100:
        attempts += 1
        a, b = rng.sample(records, 2)
        key = (a.get("source_reaction_id", ""), b.get("source_reaction_id", ""))
        if key in seen:
            continue
        seen.add(key)
        sub = clean_text(a.get("substrate_smiles", ""))
        prod = clean_text(b.get("product_smiles", ""))
        hits = []
        max_conf = 0.0
        for _, r in candidates.fillna("").iterrows():
            smarts = clean_text(r.get("reaction_smarts", ""))
            if smarts and sub and prod and rule_applies_to_pair(smarts, sub, prod):
                hits.append(r.get("rule_id", ""))
                try:
                    max_conf = max(max_conf, float(r.get("final_rule_confidence", 0) or 0))
                except Exception:
                    pass
        rows.append({
            "decoy_id": f"DECOY_{len(rows)+1:06d}",
            "substrate_reaction_id": a.get("source_reaction_id", ""),
            "product_reaction_id": b.get("source_reaction_id", ""),
            "matched_rule_count": len(hits),
            "max_rule_confidence": max_conf,
            "matched_rule_ids": ";".join(hits),
        })
    return pd.DataFrame(rows, columns=cols)
