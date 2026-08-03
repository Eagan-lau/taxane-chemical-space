from __future__ import annotations

from pathlib import Path

import pandas as pd

from .schema import FAMILY_EVIDENCE_COLUMNS
from .utils import clean_text, read_yaml


def interproscan_to_family_evidence(
    interpro_tsv: str | Path,
    family_map_yaml: str | Path,
    output_tsv: str | Path,
    *,
    source_database: str = "ProteinDomainEvidence",
    source_reaction_id_from_seq_id: bool = False,
) -> pd.DataFrame:
    """Convert InterProScan TSV output to the external family-evidence schema.

    This utility uses a user-supplied YAML map from domain signatures to family names. It does not contain
    built-in enzyme-family motifs or reaction-type priors.
    """
    df = pd.read_csv(interpro_tsv, sep="\t", header=None, dtype=str, keep_default_na=False).fillna("")
    fmap = read_yaml(family_map_yaml).get("families", {}) or {}
    rows = []
    for _, row in df.iterrows():
        seq_id = clean_text(row.iloc[0]) if len(row) > 0 else ""
        sig_acc = clean_text(row.iloc[4]) if len(row) > 4 else ""
        sig_desc = clean_text(row.iloc[5]) if len(row) > 5 else ""
        ipr_acc = clean_text(row.iloc[11]) if len(row) > 11 else ""
        ipr_desc = clean_text(row.iloc[12]) if len(row) > 12 else ""
        text = f"{sig_acc} {sig_desc} {ipr_acc} {ipr_desc}".lower()
        for family, rule in fmap.items():
            accessions = {clean_text(x) for x in (rule.get("signature_accessions", []) or []) if clean_text(x)}
            patterns = [clean_text(x).lower() for x in (rule.get("signature_patterns", []) or []) if clean_text(x)]
            hit = sig_acc in accessions or ipr_acc in accessions or any(p in text for p in patterns)
            if hit:
                rows.append({
                    "evidence_id": f"DOMAIN_{len(rows)+1:09d}",
                    "match_type": "protein_id",
                    "match_value": seq_id,
                    "source_database": source_database,
                    "source_reaction_id": seq_id if source_reaction_id_from_seq_id else "",
                    "protein_id": seq_id,
                    "primary_family": clean_text(family),
                    "family": clean_text(family),
                    "family_role": clean_text(rule.get("family_role", "primary")),
                    "evidence_type": f"InterProScan:{sig_acc or ipr_acc}",
                    "confidence": clean_text(rule.get("confidence", 0.8)),
                    "evidence_source": source_database,
                    "provenance": f"{Path(interpro_tsv).name};{Path(family_map_yaml).name}",
                })
    out = pd.DataFrame(rows)
    for c in FAMILY_EVIDENCE_COLUMNS:
        if c not in out.columns:
            out[c] = ""
    out = out[FAMILY_EVIDENCE_COLUMNS].drop_duplicates()
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_tsv, sep="\t", index=False)
    return out
