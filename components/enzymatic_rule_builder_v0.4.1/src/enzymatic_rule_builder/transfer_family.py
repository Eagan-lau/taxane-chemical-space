from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .chem import heavy_atom_count, reaction_smarts_valid, replay_reaction_smarts_on_pair
from .schema import TEMPLATE_COLUMNS
from .source_layers import normalize_evidence_layers
from .utils import clean_text, ensure_dir, join_values, sha256_text, split_multi_value, write_json

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    RDKIT_AVAILABLE = True
except Exception:  # pragma: no cover
    Chem = None
    rdFMCS = None
    RDKIT_AVAILABLE = False


TRANSFER_FAMILY_REPORT_COLUMNS = [
    "record_id",
    "source_database",
    "source_reaction_id",
    "evidence_layer",
    "transfer_family",
    "transfer_subtype",
    "donor_class",
    "acceptor_atom_class",
    "reaction_delta_fingerprint",
    "main_substrate_smiles",
    "main_product_smiles",
    "reaction_smarts",
    "template_hash",
    "extraction_status",
    "replay_pass",
    "replay_note",
    "skip_reason",
]


def _text(row: dict[str, Any], key: str) -> str:
    return clean_text(row.get(key, ""))


def _load_delta(row: dict[str, Any]) -> dict[str, Any]:
    text = _text(row, "reaction_delta_json")
    if not text:
        return {}
    try:
        data = json.loads(text)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _atom_delta(row: dict[str, Any]) -> tuple[dict[str, int], int, float]:
    delta = _load_delta(row)
    atoms = delta.get("atom_delta", {}) if isinstance(delta, dict) else {}
    atoms = atoms if isinstance(atoms, dict) else {}
    parsed = {str(k): int(v or 0) for k, v in atoms.items()}
    heavy = int(delta.get("heavy_atom_delta", 0) or 0) if isinstance(delta, dict) else 0
    mass = float(delta.get("exact_mass_delta", 0.0) or 0.0) if isinstance(delta, dict) else 0.0
    return parsed, heavy, mass


def _close(value: float, target: float, tolerance: float = 1.0) -> bool:
    return abs(float(value) - float(target)) <= tolerance


def _sugar_subtype(atoms: dict[str, int], mass: float) -> str:
    c = atoms.get("C", 0)
    h = atoms.get("H", 0)
    o = atoms.get("O", 0)
    if c == 5 and 7 <= h <= 9 and o == 4 and _close(mass, 132.0423, 2.0):
        return "O_pentosylation_or_xylosylation"
    if c == 6 and 9 <= h <= 11 and o == 5 and _close(mass, 162.0528, 2.0):
        return "O_hexosylation"
    if c == 6 and 7 <= h <= 9 and o == 6 and _close(mass, 176.0321, 2.0):
        return "O_glucuronidation_or_uronic_acid_glycosylation"
    if c == 6 and 9 <= h <= 11 and o == 4 and _close(mass, 146.0579, 2.0):
        return "O_deoxyhexosylation"
    if c >= 5 and o >= 4 and 120 <= mass <= 220:
        return "O_glycosylation"
    return ""


def _acyl_subtype(atoms: dict[str, int], mass: float) -> str:
    c = atoms.get("C", 0)
    h = atoms.get("H", 0)
    o = atoms.get("O", 0)
    if c == 2 and 1 <= h <= 3 and o == 1 and _close(mass, 42.0106, 2.0):
        return "O_acetylation"
    if c == 7 and 3 <= h <= 5 and o == 1 and _close(mass, 104.0262, 3.0):
        return "O_benzoylation"
    if c == 9 and 5 <= h <= 7 and o == 1 and _close(mass, 130.0419, 3.0):
        return "O_cinnamoylation_or_hydroxycinnamoylation"
    if c >= 2 and o >= 1 and 35 <= mass <= 180:
        return "O_acylation"
    return ""


def _methyl_subtype(atoms: dict[str, int], mass: float) -> str:
    if atoms.get("C", 0) == 1 and 1 <= atoms.get("H", 0) <= 3 and _close(mass, 14.0157, 2.0):
        return "O_methylation"
    return ""


def _phosphoryl_subtype(atoms: dict[str, int], mass: float) -> str:
    if atoms.get("P", 0) == 1 and atoms.get("O", 0) >= 3 and 75 <= mass <= 100:
        return "O_phosphorylation"
    return ""


def _infer_transfer_family(row: dict[str, Any]) -> tuple[str, str, str]:
    atoms, heavy, mass = _atom_delta(row)
    if heavy <= 0 or mass <= 0:
        return "", "", "non_addition_delta"

    donor_text = " ".join(
        [
            _text(row, "donor_class"),
            _text(row, "cofactor_or_donor_class"),
            _text(row, "external_participant_roles"),
            _text(row, "transferred_group"),
            _text(row, "transferred_group_class"),
        ]
    ).lower()
    if not donor_text:
        return "", "", "missing_donor_or_transfer_annotation"

    sugar = _sugar_subtype(atoms, mass)
    if sugar and any(token in donor_text for token in ["udp_sugar", "nucleotide_sugar", "glycosyl", "sugar"]):
        return "glycosylation", sugar, ""

    acyl = _acyl_subtype(atoms, mass)
    if acyl and any(token in donor_text for token in ["coa", "acyl", "acetyl", "benzoyl", "cinnamoyl"]):
        return "acylation", acyl, ""

    methyl = _methyl_subtype(atoms, mass)
    if methyl and any(token in donor_text for token in ["sam", "methyl"]):
        return "methylation", methyl, ""

    phosphoryl = _phosphoryl_subtype(atoms, mass)
    if phosphoryl and any(token in donor_text for token in ["atp", "nucleotide", "phosphoryl", "phosphate"]):
        return "phosphorylation", phosphoryl, ""

    return "", "", "donor_annotation_delta_mismatch"


def _acceptor_reactant_smarts(element: str, acceptor_class: str) -> str:
    elem = clean_text(element)
    cls = clean_text(acceptor_class).lower()
    if elem == "O":
        if "carboxyl" in cls:
            return "[O;H1,-1:1]"
        return "[O;H1;+0:1]"
    if elem == "N":
        return "[N;H1,H2,H3;+0:1]"
    if elem == "S":
        return "[S;H1;+0:1]"
    return ""


def _acceptor_type_prefix(element: str, acceptor_class: str) -> str:
    elem = clean_text(element)
    cls = clean_text(acceptor_class).lower()
    if elem == "N" or cls.startswith("amine"):
        return "N"
    if elem == "S" or cls.startswith("thiol"):
        return "S"
    return "O"


def _product_fragment_from_main_pair(
    substrate_smiles: str,
    product_smiles: str,
    *,
    acceptor_class: str = "",
    mcs_timeout: int = 10,
) -> tuple[str, str, str]:
    if not RDKIT_AVAILABLE:
        return "", "", "rdkit_unavailable"
    substrate = Chem.MolFromSmiles(clean_text(substrate_smiles))
    product = Chem.MolFromSmiles(clean_text(product_smiles))
    if substrate is None or product is None:
        return "", "", "invalid_main_pair_smiles"
    if product.GetNumHeavyAtoms() <= substrate.GetNumHeavyAtoms():
        return "", "", "product_not_larger_than_substrate"

    mcs = rdFMCS.FindMCS(
        [substrate, product],
        timeout=int(mcs_timeout),
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchValences=True,
    )
    if mcs.canceled or mcs.numAtoms <= 0:
        return "", "", "mcs_failed_or_timed_out"
    query = Chem.MolFromSmarts(mcs.smartsString)
    if query is None:
        return "", "", "mcs_query_failed"
    substrate_match = substrate.GetSubstructMatch(query)
    product_match = product.GetSubstructMatch(query)
    if not substrate_match or not product_match or len(substrate_match) != len(product_match):
        return "", "", "mcs_match_failed"

    old_product_atoms = set(product_match)
    product_to_substrate = {pidx: sidx for sidx, pidx in zip(substrate_match, product_match)}
    attachments: list[tuple[int, int]] = []
    for pidx in sorted(old_product_atoms):
        patom = product.GetAtomWithIdx(pidx)
        satom = substrate.GetAtomWithIdx(product_to_substrate[pidx])
        if patom.GetSymbol() not in {"O", "N", "S"} or patom.GetSymbol() != satom.GetSymbol():
            continue
        for nbr in patom.GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx not in old_product_atoms:
                attachments.append((pidx, nidx))
    if len(attachments) != 1:
        return "", "", "ambiguous_or_missing_single_attachment"

    attach_idx, new_neighbor_idx = attachments[0]
    attach_atom = product.GetAtomWithIdx(attach_idx)
    reactant_smarts = _acceptor_reactant_smarts(attach_atom.GetSymbol(), acceptor_class)
    if not reactant_smarts:
        return "", "", "unsupported_acceptor_atom"

    stack = [new_neighbor_idx]
    new_atoms: set[int] = set()
    while stack:
        idx = stack.pop()
        if idx in new_atoms or idx in old_product_atoms:
            continue
        new_atoms.add(idx)
        for nbr in product.GetAtomWithIdx(idx).GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx != attach_idx and nidx not in old_product_atoms:
                stack.append(nidx)
    if not (1 <= len(new_atoms) <= 30):
        return "", "", "transferred_fragment_size_out_of_range"

    for atom in product.GetAtoms():
        atom.SetAtomMapNum(0)
    attach_atom.SetAtomMapNum(1)
    atoms_to_use = sorted(new_atoms | {attach_idx})
    try:
        product_fragment = Chem.MolFragmentToSmiles(
            product,
            atomsToUse=atoms_to_use,
            rootedAtAtom=attach_idx,
            canonical=False,
            isomericSmiles=True,
        )
    finally:
        for atom in product.GetAtoms():
            atom.SetAtomMapNum(0)
    if ":1" not in product_fragment:
        return "", "", "mapped_acceptor_missing_from_product_fragment"
    return f"{reactant_smarts}>>{product_fragment}", attach_atom.GetSymbol(), "ok"


def _directional_reaction_type(family: str, subtype: str, acceptor_prefix: str) -> tuple[str, str]:
    st = clean_text(subtype)
    if family == "glycosylation":
        return f"{acceptor_prefix}_glycosylation", st
    if family == "acylation":
        return st.replace("O_", f"{acceptor_prefix}_", 1) if st.startswith("O_") else f"{acceptor_prefix}_acylation", st
    if family == "methylation":
        return st.replace("O_", f"{acceptor_prefix}_", 1) if st.startswith("O_") else f"{acceptor_prefix}_methylation", st
    if family == "phosphorylation":
        return st.replace("O_", f"{acceptor_prefix}_", 1) if st.startswith("O_") else f"{acceptor_prefix}_phosphorylation", st
    return family, st


def _report(row: dict[str, Any], **updates: Any) -> dict[str, Any]:
    out = {
        "record_id": _text(row, "record_id"),
        "source_database": _text(row, "source_database"),
        "source_reaction_id": _text(row, "source_reaction_id"),
        "evidence_layer": _text(row, "evidence_layer"),
        "transfer_family": "",
        "transfer_subtype": "",
        "donor_class": _text(row, "donor_class") or _text(row, "cofactor_or_donor_class"),
        "acceptor_atom_class": _text(row, "acceptor_atom_class"),
        "reaction_delta_fingerprint": _text(row, "reaction_delta_fingerprint"),
        "main_substrate_smiles": _text(row, "main_substrate_smiles"),
        "main_product_smiles": _text(row, "main_product_smiles"),
        "reaction_smarts": "",
        "template_hash": "",
        "extraction_status": "",
        "replay_pass": "",
        "replay_note": "",
        "skip_reason": "",
    }
    out.update({k: clean_text(v) for k, v in updates.items()})
    return out


def build_transfer_family_consensus_templates(
    main_pairs: pd.DataFrame,
    *,
    output_dir: str | Path | None = None,
    mcs_timeout: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Build replay-validated donor-transfer family SMARTS from projected main pairs.

    The input rows must already contain main substrate/product projection fields.
    This function does not use enzyme-family hard-coding: it requires explicit
    donor/transfer annotation plus a compatible structural delta, extracts the
    transferred product fragment from each exact main pair, and keeps only SMARTS
    that replay the source main pair.
    """
    if main_pairs.empty:
        empty = pd.DataFrame(columns=TEMPLATE_COLUMNS)
        report = pd.DataFrame(columns=TRANSFER_FAMILY_REPORT_COLUMNS)
        summary = {"enabled": True, "input_main_pairs": 0, "candidate_rows": 0, "released_template_rows": 0}
        if output_dir:
            out = ensure_dir(output_dir)
            report.to_csv(out / "transfer_family_consensus.audit.tsv", sep="\t", index=False)
            write_json(out / "transfer_family_consensus.summary.json", summary)
        return empty, report, summary

    candidates: list[dict[str, Any]] = []
    reports: list[dict[str, Any]] = []
    for _, series in main_pairs.fillna("").iterrows():
        row = series.to_dict()
        substrate = _text(row, "main_substrate_smiles")
        product = _text(row, "main_product_smiles")
        if not substrate or not product:
            reports.append(_report(row, extraction_status="skipped", skip_reason="missing_main_pair"))
            continue
        family, subtype, reason = _infer_transfer_family(row)
        if not family:
            if reason != "missing_donor_or_transfer_annotation":
                reports.append(_report(row, extraction_status="skipped", skip_reason=reason))
            continue

        smarts, acceptor_element, status = _product_fragment_from_main_pair(
            substrate,
            product,
            acceptor_class=_text(row, "acceptor_atom_class"),
            mcs_timeout=mcs_timeout,
        )
        if status != "ok":
            reports.append(
                _report(
                    row,
                    transfer_family=family,
                    transfer_subtype=subtype,
                    extraction_status="skipped",
                    skip_reason=status,
                )
            )
            continue
        ok, note = reaction_smarts_valid(smarts)
        if not ok:
            reports.append(
                _report(
                    row,
                    transfer_family=family,
                    transfer_subtype=subtype,
                    reaction_smarts=smarts,
                    extraction_status="skipped",
                    skip_reason=f"invalid_reaction_smarts:{note}",
                )
            )
            continue
        replay_ok, replay_note = replay_reaction_smarts_on_pair(smarts, substrate, product, max_products=1000)
        if not replay_ok:
            reports.append(
                _report(
                    row,
                    transfer_family=family,
                    transfer_subtype=subtype,
                    reaction_smarts=smarts,
                    extraction_status="skipped",
                    replay_pass="false",
                    replay_note=replay_note,
                    skip_reason="replay_failed_against_source_main_pair",
                )
            )
            continue

        acceptor_prefix = _acceptor_type_prefix(acceptor_element, _text(row, "acceptor_atom_class"))
        reaction_type, reaction_subtype = _directional_reaction_type(family, subtype, acceptor_prefix)
        template_hash = sha256_text(smarts)[:20]
        out = dict(row)
        out["reaction_smarts"] = smarts
        out["reaction_type_source"] = reaction_type
        out["reaction_subtype_source"] = reaction_subtype
        out["template_id"] = f"TFAM_TPL_{len(candidates)+1:09d}"
        out["template_hash"] = template_hash
        out["template_origin"] = "transfer_family_consensus_smarts"
        out["template_scope"] = "generalized_template"
        out["template_generalization"] = "role_aware_main_pair_transfer_family_consensus"
        out["source_reaction_smarts"] = ""
        out["reverse_template_hash"] = ""
        out["predictive_rule_use"] = "true"
        out["anchor_edge_use"] = "false"
        out["template_extraction_status"] = "transfer_family_consensus_projected_from_main_pair"
        out["template_qc_status"] = ""
        out["template_qc_note"] = join_values([
            "transfer_family_replay_validated_against_source_main_pair",
            replay_note,
        ])
        out["direction_handling"] = "role_inferred_forward_transfer_from_external_donor"
        out["direction_variant"] = "forward"
        out["normalized_direction"] = "substrate_to_product"
        out["direction_qc_status"] = "direction_qc_ok"
        out["direction_qc_note"] = (
            "External donor on the reactant side plus transferred-group gain in the main product "
            "supports the emitted substrate-to-product transfer direction."
        )
        out["reaction_representation_scope"] = "role_aware_main_pair_transfer_family_consensus"
        out["abstracted_from_exact_reaction"] = "false"
        out["derived_from_exact_anchor"] = "false"
        out["abstracted_smarts_applies_to_original_pair"] = "true"
        out["exact_abstraction_qc_status"] = "pass"
        out["consensus_generation_mode"] = "transfer_family_main_pair_fragment_consensus"
        candidates.append(out)
        reports.append(
            _report(
                row,
                transfer_family=family,
                transfer_subtype=subtype,
                reaction_smarts=smarts,
                template_hash=template_hash,
                extraction_status="released",
                replay_pass="true",
                replay_note=replay_note,
            )
        )

    if not candidates:
        empty = pd.DataFrame(columns=TEMPLATE_COLUMNS)
        report = pd.DataFrame(reports, columns=TRANSFER_FAMILY_REPORT_COLUMNS).fillna("")
        summary = {
            "enabled": True,
            "input_main_pairs": int(len(main_pairs)),
            "candidate_rows": 0,
            "released_template_rows": 0,
            "unique_reaction_smarts": 0,
            "rdkit_available": bool(RDKIT_AVAILABLE),
        }
        if output_dir:
            out = ensure_dir(output_dir)
            report.to_csv(out / "transfer_family_consensus.audit.tsv", sep="\t", index=False)
            write_json(out / "transfer_family_consensus.summary.json", summary)
        return empty, report, summary

    df = pd.DataFrame(candidates).fillna("")
    for template_hash, group in df.groupby("template_hash", dropna=False):
        idx = group.index
        source_tokens = [x for value in group.get("source_database", pd.Series(dtype=str)).astype(str) for x in split_multi_value(value)]
        layer_support = normalize_evidence_layers(join_values(group.get("evidence_layer", pd.Series(dtype=str)).astype(str).tolist()))
        group_id = "TRANSFER_FAMILY_CLUSTER_" + sha256_text(str(template_hash))[:20]
        df.loc[idx, "consensus_group_id"] = group_id
        df.loc[idx, "consensus_evidence_rows"] = str(int(len(group)))
        df.loc[idx, "consensus_source_database_count"] = str(int(len(set(source_tokens))))
        df.loc[idx, "consensus_evidence_layer_support"] = layer_support
        df.loc[idx, "consensus_qc_status"] = "passed_transfer_family_replay_validation"
        df.loc[idx, "consensus_representative_rule_ids"] = join_values(group.get("record_id", pd.Series(dtype=str)).astype(str).tolist())
        df.loc[idx, "consensus_supporting_reaction_types"] = join_values(group.get("reaction_type_source", pd.Series(dtype=str)).astype(str).tolist())

    for col in TEMPLATE_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    templates = df[TEMPLATE_COLUMNS].copy()
    report = pd.DataFrame(reports, columns=TRANSFER_FAMILY_REPORT_COLUMNS).fillna("")
    summary = {
        "enabled": True,
        "input_main_pairs": int(len(main_pairs)),
        "audit_rows": int(len(report)),
        "candidate_rows": int(len(candidates)),
        "released_template_rows": int(len(templates)),
        "unique_reaction_smarts": int(templates["reaction_smarts"].nunique()) if len(templates) else 0,
        "released_by_reaction_type": templates["reaction_type_source"].astype(str).value_counts().to_dict() if len(templates) else {},
        "rdkit_available": bool(RDKIT_AVAILABLE),
        "policy": "role-aware donor rows; compatible transferred-group delta; single attachment; SMARTS replay must reproduce source main pair",
    }
    if output_dir:
        out = ensure_dir(output_dir)
        report.to_csv(out / "transfer_family_consensus.audit.tsv", sep="\t", index=False)
        templates.to_csv(out / "transfer_family_consensus.templates.raw.tsv", sep="\t", index=False)
        write_json(out / "transfer_family_consensus.summary.json", summary)
    return templates, report, summary
