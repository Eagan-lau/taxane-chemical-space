from __future__ import annotations

import json
from typing import Iterable

from .utils import clean_text, join_values

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    RDKIT_AVAILABLE = True
except Exception:  # pragma: no cover
    Chem = None
    AllChem = None
    Descriptors = None
    rdMolDescriptors = None
    RDKIT_AVAILABLE = False


def canonical_smiles(smiles: str, *, isomeric: bool = True) -> str:
    s = clean_text(smiles)
    if not s:
        return ""
    if not RDKIT_AVAILABLE:
        return s
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)


def mol_from_smiles(smiles: str):
    if not RDKIT_AVAILABLE:
        return None
    s = clean_text(smiles)
    if not s:
        return None
    return Chem.MolFromSmiles(s)


def heavy_atom_count(smiles: str) -> int:
    mol = mol_from_smiles(smiles)
    if mol is not None:
        return int(mol.GetNumHeavyAtoms())
    # Conservative fallback for environments without RDKit.
    import re
    return len(re.findall(r"Br|Cl|Si|Se|[BCNOFPSIbcnops]", clean_text(smiles)))


def split_components(side_smiles: str) -> list[str]:
    s = clean_text(side_smiles)
    return [part for part in s.split(".") if clean_text(part)] if s else []


def canonicalize_components(components: Iterable[str]) -> list[str]:
    out = []
    for comp in components:
        can = canonical_smiles(comp)
        if can:
            out.append(can)
    return sorted(out)


def parse_reaction_smiles(rxn: str) -> tuple[list[str], list[str]]:
    s = clean_text(rxn)
    if not s:
        return [], []
    if ">>" in s:
        left, right = s.split(">>", 1)
    elif s.count(">") == 2:
        left, _agents, right = s.split(">", 2)
    else:
        return [], []
    return split_components(left), split_components(right)


def build_reaction_smiles(reactants: Iterable[str], products: Iterable[str]) -> str:
    left = ".".join(canonicalize_components(reactants))
    right = ".".join(canonicalize_components(products))
    return f"{left}>>{right}" if left and right else ""


def canonical_reaction_smiles(rxn: str) -> str:
    reactants, products = parse_reaction_smiles(rxn)
    return build_reaction_smiles(reactants, products)


def split_reaction_transform(transform: str) -> tuple[str, str, str]:
    """Split a reaction SMILES/SMARTS string into left, agent, right fields."""
    s = clean_text(transform)
    if not s:
        return "", "", ""
    if ">>" in s:
        left, right = s.split(">>", 1)
        return clean_text(left), "", clean_text(right)
    if s.count(">") == 2:
        left, agents, right = s.split(">", 2)
        return clean_text(left), clean_text(agents), clean_text(right)
    return "", "", ""


def reverse_reaction_transform(transform: str) -> str:
    """Reverse a directional reaction SMILES/SMARTS transform.

    This is a syntactic direction-correction helper. It does not imply that the
    biochemical reaction is reversible or that the same enzyme catalyzes both
    directions.
    """
    left, agents, right = split_reaction_transform(transform)
    if not left or not right:
        return ""
    if agents:
        return f"{right}>{agents}>{left}"
    return f"{right}>>{left}"


# Backward-compatible aliases used by older code/tests.
def reverse_reaction_text(reaction: str) -> str:
    return reverse_reaction_transform(reaction)


def reverse_reaction_smarts(smarts: str) -> str:
    return reverse_reaction_transform(smarts)


def reverse_reaction_smiles(rxn: str) -> str:
    reactants, products = parse_reaction_smiles(rxn)
    return build_reaction_smiles(products, reactants)


def reaction_smarts_valid(smarts: str) -> tuple[bool, str]:
    s = clean_text(smarts)
    if not s:
        return False, "empty"
    if not RDKIT_AVAILABLE:
        return False, "rdkit_unavailable"
    last = "invalid"
    for use_smiles in (False, True):
        try:
            rxn = AllChem.ReactionFromSmarts(s, useSmiles=use_smiles)
            if rxn is not None and rxn.GetNumReactantTemplates() > 0 and rxn.GetNumProductTemplates() > 0:
                return True, "ok_useSmiles" if use_smiles else "ok"
        except Exception as exc:
            last = str(exc)
    return False, last


def rule_applies_to_pair(reaction_smarts: str, substrate_smiles: str, product_smiles: str, max_products: int = 500) -> bool:
    if not RDKIT_AVAILABLE:
        return False
    smarts = clean_text(reaction_smarts)
    if not smarts:
        return False
    substrate = mol_from_smiles(substrate_smiles)
    target = canonical_smiles(product_smiles)
    if substrate is None or not target:
        return False
    rxn = None
    for use_smiles in (False, True):
        try:
            rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=use_smiles)
            if rxn is not None:
                break
        except Exception:
            pass
    if rxn is None or rxn.GetNumReactantTemplates() != 1:
        return False
    try:
        outcomes = rxn.RunReactants((substrate,))
    except Exception:
        return False
    checked = 0
    for outcome in outcomes:
        for mol in outcome:
            checked += 1
            if checked > max_products:
                return False
            try:
                Chem.SanitizeMol(mol)
                if Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) == target:
                    return True
            except Exception:
                continue
    return False


def _atom_counts(smiles: str) -> dict[str, int]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {}
    formula = rdMolDescriptors.CalcMolFormula(mol)
    import re
    counts: dict[str, int] = {}
    for elem, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[elem] = counts.get(elem, 0) + int(num or 1)
    return counts


def reaction_delta(substrate_smiles: str, product_smiles: str) -> tuple[str, str]:
    sub = canonical_smiles(substrate_smiles)
    prod = canonical_smiles(product_smiles)
    if not sub or not prod or not RDKIT_AVAILABLE:
        return "", "{}"
    smol, pmol = mol_from_smiles(sub), mol_from_smiles(prod)
    if smol is None or pmol is None:
        return "", "{}"
    s_counts, p_counts = _atom_counts(sub), _atom_counts(prod)
    elems = sorted(set(s_counts) | set(p_counts))
    atom_delta = {e: p_counts.get(e, 0) - s_counts.get(e, 0) for e in elems if p_counts.get(e, 0) - s_counts.get(e, 0) != 0}
    mass_delta = float(Descriptors.ExactMolWt(pmol) - Descriptors.ExactMolWt(smol))
    heavy_delta = int(pmol.GetNumHeavyAtoms() - smol.GetNumHeavyAtoms())
    ring_delta = int(rdMolDescriptors.CalcNumRings(pmol) - rdMolDescriptors.CalcNumRings(smol))
    payload = {
        "canonical_substrate_smiles": sub,
        "canonical_product_smiles": prod,
        "atom_delta": atom_delta,
        "heavy_atom_delta": heavy_delta,
        "exact_mass_delta": round(mass_delta, 6),
        "ring_count_delta": ring_delta,
    }
    sig_parts = [f"{k}{v:+d}" for k, v in atom_delta.items()]
    sig_parts.append(f"HA{heavy_delta:+d}")
    if ring_delta:
        sig_parts.append(f"R{ring_delta:+d}")
    sig_parts.append(f"M{mass_delta:+.3f}")
    return "|".join(sig_parts), json.dumps(payload, ensure_ascii=False, sort_keys=True)


def replay_reaction_smarts_on_pair(reaction_smarts: str, substrate_smiles: str, product_smiles: str, max_products: int = 500) -> tuple[bool, str]:
    """Validate that a directional reaction SMARTS reproduces an exact substrate-product pair."""
    if not RDKIT_AVAILABLE:
        return False, "rdkit_unavailable"
    if not clean_text(reaction_smarts):
        return False, "empty_reaction_smarts"
    if not clean_text(substrate_smiles) or not clean_text(product_smiles):
        return False, "missing_substrate_or_product"
    if rule_applies_to_pair(reaction_smarts, substrate_smiles, product_smiles, max_products=max_products):
        return True, "replay_product_match"
    return False, "replay_product_not_generated"
