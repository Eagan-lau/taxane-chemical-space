from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass
from typing import Any


FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
DELTA_TOKEN = re.compile(r"^([A-Z][a-z]?)([+-]\d+)$")
NON_ELEMENT_DELTA_FIELDS = {"HA", "R", "M"}


def require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as exc:
        raise RuntimeError(
            "RDKit is required. Run this package with the configured aimd interpreter."
        ) from exc
    return Chem, AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors, rdMolStandardize


def stable_hash(text: str, length: int = 20) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:length]


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def falsey(value: Any) -> bool:
    return str(value).strip().lower() in {"0", "false", "f", "no", "n"}


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or str(value).strip().lower() in {"", "nan", "none"}:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int | None = None) -> int | None:
    try:
        if value is None or str(value).strip().lower() in {"", "nan", "none"}:
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def parse_formula(formula: str) -> dict[str, int]:
    result: dict[str, int] = {}
    for element, count in FORMULA_TOKEN.findall(str(formula).strip()):
        result[element] = result.get(element, 0) + int(count or 1)
    return result


def formula_delta(source_formula: str, product_formula: str) -> dict[str, int]:
    source = parse_formula(source_formula)
    product = parse_formula(product_formula)
    elements = set(source) | set(product)
    return {
        element: product.get(element, 0) - source.get(element, 0)
        for element in sorted(elements)
        if product.get(element, 0) != source.get(element, 0)
    }


def parse_reaction_delta_fingerprint(value: str) -> dict[str, int]:
    result: dict[str, int] = {}
    for token in str(value).split("|"):
        match = DELTA_TOKEN.match(token.strip())
        if match and match.group(1) not in NON_ELEMENT_DELTA_FIELDS:
            result[match.group(1)] = int(match.group(2))
    return result


def atom_delta_matches(
    observed: dict[str, int],
    expected: dict[str, int],
    *,
    hydrogen_tolerance: int = 2,
) -> bool:
    if not expected:
        return True
    elements = set(observed) | set(expected)
    for element in elements:
        tolerance = hydrogen_tolerance if element == "H" else 0
        if abs(observed.get(element, 0) - expected.get(element, 0)) > tolerance:
            return False
    return True


@dataclass(frozen=True)
class ReactionTemplateMetrics:
    compile_status: str
    compile_error: str
    n_reactants: int
    n_products: int
    reactant_atoms: int
    product_atoms: int
    mapped_reactant_atoms: int
    mapped_product_atoms: int
    generic_reactant_atoms: int
    reactant_pattern_fp_hex: str
    changed_mapped_atom_count: int
    inferred_independent_reaction_centers: int
    structural_element_delta: str
    reaction_edit_signature: str
    mapped_atom_retention: float
    mapping_coverage: float


def _template_element_counts(mol) -> dict[str, int]:
    result: dict[str, int] = {}
    for atom in mol.GetAtoms():
        atomic_number = atom.GetAtomicNum()
        if atomic_number <= 1:
            continue
        symbol = atom.GetSymbol()
        result[symbol] = result.get(symbol, 0) + 1
    return result


def _element_delta_text(reactant, product) -> str:
    left = _template_element_counts(reactant)
    right = _template_element_counts(product)
    tokens = []
    for element in sorted(set(left) | set(right)):
        delta = right.get(element, 0) - left.get(element, 0)
        if delta:
            tokens.append(f"{element}{delta:+d}")
    return "|".join(tokens)


def _mapped_bonds(mol) -> dict[tuple[int, int], str]:
    result: dict[tuple[int, int], str] = {}
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtom().GetAtomMapNum()
        end = bond.GetEndAtom().GetAtomMapNum()
        if begin > 0 and end > 0:
            key = tuple(sorted((begin, end)))
            result[key] = str(bond.GetBondType())
    return result


def _changed_mapped_atoms(reactant, product) -> set[int]:
    left = {
        atom.GetAtomMapNum(): atom
        for atom in reactant.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    right = {
        atom.GetAtomMapNum(): atom
        for atom in product.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    changed = set(left) ^ set(right)
    for map_number in set(left) & set(right):
        left_atom = left[map_number]
        right_atom = right[map_number]
        if (
            left_atom.GetAtomicNum(),
            left_atom.GetFormalCharge(),
            str(left_atom.GetChiralTag()),
        ) != (
            right_atom.GetAtomicNum(),
            right_atom.GetFormalCharge(),
            str(right_atom.GetChiralTag()),
        ):
            changed.add(map_number)

        left_unmapped = sorted(
            neighbor.GetAtomicNum()
            for neighbor in left_atom.GetNeighbors()
            if neighbor.GetAtomMapNum() <= 0
        )
        right_unmapped = sorted(
            neighbor.GetAtomicNum()
            for neighbor in right_atom.GetNeighbors()
            if neighbor.GetAtomMapNum() <= 0
        )
        if left_unmapped != right_unmapped:
            changed.add(map_number)

    left_bonds = _mapped_bonds(reactant)
    right_bonds = _mapped_bonds(product)
    for pair in set(left_bonds) | set(right_bonds):
        if left_bonds.get(pair) != right_bonds.get(pair):
            changed.update(pair)
    return changed


def _reaction_edit_signature(reactant, product, changed: set[int]) -> str:
    left_atoms = {
        atom.GetAtomMapNum(): atom
        for atom in reactant.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    right_atoms = {
        atom.GetAtomMapNum(): atom
        for atom in product.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    atom_transitions: Counter[str] = Counter()
    neighbor_delta_total: Counter[str] = Counter()
    chiral_change_count = 0
    for map_number in sorted(changed):
        left = left_atoms.get(map_number)
        right = right_atoms.get(map_number)

        def atom_state(atom) -> tuple[str, int]:
            if atom is None:
                return ("missing", 0)
            return (atom.GetSymbol(), atom.GetFormalCharge())

        def unmapped_counts(atom) -> Counter[str]:
            if atom is None:
                return Counter()
            return Counter(
                neighbor.GetSymbol() if neighbor.GetAtomicNum() else "*"
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomMapNum() <= 0
            )

        left_unmapped = unmapped_counts(left)
        right_unmapped = unmapped_counts(right)
        for symbol in set(left_unmapped) | set(right_unmapped):
            neighbor_delta_total[symbol] += (
                right_unmapped.get(symbol, 0) - left_unmapped.get(symbol, 0)
            )
        left_state = atom_state(left)
        right_state = atom_state(right)
        if left_state != right_state:
            atom_transitions[f"{left_state}>{right_state}"] += 1
        if (
            left is not None
            and right is not None
            and str(left.GetChiralTag()) != str(right.GetChiralTag())
        ):
            chiral_change_count += 1

    left_bonds = _mapped_bonds(reactant)
    right_bonds = _mapped_bonds(product)
    bond_edits = Counter(
        f"{left_bonds.get(pair, 'none')}>{right_bonds.get(pair, 'none')}"
        for pair in set(left_bonds) | set(right_bonds)
        if left_bonds.get(pair) != right_bonds.get(pair)
    )
    non_chiral_edit = bool(
        atom_transitions
        or bond_edits
        or any(neighbor_delta_total.values())
        or _element_delta_text(reactant, product)
    )
    payload = {
        "atom_transitions": sorted(atom_transitions.items()),
        "bond_edits": sorted(bond_edits.items()),
        "element_delta": _element_delta_text(reactant, product),
        "neighbor_delta": sorted(
            (key, value) for key, value in neighbor_delta_total.items() if value
        ),
        "chiral_only_changes": chiral_change_count if not non_chiral_edit else 0,
    }
    return stable_hash(
        json.dumps(payload, sort_keys=True, separators=(",", ":")),
        length=24,
    )


def _independent_center_count(reactant, product, changed: set[int]) -> int:
    if not changed:
        return 0
    adjacency = {map_number: set() for map_number in changed}
    for mol in (reactant, product):
        for bond in mol.GetBonds():
            begin = bond.GetBeginAtom().GetAtomMapNum()
            end = bond.GetEndAtom().GetAtomMapNum()
            if begin in changed and end in changed:
                adjacency[begin].add(end)
                adjacency[end].add(begin)
    remaining = set(changed)
    components = 0
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
    return components


def reaction_template_metrics(reaction_smarts: str) -> ReactionTemplateMetrics:
    Chem, AllChem, *_rest = require_rdkit()
    try:
        reaction = AllChem.ReactionFromSmarts(str(reaction_smarts))
        if reaction is None:
            raise ValueError("ReactionFromSmarts returned None")
        reaction.Initialize()
        n_reactants = int(reaction.GetNumReactantTemplates())
        n_products = int(reaction.GetNumProductTemplates())
        reactant_atoms = 0
        product_atoms = 0
        mapped_reactant_atoms = 0
        mapped_product_atoms = 0
        generic_reactant_atoms = 0
        fp_hex = ""
        if n_reactants:
            reactant = reaction.GetReactantTemplate(0)
            reactant_atoms = int(reactant.GetNumAtoms())
            mapped_reactant_atoms = sum(
                1 for atom in reactant.GetAtoms() if atom.GetAtomMapNum() > 0
            )
            generic_reactant_atoms = sum(
                1 for atom in reactant.GetAtoms() if atom.GetAtomicNum() == 0
            )
            fp = Chem.PatternFingerprint(reactant, fpSize=2048)
            fp_hex = fp.ToBitString()
        if n_products:
            product = reaction.GetProductTemplate(0)
            product_atoms = int(product.GetNumAtoms())
            mapped_product_atoms = sum(
                1 for atom in product.GetAtoms() if atom.GetAtomMapNum() > 0
            )
        changed = (
            _changed_mapped_atoms(reactant, product)
            if n_reactants == 1 and n_products == 1
            else set()
        )
        center_count = (
            _independent_center_count(reactant, product, changed)
            if n_reactants == 1 and n_products == 1
            else 0
        )
        structural_delta = (
            _element_delta_text(reactant, product)
            if n_reactants == 1 and n_products == 1
            else ""
        )
        signature = (
            _reaction_edit_signature(reactant, product, changed)
            if n_reactants == 1 and n_products == 1 and changed
            else ""
        )
        reactant_map_numbers = {
            atom.GetAtomMapNum()
            for atom in reactant.GetAtoms()
            if atom.GetAtomMapNum() > 0
        } if n_reactants == 1 else set()
        product_map_numbers = {
            atom.GetAtomMapNum()
            for atom in product.GetAtoms()
            if atom.GetAtomMapNum() > 0
        } if n_products == 1 else set()
        retention = (
            len(reactant_map_numbers & product_map_numbers)
            / len(reactant_map_numbers)
            if reactant_map_numbers
            else 0.0
        )
        mapping_coverage = (
            len(reactant_map_numbers) / reactant_atoms if reactant_atoms else 0.0
        )
        return ReactionTemplateMetrics(
            compile_status="ok",
            compile_error="",
            n_reactants=n_reactants,
            n_products=n_products,
            reactant_atoms=reactant_atoms,
            product_atoms=product_atoms,
            mapped_reactant_atoms=mapped_reactant_atoms,
            mapped_product_atoms=mapped_product_atoms,
            generic_reactant_atoms=generic_reactant_atoms,
            reactant_pattern_fp_hex=fp_hex,
            changed_mapped_atom_count=len(changed),
            inferred_independent_reaction_centers=center_count,
            structural_element_delta=structural_delta,
            reaction_edit_signature=signature,
            mapped_atom_retention=round(retention, 6),
            mapping_coverage=round(mapping_coverage, 6),
        )
    except Exception as exc:
        return ReactionTemplateMetrics(
            compile_status="failed",
            compile_error=str(exc),
            n_reactants=0,
            n_products=0,
            reactant_atoms=0,
            product_atoms=0,
            mapped_reactant_atoms=0,
            mapped_product_atoms=0,
            generic_reactant_atoms=0,
            reactant_pattern_fp_hex="",
            changed_mapped_atom_count=0,
            inferred_independent_reaction_centers=0,
            structural_element_delta="",
            reaction_edit_signature="",
            mapped_atom_retention=0.0,
            mapping_coverage=0.0,
        )
