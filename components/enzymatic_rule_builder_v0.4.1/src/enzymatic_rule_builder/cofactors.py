from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

from .chem import build_reaction_smiles, canonical_smiles, heavy_atom_count, parse_reaction_smiles
from .direction import direction_qc_from_handling, source_direction_mode
from .utils import clean_text, join_values, read_yaml, split_multi_value


@lru_cache(maxsize=128)
def _smarts_query(smarts: str):
    try:
        from rdkit import Chem  # type: ignore
    except Exception:  # pragma: no cover - optional dependency
        return None
    try:
        return Chem.MolFromSmarts(smarts)
    except Exception:
        return None


@lru_cache(maxsize=200000)
def _mol_has_substructure(smiles: str, smarts: str) -> bool:
    try:
        from rdkit import Chem  # type: ignore
    except Exception:  # pragma: no cover - optional dependency
        return False
    query = _smarts_query(smarts)
    if query is None:
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
        return bool(mol and mol.HasSubstructMatch(query))
    except Exception:
        return False


def _looks_like_coa_carrier(smiles: str) -> bool:
    s = clean_text(smiles)
    if not s:
        return False
    # CoA and acyl-CoA are large, sulfur/phosphate-containing carriers.  This
    # deliberately avoids silently deleting every sulfur compound; it only marks
    # CoA-like external participants for main-pair projection.
    return ("S" in s and "P" in s and heavy_atom_count(s) >= 25)


def _looks_like_acyl_coa_donor(smiles: str) -> bool:
    s = clean_text(smiles)
    if not s or "S" not in s:
        return False
    return _looks_like_coa_carrier(s) and _mol_has_substructure(s, "[#6](=[#8])-[#16]")


def _classify_acyl_donor(smiles: str) -> tuple[str, str]:
    """Return donor class and transferred group for common acyl-CoA donors."""
    s = clean_text(smiles)
    if _mol_has_substructure(s, "[#6](=[#8])-[#6]=[#6]-[c]") or _mol_has_substructure(
        s, "[c]-[#6]=[#6]-[#6](=[#8])-[#16]"
    ):
        return "cinnamoyl-CoA_or_hydroxycinnamoyl-CoA", "cinnamoyl_or_hydroxycinnamoyl"
    if _mol_has_substructure(s, "[c]-[#6](=[#8])-[#16]"):
        return "benzoyl-CoA_or_aromatic_acyl-CoA", "benzoyl_or_aromatic_acyl"
    if _mol_has_substructure(s, "[#6H3]-[#6](=[#8])-[#16]"):
        return "acetyl-CoA", "acetyl"
    if _mol_has_substructure(s, "[#6](=[#8])-[#16]"):
        return "acyl-CoA", "acyl"
    return "external_acyl_donor", "acyl"


def _phosphorus_count(smiles: str) -> int:
    return clean_text(smiles).count("P")


def _sulfur_count(smiles: str) -> int:
    return clean_text(smiles).count("S")


def _looks_like_nucleotide_carrier(smiles: str) -> bool:
    s = clean_text(smiles)
    if not s:
        return False
    return _phosphorus_count(s) >= 2 and heavy_atom_count(s) >= 20 and _sulfur_count(s) == 0


def _looks_like_sam_or_sah(smiles: str) -> bool:
    s = clean_text(smiles)
    if not s:
        return False
    return _sulfur_count(s) >= 1 and _phosphorus_count(s) == 0 and 15 <= heavy_atom_count(s) <= 50


def _classify_transfer_donor(smiles: str) -> dict[str, str] | None:
    if _looks_like_acyl_coa_donor(smiles):
        donor_class, transferred_group = _classify_acyl_donor(smiles)
        return {
            "donor_class": donor_class,
            "transferred_group": transferred_group,
            "transferred_group_class": transferred_group,
            "leaving_group_class": "CoA",
            "right_external_product": "CoA",
            "projection_method": "role_aware_acyl_CoA_acceptor_projection",
        }
    if _looks_like_nucleotide_carrier(smiles):
        p_count = _phosphorus_count(smiles)
        if p_count >= 3 and heavy_atom_count(smiles) < 45:
            return {
                "donor_class": "ATP_or_nucleotide_triphosphate",
                "transferred_group": "phosphoryl",
                "transferred_group_class": "phosphoryl",
                "leaving_group_class": "ADP_or_nucleotide_diphosphate",
                "right_external_product": "ADP_or_nucleotide_diphosphate",
                "projection_method": "role_aware_nucleotide_phosphoryl_donor_projection",
            }
        return {
            "donor_class": "UDP_sugar_or_nucleotide_sugar",
            "transferred_group": "glycosyl",
            "transferred_group_class": "glycosyl",
            "leaving_group_class": "UDP_or_nucleotide_diphosphate",
            "right_external_product": "UDP_or_nucleotide_diphosphate",
            "projection_method": "role_aware_nucleotide_sugar_acceptor_projection",
        }
    if _looks_like_sam_or_sah(smiles):
        return {
            "donor_class": "SAM_or_sulfur_methyl_donor",
            "transferred_group": "methyl",
            "transferred_group_class": "methyl",
            "leaving_group_class": "SAH_or_sulfur_carrier",
            "right_external_product": "SAH_or_sulfur_carrier",
            "projection_method": "role_aware_SAM_methyl_acceptor_projection",
        }
    return None


def _looks_like_external_carrier(smiles: str) -> bool:
    return (
        _looks_like_coa_carrier(smiles)
        or _looks_like_nucleotide_carrier(smiles)
        or _looks_like_sam_or_sah(smiles)
        or heavy_atom_count(smiles) <= 3
    )


def _acceptor_atom_class(smiles: str) -> str:
    s = clean_text(smiles)
    if _mol_has_substructure(s, "[N;H1,H2,H3]"):
        return "amine_N_acceptor"
    if _mol_has_substructure(s, "[O;H1]-[#6]"):
        return "alcohol_O_acceptor"
    if _mol_has_substructure(s, "[S;H1]"):
        return "thiol_S_acceptor"
    return "unclassified_acceptor"


def _project_transfer_donor_main_pair(
    reactants_can: list[str],
    products_can: list[str],
    registry_entries: dict[str, dict[str, str]],
) -> dict[str, str] | None:
    donor_specs = [(x, _classify_transfer_donor(x)) for x in reactants_can]
    donor_specs = [(x, spec) for x, spec in donor_specs if spec]
    donors = [x for x, _ in donor_specs]
    if not donors:
        return None
    donor, donor_spec = max(donor_specs, key=lambda item: heavy_atom_count(item[0]))
    assert donor_spec is not None
    donor_class = donor_spec["donor_class"]
    transferred_group = donor_spec["transferred_group"]
    reactant_aux = set(donors)
    product_aux = {x for x in products_can if _looks_like_external_carrier(x) or x in registry_entries}
    reactant_aux.update(x for x in reactants_can if x in registry_entries and x != donor)

    acceptor_candidates = [
        x for x in reactants_can
        if x not in reactant_aux and not _looks_like_external_carrier(x) and heavy_atom_count(x) > 1
    ]
    product_candidates = [
        x for x in products_can
        if x not in product_aux and not _looks_like_external_carrier(x) and heavy_atom_count(x) > 1
    ]
    if not acceptor_candidates or not product_candidates:
        return None

    product = max(product_candidates, key=heavy_atom_count)

    def acceptor_score(candidate: str) -> tuple[int, int, int]:
        cls = _acceptor_atom_class(candidate)
        acceptor_bonus = 2 if cls == "amine_N_acceptor" else (1 if cls == "alcohol_O_acceptor" else 0)
        diff = heavy_atom_count(product) - heavy_atom_count(candidate)
        plausible_addition = 1 if 1 <= diff <= 30 else 0
        return (acceptor_bonus, plausible_addition, heavy_atom_count(candidate))

    acceptor = max(acceptor_candidates, key=acceptor_score)
    if heavy_atom_count(product) <= heavy_atom_count(acceptor):
        return None

    return {
        "main_substrate_smiles": acceptor,
        "main_product_smiles": product,
        "donor_smiles": donor,
        "removed_participants": join_values([x for x in reactants_can + products_can if x not in {acceptor, product}]),
        "external_participant_roles": join_values([
            f"{donor_class}:left_required_external_participant",
            f"{donor_spec['right_external_product']}:right_external_product" if product_aux else "",
        ]),
        "cofactor_or_donor_class": donor_class,
        "donor_class": donor_class,
        "acceptor_atom_class": _acceptor_atom_class(acceptor),
        "transferred_group": transferred_group,
        "transferred_group_class": donor_spec["transferred_group_class"],
        "leaving_group_class": donor_spec["leaving_group_class"],
        "main_pair_projection_method": donor_spec["projection_method"],
        "main_pair_projection_note": (
            "External transfer donor was treated as an external participant; main pair is "
            "the acceptor scaffold to the transferred-group product."
        ),
    }


def _direction_fields(direction_handling: str, direction_variant: str) -> dict[str, str]:
    qc, note = direction_qc_from_handling(direction_handling, direction_variant)
    return {
        "normalized_direction": "substrate_to_product",
        "direction_qc_status": qc,
        "direction_qc_note": note,
    }


class CofactorRegistry:
    """External participant registry for donor/cofactor/carrier handling.

    The registry is deliberately loaded from an auditable YAML file rather than
    hard-coded in Python. The recommended production workflow is to build this
    YAML from databases with `enzymatic-rules build-participant-registry`, review
    it, then pass it to `build-rules` via `--cofactor-yaml`.
    """

    def __init__(self, entries: dict[str, dict[str, str] | str] | None = None) -> None:
        norm: dict[str, dict[str, str]] = {}
        for k, v in (entries or {}).items():
            if isinstance(v, dict):
                vv = dict(v)
                vv.setdefault("label", vv.get("class") or vv.get("role_class") or "external_participant")
                vv.setdefault("role_class", vv.get("role") or vv.get("label") or "external_participant")
                vv.setdefault("confidence", "0.80")
                norm[k] = vv
            else:
                norm[k] = {"label": str(v), "role_class": str(v), "confidence": "0.80", "transferred_group": "", "leaving_group_class": ""}
        self.entries = norm

    @classmethod
    def from_yaml(cls, path: str | Path | None) -> "CofactorRegistry":
        if not path:
            return cls({})
        data = read_yaml(path)
        entries: dict[str, dict[str, str]] = {}

        def add_entry(label: str, role: str, smi: str, row: dict[str, Any] | None = None) -> None:
            can = canonical_smiles(smi)
            if not can:
                return
            r = row or {}
            entries[can] = {
                "class": clean_text(label) or "external_participant",
                "role": clean_text(role) or clean_text(r.get("role_class", "")) or "external_participant",
                "confidence": clean_text(r.get("confidence", r.get("role_confidence", ""))) or "",
                "transferred_group": clean_text(r.get("transferred_group", "")),
                "leaving_group_class": clean_text(r.get("leaving_group", r.get("leaving_group_class", ""))),
                "provenance": clean_text(r.get("provenance", "")),
            }

        for row in data.get("participants", []) or []:
            label = clean_text(row.get("class") or row.get("registry_class") or row.get("name") or "external_participant")
            role = clean_text(row.get("role") or row.get("role_class") or "external_participant")
            for smi in row.get("smiles", []) or []:
                add_entry(label, role, smi, row)
        for label0, row in (data.get("classes", {}) or {}).items():
            if not isinstance(row, dict):
                continue
            label = clean_text(row.get("class")) or clean_text(label0) or "external_participant"
            role = clean_text(row.get("role")) or clean_text(row.get("role_class")) or "external_participant"
            for smi in row.get("smiles", []) or []:
                add_entry(label, role, smi, row)
        return cls(entries)

    def labels_for(self, components: list[str]) -> str:
        labels = []
        for comp in components:
            can = canonical_smiles(comp)
            if can in self.entries:
                labels.append(self.entries[can].get("label", "external_participant"))
        return join_values(labels)

    def roles_for(self, components: list[str]) -> str:
        roles = []
        for comp in components:
            can = canonical_smiles(comp)
            if can in self.entries:
                e = self.entries[can]
                roles.append(f"{e.get('class','external_participant')}:{e.get('role','external_participant')}")
        return join_values(roles)

    def transferred_groups_for(self, components: list[str]) -> str:
        vals = []
        for comp in components:
            can = canonical_smiles(comp)
            if can in self.entries and self.entries[can].get("transferred_group"):
                vals.append(self.entries[can]["transferred_group"])
        return join_values(vals)

    def leaving_groups_for(self, components: list[str]) -> str:
        vals = []
        for comp in components:
            can = canonical_smiles(comp)
            if can in self.entries and self.entries[can].get("leaving_group_class"):
                vals.append(self.entries[can]["leaving_group_class"])
        return join_values(vals)

    def _split_core_aux(self, components: list[str]) -> tuple[list[str], list[str]]:
        core: list[str] = []
        aux: list[str] = []
        for comp in components:
            if comp in self.entries:
                aux.append(comp)
            else:
                core.append(comp)
        return core, aux

    def strip_reaction(
        self,
        reaction_smiles: str,
        substrate_smiles: str = "",
        product_smiles: str = "",
        source_cofactor_class: str = "",
        source_direction: str = "",
        is_reversible: str = "",
    ) -> dict[str, str]:
        sub = canonical_smiles(substrate_smiles)
        prod = canonical_smiles(product_smiles)
        direction_mode = source_direction_mode(source_direction, is_reversible)
        # Explicit substrate/product columns are curated main-pair fields.
        if sub and prod:
            delta_rxn = build_reaction_smiles([sub], [prod])
            return {
                "main_substrate_smiles": sub,
                "main_product_smiles": prod,
                "canonical_reaction_smiles": delta_rxn,
                "canonical_substrate_smiles": sub,
                "canonical_product_smiles": prod,
                "removed_participants": "",
                "external_participant_roles": "",
                "participant_role_confidence": "",
                "reaction_representation_scope": "source_curated_node_transformation",
                "transferred_group": "",
                "leaving_group_class": "",
                "cofactor_or_donor_class": clean_text(source_cofactor_class),
                "main_pair_method": "source_substrate_product_columns",
                "main_pair_confidence": "1.0",
                "direction_handling": "source_columns_as_curated_main_pair",
                "reversible_group_id": "",
                "direction_variant": "forward",
                **_direction_fields("source_columns_as_curated_main_pair", "forward"),
            }
        reactants, products = parse_reaction_smiles(reaction_smiles)
        if not reactants or not products:
            return {
                "main_substrate_smiles": "", "main_product_smiles": "", "canonical_reaction_smiles": "",
                "canonical_substrate_smiles": "", "canonical_product_smiles": "", "removed_participants": "",
                "external_participant_roles": "", "participant_role_confidence": "", "reaction_representation_scope": "unknown",
                "transferred_group": "", "leaving_group_class": "",
                "cofactor_or_donor_class": clean_text(source_cofactor_class), "main_pair_method": "missing", "main_pair_confidence": "0.0",
                "direction_handling": "missing", "reversible_group_id": "", "direction_variant": "",
                **_direction_fields("missing", ""),
            }
        reactants_can = [canonical_smiles(x) for x in reactants]
        products_can = [canonical_smiles(x) for x in products]
        reactants_can = [x for x in reactants_can if x]
        products_can = [x for x in products_can if x]
        if not reactants_can or not products_can:
            return {
                "main_substrate_smiles": "", "main_product_smiles": "", "canonical_reaction_smiles": "",
                "canonical_substrate_smiles": "", "canonical_product_smiles": "", "removed_participants": "",
                "external_participant_roles": "", "participant_role_confidence": "", "reaction_representation_scope": "unknown",
                "transferred_group": "", "leaving_group_class": "",
                "cofactor_or_donor_class": clean_text(source_cofactor_class), "main_pair_method": "parse_failed", "main_pair_confidence": "0.0",
                "direction_handling": "parse_failed", "reversible_group_id": "", "direction_variant": "",
                **_direction_fields("parse_failed", ""),
            }

        role_projection = _project_transfer_donor_main_pair(reactants_can, products_can, self.entries)
        if role_projection:
            main_sub = role_projection["main_substrate_smiles"]
            main_prod = role_projection["main_product_smiles"]
            if direction_mode == "source_reverse":
                main_sub, main_prod = main_prod, main_sub
                direction_handling = "reversed_from_source"
                direction_variant = "forward_after_source_reverse_correction"
            elif direction_mode == "source_reversible":
                direction_handling = "source_reversible_main_pair_left_to_right"
                direction_variant = "forward_member"
            elif direction_mode == "source_forward":
                direction_handling = "kept_forward"
                direction_variant = "forward"
            else:
                direction_handling = "unknown_direction_kept_left_to_right"
                direction_variant = "left_to_right"
            rxn = build_reaction_smiles([main_sub], [main_prod])
            labels = split_multi_value(source_cofactor_class) + split_multi_value(role_projection["cofactor_or_donor_class"])
            return {
                "main_substrate_smiles": main_sub,
                "main_product_smiles": main_prod,
                "canonical_reaction_smiles": rxn,
                "canonical_substrate_smiles": main_sub,
                "canonical_product_smiles": main_prod,
                "removed_participants": role_projection["removed_participants"],
                "external_participant_roles": role_projection["external_participant_roles"],
                "participant_role_confidence": "0.85",
                "reaction_representation_scope": "role_aware_main_pair_projection_with_external_participant_annotation",
                "transferred_group": role_projection["transferred_group"],
                "leaving_group_class": role_projection["leaving_group_class"],
                "cofactor_or_donor_class": join_values(labels),
                "donor_class": role_projection["donor_class"],
                "acceptor_atom_class": role_projection["acceptor_atom_class"],
                "transferred_group_class": role_projection["transferred_group_class"],
                "main_pair_projection_method": role_projection["main_pair_projection_method"],
                "main_pair_projection_note": role_projection["main_pair_projection_note"],
                "main_pair_method": role_projection["main_pair_projection_method"],
                "main_pair_confidence": "0.85",
                "direction_handling": direction_handling,
                "reversible_group_id": "",
                "direction_variant": direction_variant,
                **_direction_fields(direction_handling, direction_variant),
            }

        reactant_core, reactant_aux = self._split_core_aux(reactants_can)
        product_core, product_aux = self._split_core_aux(products_can)
        used_registry = bool(reactant_aux or product_aux)
        candidate_reactants = reactant_core or reactants_can
        candidate_products = product_core or products_can

        main_sub = max(candidate_reactants, key=heavy_atom_count)
        main_prod = max(candidate_products, key=heavy_atom_count)

        if direction_mode == "source_reverse":
            main_sub, main_prod = main_prod, main_sub
            direction_handling = "reversed_from_source"
            direction_variant = "forward_after_source_reverse_correction"
        elif direction_mode == "source_reversible":
            direction_handling = "source_reversible_main_pair_left_to_right"
            direction_variant = "forward_member"
        elif direction_mode == "source_forward":
            direction_handling = "kept_forward"
            direction_variant = "forward"
        else:
            direction_handling = "unknown_direction_kept_left_to_right"
            direction_variant = "left_to_right"

        removed = [c for c in reactants_can + products_can if c not in {main_sub, main_prod}]
        aux = reactant_aux + product_aux
        labels = split_multi_value(source_cofactor_class) + split_multi_value(self.labels_for(aux))
        roles = self.roles_for(aux)
        rxn = build_reaction_smiles([main_sub], [main_prod])
        method = "registry_stripped_largest_core_pair" if used_registry else "largest_component_pair"
        if direction_mode == "source_reverse":
            method += "_source_reverse_corrected"
        confidence = "0.90" if used_registry and reactant_core and product_core else ("0.65" if removed else "0.85")
        return {
            "main_substrate_smiles": main_sub,
            "main_product_smiles": main_prod,
            "canonical_reaction_smiles": rxn,
            "canonical_substrate_smiles": main_sub,
            "canonical_product_smiles": main_prod,
            "removed_participants": join_values(removed),
            "external_participant_roles": roles,
            "participant_role_confidence": "0.90" if used_registry else "",
            "reaction_representation_scope": "node_transformation_with_external_participant_annotation" if used_registry else "node_transformation_largest_component_pair",
            "transferred_group": self.transferred_groups_for(aux),
            "leaving_group_class": self.leaving_groups_for(aux),
            "cofactor_or_donor_class": join_values(labels),
            "main_pair_method": method,
            "main_pair_confidence": confidence,
            "direction_handling": direction_handling,
            "reversible_group_id": "",
            "direction_variant": direction_variant,
            **_direction_fields(direction_handling, direction_variant),
        }
