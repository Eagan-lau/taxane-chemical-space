from pathlib import Path

import pandas as pd
import yaml

from enzymatic_rule_builder.cofactors import CofactorRegistry
from enzymatic_rule_builder.family import annotate_families


def test_cofactor_registry_accepts_classes_schema(tmp_path):
    p = tmp_path / "cofactors.yaml"
    p.write_text(yaml.safe_dump({"classes": {"water": {"smiles": ["O"]}}}))
    reg = CofactorRegistry.from_yaml(p)
    assert reg.labels_for(["O"]) == "water"


def test_family_evidence_external_match_only():
    fe = pd.DataFrame([
        {"match_type": "ec_prefix", "match_value": "1.14", "primary_family": "CYP450", "confidence": "0.9", "evidence_source": "external_demo"}
    ])
    primary, secondary, mode, evidence = annotate_families({"candidate_ec_numbers": "1.14.14.-"}, fe)
    assert primary == "CYP450"
    assert secondary == ""
    assert mode == "external_evidence"
    assert "external_demo" in evidence


def test_family_ec_prefix_does_not_overmatch():
    fe = pd.DataFrame([
        {"match_type": "ec_prefix", "match_value": "1.1", "primary_family": "ADH", "confidence": "0.9", "evidence_source": "external_demo"}
    ])
    primary, secondary, mode, evidence = annotate_families({"candidate_ec_numbers": "1.14.14.176"}, fe)
    assert primary == ""
    assert mode == "none"


def test_cofactor_registry_strips_registered_donor_before_largest_pair():
    from enzymatic_rule_builder.chem import canonical_smiles

    substrate = "CCO"
    product = "CCOC(C)=O"
    donor = "CCCCCCCCCCCCCCCC(=O)SC"
    coproduct = "CCCCCCCCCCCCCCCS"
    rxn = f"{substrate}.{donor}>>{product}.{coproduct}"
    reg = CofactorRegistry({
        canonical_smiles(donor): "acyl_donor",
        canonical_smiles(coproduct): "thiol_coproduct",
    })
    out = reg.strip_reaction(rxn)
    assert out["main_substrate_smiles"] == canonical_smiles(substrate)
    assert out["main_product_smiles"] == canonical_smiles(product)
    assert out["main_pair_method"] == "registry_stripped_largest_core_pair"
    assert "acyl_donor" in out["cofactor_or_donor_class"]
