from pathlib import Path

import pandas as pd
import pytest
import yaml

from enzymatic_rule_builder.cofactors import CofactorRegistry
from enzymatic_rule_builder.normalization import add_main_pairs
from enzymatic_rule_builder.pipeline import run_build_all
from enzymatic_rule_builder.rules import build_rule_library
from enzymatic_rule_builder.templates import build_templates, qc_templates, deduplicate_templates


def _one_source(**overrides):
    row = {
        "record_id": "SRC:1",
        "source_database": "RetroRules",
        "evidence_layer": "T1_Bio_Core",
        "source_reaction_id": "R1",
        "source_file": "demo",
        "parser_name": "test",
        "reaction_smiles": "CCO>>CC=O",
        "substrate_smiles": "",
        "product_smiles": "",
        "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
        "reaction_equation": "",
        "direction": "forward",
        "is_reversible": "false",
        "ec_numbers": "1.1.1.1",
        "template_ec_candidates": "1.1.1.1",
        "database_ec_candidates": "",
        "ec_prior_candidates": "",
        "rhea_ids": "",
        "kegg_ids": "",
        "metanetx_ids": "",
        "reaction_type_source": "alcohol_oxidation",
        "reaction_subtype_source": "",
        "cofactor_or_donor_class": "",
        "enzyme_name": "",
        "protein_ids": "",
        "protein_sequence": "",
        "source_evidence_text": "",
        "raw_row_index": "0",
        "source_row_hash": "h",
    }
    row.update(overrides)
    return pd.DataFrame([row])


def test_reverse_source_direction_reverses_smarts_and_main_pair():
    # Source stores B -> A but says that orientation is reverse relative to the
    # builder's substrate-to-product convention. Both exact pair and SMARTS are
    # corrected before release.
    src = _one_source(
        reaction_smiles="CC=O>>CCO",
        reaction_smarts="[C:1]=[O:2]>>[C:1][OH:2]",
        direction="reverse",
    )
    main = add_main_pairs(src)
    assert main.loc[0, "canonical_reaction_smiles"] == "CCO>>CC=O"
    assert main.loc[0, "direction_handling"] == "reversed_from_source"
    templates = deduplicate_templates(qc_templates(build_templates(main)))
    rules = build_rule_library(templates, pd.DataFrame())
    assert len(rules) == 1
    assert rules.loc[0, "reaction_smarts"] == "[C:1][OH:2]>>[C:1]=[O:2]"
    assert rules.loc[0, "direction_handling"] == "reversed_from_source"
    assert rules.loc[0, "direction_evidence_type"] == "source_reverse_corrected"


def test_reversible_source_splits_into_two_directional_rules():
    src = _one_source(direction="reversible", is_reversible="true")
    main = add_main_pairs(src)
    templates = deduplicate_templates(qc_templates(build_templates(main)))
    rules = build_rule_library(templates, pd.DataFrame())
    assert set(rules["direction_handling"]) == {"split_reversible_forward", "split_reversible_reverse"}
    assert set(rules["reaction_smarts"]) == {
        "[C:1][OH:2]>>[C:1]=[O:2]",
        "[C:1]=[O:2]>>[C:1][OH:2]",
    }
    assert rules["reverse_transform_available"].astype(str).str.lower().eq("true").all()
    reverse_row = rules[rules["direction_handling"] == "split_reversible_reverse"].iloc[0]
    assert reverse_row["example_substrate_smiles"] == "CC=O"
    assert reverse_row["example_product_smiles"] == "CCO"


def test_cofactor_registry_removes_large_donor_before_largest_pair_selection(tmp_path: Path):
    donor = "CCCCCCCCCCCCCCCCCCCC(=O)SCC"
    coproduct = "CCCCCCCCCCCCCCCCCCCCS"
    substrate = "CCO"
    product = "CCOC(C)=O"
    rxn = f"{substrate}.{donor}>>{product}.{coproduct}"
    cofactor_yaml = tmp_path / "cofactors.yaml"
    cofactor_yaml.write_text(
        yaml.safe_dump({"classes": {"acyl_donor": {"smiles": [donor]}, "thiol_coproduct": {"smiles": [coproduct]}}}),
        encoding="utf-8",
    )
    reg = CofactorRegistry.from_yaml(cofactor_yaml)
    stripped = reg.strip_reaction(rxn)
    assert stripped["main_substrate_smiles"] == "CCO"
    assert stripped["main_product_smiles"] == "CCOC(C)=O"
    assert stripped["main_pair_method"] == "registry_stripped_largest_core_pair"
    assert "acyl_donor" in stripped["cofactor_or_donor_class"]
    assert "thiol_coproduct" in stripped["cofactor_or_donor_class"]


def test_require_core_rules_fails_when_all_smarts_exists_but_core_is_empty(tmp_path: Path):
    table = tmp_path / "templates.tsv"
    pd.DataFrame([
        {
            "source_reaction_id": "R1",
            "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
            "reaction_smiles": "CCO>>CC=O",
            "ec_numbers": "1.1.1.1",
            "direction": "forward",
            "reaction_type_source": "alcohol_oxidation",
        }
    ]).to_csv(table, sep="\t", index=False)
    manifest = tmp_path / "manifest.yaml"
    manifest.write_text(
        yaml.safe_dump({
            "datasets": {
                "demo": {
                    "enabled": True,
                    "parser": "generic_reaction_table",
                    "path": str(table),
                    "delimiter": "\t",
                    "source_database": "DemoT2",
                    "evidence_layer": "T2_Bio_Extended",
                }
            },
            "family_evidence": {"enabled": False, "path": ""},
        }),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="strict-core"):
        run_build_all(manifest=manifest, output_root=tmp_path / "out", require_core_rules=True)
    summary = yaml.safe_load((tmp_path / "out" / "build_summary.json").read_text())
    assert summary["reaction_smarts_rules_all"] > 0
    assert summary["reaction_smarts_rules_core"] == 0
    assert summary["core_release_empty_warning"] is True
