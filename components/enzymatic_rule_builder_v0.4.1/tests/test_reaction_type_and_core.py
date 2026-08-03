import pandas as pd

from enzymatic_rule_builder.reaction_type import classify_reaction_type
from enzymatic_rule_builder.rules import build_rule_library


def test_structural_delta_specific_acetylation():
    rt, subtype, mode = classify_reaction_type({
        "reaction_delta_json": '{"atom_delta":{"C":2,"H":2,"O":1},"heavy_atom_delta":3,"exact_mass_delta":42.010565,"ring_count_delta":0}',
        "reaction_delta_fingerprint": "C+2|H+2|O+1|HA+3|M+42.011",
    })
    assert rt == "acetylation_or_deacetylation_like_acyl_transfer"
    assert mode == "structural_delta_specific"
    assert "C+2" in subtype


def test_t1_strict_ec_rule_can_be_core_without_family_table():
    templates = pd.DataFrame([
        {
            "template_hash": "abc",
            "template_scope": "generalized_template",
            "predictive_rule_use": "true",
            "anchor_edge_use": "false",
            "template_qc_status": "ok",
            "template_qc_note": "",
            "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
            "canonical_reaction_smiles": "CO>>C=O",
            "main_substrate_smiles": "CO",
            "main_product_smiles": "C=O",
            "reaction_delta_fingerprint": "",
            "reaction_delta_json": "",
            "source_database": "RetroRules",
            "evidence_layer": "T1_Bio_Core",
            "source_reaction_id": "RetroRules:R1",
            "record_id": "R1",
            "ec_numbers": "1.1.1.1",
            "template_ec_candidates": "1.1.1.1",
            "database_ec_candidates": "",
            "ec_prior_candidates": "",
            "rhea_ids": "RHEA:12345",
            "kegg_ids": "",
            "metanetx_ids": "",
            "reaction_type_source": "alcohol_oxidation",
            "reaction_subtype_source": "",
            "direction": "forward",
            "cofactor_or_donor_class": "",
            "enzyme_name": "",
            "protein_ids": "",
            "source_evidence_text": "retrorules_sqlite_join",
            "source_file": "demo",
        }
    ])
    rules = build_rule_library(templates, pd.DataFrame())
    assert len(rules) == 1
    assert rules.loc[0, "family_assignment_mode"] == "ec_supported_family_unassigned"
    assert bool(rules.loc[0, "strict_core_use"])
