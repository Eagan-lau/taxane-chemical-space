from pathlib import Path

import pandas as pd

from enzymatic_rule_builder.granularity import annotate_biochemical_step_granularity
from enzymatic_rule_builder.participant_registry import build_participant_registry_from_sources
from enzymatic_rule_builder.raw import build_raw_source_manifest
from enzymatic_rule_builder.rules import build_rule_library


def test_single_smarts_application_and_possible_composite_annotation():
    row = {
        "template_scope": "generalized_template",
        "reaction_smarts": "[C:1]>>[C:1]O",
        "reaction_type_source": "hydroxylation+acetylation",
        "reaction_delta_fingerprint": "O+1|C+2|HA+3|M+58.000",
    }
    g = annotate_biochemical_step_granularity(row, "hydroxylation+acetylation", "")
    assert g["rule_application_unit"] == "single_smarts_application"
    assert g["biochemical_step_granularity"] == "possible_composite_step"
    assert g["composite_rule_flag"] == "true"


def test_family_evidence_absence_does_not_block_core_for_high_quality_rule():
    templates = pd.DataFrame([
        {
            "template_hash": "h_ec_broad",
            "template_scope": "generalized_template",
            "predictive_rule_use": "true",
            "anchor_edge_use": "false",
            "template_qc_status": "ok",
            "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
            "canonical_reaction_smiles": "CCO>>CC=O",
            "main_substrate_smiles": "CCO",
            "main_product_smiles": "CC=O",
            "reaction_delta_fingerprint": "H-2|M-2.016",
            "reaction_type_source": "",
            "reaction_subtype_source": "",
            "direction": "forward",
            "source_database": "Rhea",
            "evidence_layer": "T1_Bio_Core",
            "source_reaction_id": "RHEA:10001",
            "record_id": "REC1",
            "ec_numbers": "1.-.-.-",
            "template_ec_candidates": "",
            "database_ec_candidates": "",
            "ec_prior_candidates": "",
            "rhea_ids": "RHEA:10001",
            "kegg_ids": "",
            "metanetx_ids": "",
            "cofactor_or_donor_class": "",
            "enzyme_name": "",
            "protein_ids": "",
            "source_evidence_text": "high quality source EC broad annotation",
            "source_file": "rhea.tsv",
        }
    ])
    rules = build_rule_library(templates, family_evidence=pd.DataFrame())
    assert len(rules) == 1
    row = rules.iloc[0]
    assert str(row["strict_core_use"]).lower() == "true"
    assert row["family_annotation_scope"] == "not_assigned_at_rule_stage"


def test_participant_registry_is_database_derived_from_recurring_components():
    df = pd.DataFrame([
        {"reaction_smiles": "CCO.O>>CC=O.O", "source_database": "Rhea", "source_reaction_id": "R1"},
        {"reaction_smiles": "CCC.O>>CCCO.O", "source_database": "KEGG", "source_reaction_id": "R2"},
        {"reaction_smiles": "CCCC.O>>CCCCO.O", "source_database": "BioNavi", "source_reaction_id": "R3"},
    ])
    reg = build_participant_registry_from_sources(df, min_occurrence=3, max_heavy_atoms=5)
    assert len(reg) >= 1
    assert reg["role_class"].eq("small_molecule_reagent_or_byproduct").any()
    assert reg["registry_class"].eq("small_molecule").any()
    assert reg["provenance"].eq("database_derived_recurring_participant_registry").all()


def test_raw_manifest_taxol_anchors_are_auto_enabled_only_when_path_supplied(tmp_path: Path):
    root = tmp_path / "external"
    (root / "index").mkdir(parents=True)
    (root / "index" / "rhea_reaction_evidence.csv").write_text(
        "source_database,database_reaction_id,source_reaction_id,reaction_smiles,ec_numbers\n"
        "Rhea,10001,RHEA:10001,CCO>>CC=O,1.1.1.1\n",
        encoding="utf-8",
    )
    taxol = tmp_path / "taxol.csv"
    taxol.write_text("Enzyme,Substrate,Product,EC\nT,CCO,CC=O,1.1.-.-\n", encoding="utf-8")
    manifest1, disc1 = build_raw_source_manifest(
        external_db_root=root,
        output_dir=tmp_path / "m1",
        include_uspto=False,
        include_annotation_only=False,
        max_rhea=1,
    )
    assert disc1["taxol_pathway_include_requested"] is False
    manifest2, disc2 = build_raw_source_manifest(
        external_db_root=root,
        output_dir=tmp_path / "m2",
        taxol_pathway=taxol,
        include_uspto=False,
        include_annotation_only=False,
        max_rhea=1,
    )
    assert disc2["taxol_pathway_include_requested"] is True
    assert "taxol_known_pathway_anchors" in disc2["included_dataset_names"]
