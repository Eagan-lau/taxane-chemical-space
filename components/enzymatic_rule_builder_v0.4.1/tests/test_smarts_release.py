from pathlib import Path
import pandas as pd

from enzymatic_rule_builder.pipeline import run_build_all
from enzymatic_rule_builder.release import export_reaction_smarts_library


def test_taxol_only_without_abstraction_emits_no_reaction_smarts_rules(tmp_path: Path):
    summary = run_build_all(
        manifest=Path("configs/source_manifest.example.yaml"),
        output_root=tmp_path / "taxol_only",
        max_decoys=5,
        abstract_exact_reactions=False,
    )
    core = pd.read_csv(tmp_path / "taxol_only" / "04_release" / "reaction_smarts_library.T1_core.tsv", sep="\t")
    assert len(core) == 0
    assert summary["reaction_smarts_library_core"] == 0
    anchors = pd.read_csv(tmp_path / "taxol_only" / "04_release" / "curated_exact_anchor_edges.tsv", sep="\t")
    assert len(anchors) > 0


def test_demo_generalized_emits_network_ready_reaction_smarts(tmp_path: Path):
    summary = run_build_all(
        manifest=Path("examples/demo_generalized/manifest.yaml"),
        output_root=tmp_path / "demo",
        max_decoys=1,
    )
    core = pd.read_csv(tmp_path / "demo" / "04_release" / "reaction_smarts_library.T1_core.tsv", sep="\t")
    assert len(core) == 1
    assert summary["reaction_smarts_library_core"] == 1
    assert core.loc[0, "reaction_smarts"]
    assert core.loc[0, "template_scope"] == "generalized_template"
    assert core.loc[0, "template_qc_status"] == "ok"


def test_require_core_rules_fails_when_all_smarts_but_core_empty(tmp_path: Path):
    import pytest
    import yaml

    source = tmp_path / "templates.tsv"
    pd.DataFrame([
        {
            "source_reaction_id": "R_EXPANDED_ONLY",
            "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
            "reaction_type_source": "alcohol_oxidation",
            "direction": "forward",
        }
    ]).to_csv(source, sep="\t", index=False)
    manifest = tmp_path / "manifest.yaml"
    manifest.write_text(yaml.safe_dump({
        "datasets": {
            "demo": {
                "enabled": True,
                "path": str(source),
                "parser": "generic_reaction_table",
                "delimiter": "\t",
                "source_database": "BioNaviNP_BioChem",
                "evidence_layer": "T2_Bio_Extended",
            }
        }
    }), encoding="utf-8")

    summary = run_build_all(manifest=manifest, output_root=tmp_path / "ok", max_decoys=0)
    assert summary["reaction_smarts_rules_all"] == 1
    assert summary["reaction_smarts_rules_core"] == 0
    assert summary["core_release_empty_warning"] is True

    with pytest.raises(ValueError, match="No strict-core predictive reaction SMARTS rules"):
        run_build_all(manifest=manifest, output_root=tmp_path / "fail", max_decoys=0, require_core_rules=True)


def test_t1_t2_t3_releases_are_mutually_exclusive_by_priority():
    rules = pd.DataFrame([
        {
            "rule_id": "R_T1",
            "template_scope": "generalized_template",
            "predictive_rule_use": "true",
            "template_qc_status": "ok",
            "reaction_smarts": "[C:1]>>[C:1]O",
            "strict_core_use": "true",
            "expanded_use": "true",
            "exploratory_use": "true",
        },
        {
            "rule_id": "R_T2",
            "template_scope": "generalized_template",
            "predictive_rule_use": "true",
            "template_qc_status": "ok",
            "reaction_smarts": "[C:1]O>>[C:1]=O",
            "strict_core_use": "false",
            "expanded_use": "true",
            "exploratory_use": "true",
        },
        {
            "rule_id": "R_T3",
            "template_scope": "generalized_template",
            "predictive_rule_use": "true",
            "template_qc_status": "ok",
            "reaction_smarts": "[C:1]>>[C:1]C",
            "strict_core_use": "false",
            "expanded_use": "false",
            "exploratory_use": "true",
        },
    ])
    t1 = export_reaction_smarts_library(rules, tier="T1_core")
    t2 = export_reaction_smarts_library(rules, tier="T2_extended")
    t3 = export_reaction_smarts_library(rules, tier="T3_exploratory")
    all_smarts = export_reaction_smarts_library(rules, tier="all")
    t2_cumulative = export_reaction_smarts_library(rules, tier="T2_cumulative")
    t3_cumulative = export_reaction_smarts_library(rules, tier="T3_cumulative")

    assert set(t1["rule_id"]) == {"R_T1"}
    assert set(t2["rule_id"]) == {"R_T2"}
    assert set(t3["rule_id"]) == {"R_T3"}
    assert set(t1["rule_id"]).isdisjoint(set(t2["rule_id"]))
    assert set(t1["rule_id"]).isdisjoint(set(t3["rule_id"]))
    assert set(t2["rule_id"]).isdisjoint(set(t3["rule_id"]))
    assert set(all_smarts["exclusive_release_tier"]) == {"T1_only", "T2_only", "T3_only"}
    assert set(t2_cumulative["rule_id"]) == {"R_T1", "R_T2"}
    assert set(t3_cumulative["rule_id"]) == {"R_T1", "R_T2", "R_T3"}
