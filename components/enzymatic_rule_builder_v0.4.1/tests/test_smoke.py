from pathlib import Path

import pandas as pd

from enzymatic_rule_builder.pipeline import run_build_all


def test_taxol_only_without_abstraction_outputs_anchor_but_no_smarts_release(tmp_path):
    root = Path(__file__).resolve().parents[1]
    summary = run_build_all(
        manifest=root / "configs" / "source_manifest.example.yaml",
        output_root=tmp_path / "build",
        max_decoys=5,
        abstract_exact_reactions=False,
    )
    assert summary["source_records"] > 0
    assert summary["exact_anchor_rules"] > 0
    assert summary["predictive_generalized_rules"] == 0
    assert summary["reaction_smarts_rules_all"] == 0

    audit = pd.read_csv(tmp_path / "build" / "03_rules" / "general_transformation_rules.audit.tsv", sep="\t")
    assert set(audit["template_scope"]) == {"exact_anchor"}

    smarts_all = pd.read_csv(tmp_path / "build" / "04_release" / "reaction_smarts_rules.all.tsv", sep="\t")
    assert len(smarts_all) == 0

    anchors = pd.read_csv(tmp_path / "build" / "04_release" / "curated_exact_anchor_edges.tsv", sep="\t")
    assert len(anchors) > 0


def test_generalized_demo_produces_smarts_only_release(tmp_path):
    root = Path(__file__).resolve().parents[1]
    demo = root / "examples" / "demo_generalized"
    summary = run_build_all(
        manifest=demo / "manifest.yaml",
        output_root=tmp_path / "build_demo",
        max_decoys=5,
        require_smarts_rules=True,
    )
    assert summary["reaction_smarts_rules_all"] == 1
    assert summary["reaction_smarts_rules_core"] == 1
    smarts = pd.read_csv(tmp_path / "build_demo" / "04_release" / "reaction_smarts_rules.core.tsv", sep="\t")
    assert len(smarts) == 1
    assert smarts["reaction_smarts"].astype(str).str.len().min() > 0
    assert set(smarts["template_scope"]) == {"generalized_template"}
    assert smarts["predictive_rule_use"].astype(str).str.lower().eq("true").all()
