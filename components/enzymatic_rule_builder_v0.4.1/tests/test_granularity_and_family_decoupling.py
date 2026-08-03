from pathlib import Path
import pandas as pd
import yaml

from enzymatic_rule_builder.pipeline import run_build_all
from enzymatic_rule_builder.granularity import annotate_biochemical_step_granularity


def test_core_rules_emit_single_smarts_application_and_granularity(tmp_path: Path):
    summary = run_build_all(
        manifest=Path("examples/demo_generalized/manifest.yaml"),
        output_root=tmp_path / "demo",
        max_decoys=1,
    )
    assert summary["reaction_smarts_library_core"] == 1
    core = pd.read_csv(tmp_path / "demo" / "04_release" / "reaction_smarts_library.T1_core.tsv", sep="\t")
    assert core.loc[0, "rule_application_unit"] == "single_smarts_application"
    assert core.loc[0, "biochemical_step_granularity"] in {"likely_single_enzyme_step", "uncertain"}
    assert "family_annotation_available" in core.columns


def test_possible_composite_granularity_from_multiple_functional_changes():
    row = {
        "template_scope": "generalized_template",
        "reaction_smarts": "[C:1]>>[C:1][O]",
        "reaction_delta_json": '{"atom_delta":{"C":2,"O":2,"H":0},"heavy_atom_delta":4,"exact_mass_delta":58.0,"ring_count_delta":0}',
        "reaction_delta_fingerprint": "C+2|O+2|HA+4|M+58.0",
    }
    g = annotate_biochemical_step_granularity(row, "hydroxylation_or_oxygenation", "structural_delta_specific")
    assert g["rule_application_unit"] == "single_smarts_application"
    assert g["biochemical_step_granularity"] in {"possible_composite_step", "likely_single_enzyme_step"}


def test_rule_confidence_not_gated_by_family_evidence_for_t1_core(tmp_path: Path):
    source = tmp_path / "templates.tsv"
    pd.DataFrame([
        {
            "source_reaction_id": "R_NO_FAMILY",
            "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
            "reaction_type_source": "alcohol_oxidation",
            "direction": "forward",
            "ec_numbers": "1.1.-.-",
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
                "source_database": "Rhea",
                "evidence_layer": "T1_Bio_Core",
            }
        },
        "family_evidence": {"enabled": False, "path": ""},
    }), encoding="utf-8")
    summary = run_build_all(manifest=manifest, output_root=tmp_path / "out", max_decoys=0)
    assert summary["reaction_smarts_library_core"] == 1
    core = pd.read_csv(tmp_path / "out" / "04_release" / "reaction_smarts_library.T1_core.tsv", sep="\t")
    assert core.loc[0, "family_annotation_scope"] == "not_assigned_at_rule_stage"
