from pathlib import Path

import pandas as pd

from enzymatic_rule_builder.raw import build_raw_source_manifest, resolve_external_db_root
from enzymatic_rule_builder.sources import load_datasets, load_reaction_evidence_csv, load_reaction_smiles_lines


def test_reaction_smiles_lines_loader(tmp_path: Path):
    path = tmp_path / "rxns.txt"
    path.write_text("CCO>>CC=O\nnot_a_reaction\nCC(=O)O>>CC(=O)OC\n", encoding="utf-8")
    df = load_reaction_smiles_lines(
        path,
        {
            "source_database": "BioNaviNP_BioChem",
            "evidence_layer": "T2_Bio_Extended",
            "source_dataset": "mini",
        },
    )
    assert len(df) == 2
    assert set(df["reaction_smiles"]) == {"CCO>>CC=O", "CC(=O)O>>CC(=O)OC"}
    assert df["evidence_layer"].eq("T2_Bio_Extended").all()


def test_reaction_evidence_csv_loader_maps_ids(tmp_path: Path):
    path = tmp_path / "rhea_reaction_evidence.csv"
    pd.DataFrame(
        [
            {
                "source_database": "Rhea",
                "database_reaction_id": "10001",
                "source_reaction_id": "RHEA:10001",
                "reaction_smiles": "CCO>>CC=O",
                "reaction_equation": "ethanol = acetaldehyde",
                "ec_numbers": "1.1.1.1",
                "kegg_reaction_ids": "R00623",
                "cross_references": "RHEA:10001; KEGG:R00623",
                "is_balanced": "true",
                "is_transport": "false",
            }
        ]
    ).to_csv(path, index=False)
    df = load_reaction_evidence_csv(path, {"source_database": "Rhea"})
    assert len(df) == 1
    assert df.loc[0, "rhea_ids"] == "RHEA:10001"
    assert df.loc[0, "kegg_ids"] == "R00623"
    assert df.loc[0, "ec_numbers"] == "1.1.1.1"
    assert df.loc[0, "evidence_layer"] == "T1_Bio_Core"


def test_raw_manifest_uses_raw_files_not_network_ready(tmp_path: Path):
    root = tmp_path / "external"
    (root / "index").mkdir(parents=True)
    (root / "bionavi_np" / "processed").mkdir(parents=True)
    (root / "index" / "rhea_reaction_evidence.csv").write_text(
        "source_database,database_reaction_id,source_reaction_id,reaction_smiles,ec_numbers\n"
        "Rhea,10001,RHEA:10001,CCO>>CC=O,1.1.1.1\n",
        encoding="utf-8",
    )
    (root / "bionavi_np" / "processed" / "biochem_train.rdkit_valid_reactions.txt").write_text("CCO>>CC=O\n", encoding="utf-8")
    manifest, discovery = build_raw_source_manifest(
        external_db_root=root,
        output_dir=tmp_path / "manifest",
        include_uspto=False,
        include_annotation_only=False,
        max_rhea=1,
        max_bionavi_per_file=1,
    )
    assert manifest.exists()
    assert discovery["dataset_count"] == 2
    assert discovery["manifest_policy"].startswith("v0.4.1-native")
    source_df, used = load_datasets(manifest)
    assert len(used) == 2
    assert len(source_df) == 2
    assert not any("network_ready" in p for p in used)


def test_external_db_parent_resolves_latest_dated_child(tmp_path: Path):
    parent = tmp_path / "external_databases"
    old = parent / "2026-01-01"
    new = parent / "2026-06-09"
    for root in [old, new]:
        (root / "index").mkdir(parents=True)
        (root / "index" / "rhea_reaction_evidence.csv").write_text(
            "source_database,database_reaction_id,source_reaction_id,reaction_smiles,ec_numbers\n"
            "Rhea,10001,RHEA:10001,CCO>>CC=O,1.1.1.1\n",
            encoding="utf-8",
        )

    resolved, info = resolve_external_db_root(parent)
    assert resolved == new.resolve()
    assert info["root_resolution_mode"] == "latest_dated_child"

    manifest, discovery = build_raw_source_manifest(
        external_db_root=parent,
        output_dir=tmp_path / "manifest_parent",
        include_uspto=False,
        include_annotation_only=False,
        max_rhea=1,
    )
    assert discovery["requested_external_db_root"] == str(parent.resolve())
    assert discovery["resolved_external_db_root"] == str(new.resolve())
    assert discovery["taxol_pathway"]["exists"] is True
    source_df, used = load_datasets(manifest)
    assert len(source_df) == 1
    assert str(new) in used[0]
