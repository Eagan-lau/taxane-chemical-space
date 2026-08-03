from __future__ import annotations

from pathlib import Path
from typing import Any
import re

from .utils import ensure_dir, write_json, write_yaml


def _path_text(path: Path) -> str:
    return str(path.resolve())


def _file_info(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"path": _path_text(path), "exists": False}
    return {
        "path": _path_text(path),
        "exists": True,
        "bytes": int(path.stat().st_size),
    }


def _with_limit(dataset: dict[str, Any], max_records: int | None) -> dict[str, Any]:
    if max_records is not None and int(max_records) > 0:
        dataset["max_records"] = int(max_records)
    return dataset


def _looks_like_external_db_root(path: Path) -> bool:
    """Return whether a directory has the raw database layout expected by v0.4.1."""
    if not path.exists() or not path.is_dir():
        return False
    markers = [
        path / "index" / "rhea_reaction_evidence.csv",
        path / "index" / "kegg_reaction_evidence.csv",
        path / "index" / "metanetx_reaction_evidence.csv",
        path / "retrorules" / "retrorules_2019_rr01_sqlite" / "mvc.db",
        path / "bionavi_np" / "processed" / "biochem_train.rdkit_valid_reactions.txt",
    ]
    return any(p.exists() for p in markers)


def resolve_external_db_root(external_db_root: str | Path) -> tuple[Path, dict[str, Any]]:
    """Accept either a dated raw database directory or its parent directory.

    For example, both of these are valid:

    - `/path/to/external_databases/2026-06-09`
    - `/path/to/external_databases`

    When a parent directory is supplied, the latest dated child that contains
    recognized raw database files is selected.
    """
    requested = Path(external_db_root).resolve()
    info: dict[str, Any] = {
        "requested_external_db_root": _path_text(requested),
        "resolved_external_db_root": _path_text(requested),
        "root_resolution_mode": "as_provided",
        "candidate_roots": [],
    }
    if _looks_like_external_db_root(requested):
        return requested, info
    if not requested.exists() or not requested.is_dir():
        info["root_resolution_mode"] = "path_missing"
        return requested, info

    candidates = []
    for child in requested.iterdir():
        if child.is_dir() and _looks_like_external_db_root(child):
            dated = bool(re.fullmatch(r"\d{4}-\d{2}-\d{2}", child.name))
            candidates.append((dated, child.name, child.stat().st_mtime, child.resolve()))
    candidates.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
    info["candidate_roots"] = [_path_text(c[3]) for c in candidates]
    if candidates:
        resolved = candidates[0][3]
        info["resolved_external_db_root"] = _path_text(resolved)
        info["root_resolution_mode"] = "latest_dated_child" if candidates[0][0] else "latest_matching_child"
        return resolved, info
    info["root_resolution_mode"] = "no_matching_child_found"
    return requested, info


def _add_dataset(
    datasets: dict[str, dict[str, Any]],
    discovery: dict[str, Any],
    name: str,
    *,
    path: Path,
    parser: str,
    source_database: str,
    evidence_layer: str | None = None,
    max_records: int | None = None,
    extra: dict[str, Any] | None = None,
) -> None:
    info = _file_info(path)
    if not path.exists():
        discovery["skipped"].append({"dataset": name, "reason": "missing_file", **info})
        return
    row: dict[str, Any] = {
        "enabled": True,
        "parser": parser,
        "path": _path_text(path),
        "source_database": source_database,
    }
    if evidence_layer:
        row["evidence_layer"] = evidence_layer
    if extra:
        row.update(extra)
    datasets[name] = _with_limit(row, max_records)
    discovery["included"].append({
        "dataset": name,
        "parser": parser,
        "source_database": source_database,
        "evidence_layer": evidence_layer or "inferred",
        "max_records": int(max_records) if max_records is not None and int(max_records) > 0 else None,
        **info,
    })


def default_taxol_pathway_path() -> Path:
    return Path(__file__).resolve().parents[2] / "examples" / "taxol_pathway.csv"


def resolve_taxol_pathway_path(taxol_pathway: str | Path | None = None) -> Path:
    return Path(taxol_pathway).resolve() if taxol_pathway else default_taxol_pathway_path().resolve()


def build_raw_source_manifest(
    *,
    external_db_root: str | Path,
    output_dir: str | Path,
    taxol_pathway: str | Path | None = None,
    include_taxol_anchors: bool | None = None,
    include_uspto: bool = True,
    include_annotation_only: bool = True,
    use_processed_bionavi: bool = True,
    max_retrorules: int | None = None,
    max_rhea: int | None = None,
    max_bionavi_per_file: int | None = None,
    max_uspto: int | None = None,
    max_kegg: int | None = None,
    max_metanetx: int | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Create a v0.4.1-native manifest from the local raw database directory.

    The generated manifest deliberately avoids old `network_ready_*` intermediate
    files. Exact reaction sources can later be generalized by
    `--abstract-exact-reactions`; existing SMARTS sources such as RetroRules are
    imported directly.
    """
    requested_root = Path(external_db_root).resolve()
    root, root_resolution = resolve_external_db_root(external_db_root)
    taxol_path = resolve_taxol_pathway_path(taxol_pathway)
    if include_taxol_anchors is None:
        include_taxol_anchors = bool(taxol_pathway) and taxol_path.exists()
    out_dir = ensure_dir(output_dir)
    datasets: dict[str, dict[str, Any]] = {}
    discovery: dict[str, Any] = {
        **root_resolution,
        "external_db_root": _path_text(root),
        "manifest_policy": "v0.4.1-native raw database discovery; no network_ready inputs",
        "taxol_pathway": _file_info(taxol_path),
        "taxol_pathway_default": _path_text(default_taxol_pathway_path().resolve()),
        "taxol_pathway_include_requested": bool(include_taxol_anchors),
        "included": [],
        "skipped": [],
    }

    if include_taxol_anchors:
        _add_dataset(
            datasets,
            discovery,
            "taxol_known_pathway_anchors",
            path=taxol_path,
            parser="taxol_pathway_csv",
            source_database="TaxolKnownPathway_Curated",
            evidence_layer="T1_Bio_Core",
            extra={
                "curated_taxol_anchor": "true",
                "curated_pathway_name": "TaxolKnownPathway",
            },
        )
    else:
        discovery["skipped"].append({
            "dataset": "taxol_known_pathway_anchors",
            "reason": "disabled_by_user_or_missing_taxol_path",
            **_file_info(taxol_path),
        })

    _add_dataset(
        datasets,
        discovery,
        "retrorules_sqlite",
        path=root / "retrorules" / "retrorules_2019_rr01_sqlite" / "mvc.db",
        parser="retrorules_sqlite",
        source_database="RetroRules",
        max_records=max_retrorules,
    )

    _add_dataset(
        datasets,
        discovery,
        "rhea_reaction_evidence",
        path=root / "index" / "rhea_reaction_evidence.csv",
        parser="reaction_evidence_csv",
        source_database="Rhea",
        max_records=max_rhea,
    )

    bionavi_processed = root / "bionavi_np" / "processed"
    if use_processed_bionavi:
        for stem in ["biochem_train", "biochem_valid", "biochem_test"]:
            _add_dataset(
                datasets,
                discovery,
                f"bionavi_np_{stem}",
                path=bionavi_processed / f"{stem}.rdkit_valid_reactions.txt",
                parser="reaction_smiles_lines",
                source_database="BioNaviNP_BioChem",
                evidence_layer="T2_Bio_Extended",
                max_records=max_bionavi_per_file,
                extra={"source_dataset": stem},
            )
    else:
        _add_dataset(
            datasets,
            discovery,
            "bionavi_np_biochem_zip",
            path=root / "bionavi_np" / "biochem.zip",
            parser="reaction_smiles_zip",
            source_database="BioNaviNP_BioChem",
            evidence_layer="T2_Bio_Extended",
            max_records=max_bionavi_per_file,
            extra={
                "members": ["biochem/train.txt", "biochem/valid.txt", "biochem/test.txt"],
                "source_dataset": "biochem",
            },
        )

    if include_uspto:
        _add_dataset(
            datasets,
            discovery,
            "bionavi_np_uspto_npl_np_like",
            path=bionavi_processed / "uspto_npl_np_like.rdkit_valid_reactions.txt",
            parser="reaction_smiles_lines",
            source_database="BioNaviNP_USPTO_NPL",
            evidence_layer="T3_Chem_like",
            max_records=max_uspto,
            extra={"source_dataset": "uspto_npl_np_like"},
        )
    else:
        discovery["skipped"].append({
            "dataset": "bionavi_np_uspto_npl_np_like",
            "reason": "disabled",
            **_file_info(bionavi_processed / "uspto_npl_np_like.rdkit_valid_reactions.txt"),
        })

    if include_annotation_only:
        _add_dataset(
            datasets,
            discovery,
            "kegg_reaction_evidence",
            path=root / "index" / "kegg_reaction_evidence.csv",
            parser="reaction_evidence_csv",
            source_database="KEGG",
            max_records=max_kegg,
        )
        _add_dataset(
            datasets,
            discovery,
            "metanetx_reaction_evidence",
            path=root / "index" / "metanetx_reaction_evidence.csv",
            parser="reaction_evidence_csv",
            source_database="MetaNetX",
            max_records=max_metanetx,
        )
    else:
        for name, path in {
            "kegg_reaction_evidence": root / "index" / "kegg_reaction_evidence.csv",
            "metanetx_reaction_evidence": root / "index" / "metanetx_reaction_evidence.csv",
        }.items():
            discovery["skipped"].append({"dataset": name, "reason": "annotation_only_disabled", **_file_info(path)})

    manifest = {
        "build_metadata": {
            "builder": "enzymatic_rule_builder_v0.4.1",
            "input_policy": "raw_external_databases_only",
            "legacy_network_ready_inputs": "not_used",
            "requested_external_db_root": _path_text(requested_root),
            "external_db_root": _path_text(root),
            "root_resolution_mode": root_resolution.get("root_resolution_mode", "as_provided"),
            "taxol_pathway": _path_text(taxol_path),
            "taxol_pathway_included": bool(include_taxol_anchors and taxol_path.exists()),
        },
        "datasets": datasets,
        "family_evidence": {
            "enabled": False,
            "path": "",
            "delimiter": "\t",
        },
    }
    discovery["dataset_count"] = len(datasets)
    discovery["included_dataset_names"] = list(datasets.keys())
    discovery["warnings"] = []
    if not datasets:
        discovery["warnings"].append("No raw datasets were discovered. Check external_db_root.")
    if include_annotation_only:
        discovery["warnings"].append(
            "KEGG and MetaNetX rows are annotation/cross-reference evidence unless paired with computable reaction SMILES or SMARTS."
        )
    if include_uspto:
        discovery["warnings"].append("USPTO_NPL is included as T3_Chem_like exploratory evidence, not strict biochemical core evidence.")

    manifest_path = out_dir / "source_manifest.raw_v0.4.1.yaml"
    discovery_path = out_dir / "raw_database_discovery.json"
    discovery["manifest_path"] = _path_text(manifest_path)
    discovery["discovery_path"] = _path_text(discovery_path)
    write_yaml(manifest_path, manifest)
    write_json(discovery_path, discovery)
    return manifest_path, discovery
