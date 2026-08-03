from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd

from .io import ensure_dir, read_table, write_json, write_table


REQUIRED_ANALYSIS_FILES = (
    "generation_expansion_summary.tsv",
    "physicochemical_descriptor_summary.tsv",
    "nearest_G0_similarity.tsv",
    "chemical_space_fingerprint_projection.tsv",
    "functional_state_counts.tsv",
    "functional_state_transition_definitions.tsv",
    "functional_state_transition_rule_summary.tsv",
    "functional_state_transition_summary.tsv",
    "reaction_grammar_usage.tsv",
    "reaction_grammar_use_concentration.tsv",
    "convergence_and_route_multiplicity.tsv",
    "latent_bridge_candidates.tsv",
    "known_G0_pair_bridge_summary.tsv",
    "reaction_edit_landscape.tsv",
    "generation_rejection_summary.tsv",
    "chemical_space_analysis_summary.json",
)

REQUIRED_SPACE_FILES = (
    "chemical_space_nodes.tsv",
    "derivation_events.tsv",
    "rule_application_audit.tsv",
    "rejection_events.tsv",
    "generation_parent_progress.tsv",
    "taxane_reaction_grammar_space.sqlite",
    "chemical_space_build_summary.json",
    "G1_generation_summary.json",
    "G2_generation_summary.json",
    "G3_generation_summary.json",
)

REQUIRED_MANUSCRIPT_FILES = (
    "MANUSCRIPT_FINAL.md",
    "MATERIALS_AND_METHODS.md",
    "FIGURE_LEGENDS.md",
    "SUPPLEMENTARY_INFORMATION.md",
    "REFERENCES.md",
)

FORBIDDEN_MANUSCRIPT_PATTERNS = {
    "EC annotation": re.compile(r"\bEC(?:\s+number|\s+annotation)?\b"),
    "genome": re.compile(r"\bgenom(?:e|ic|ics)\b", re.IGNORECASE),
    "protein": re.compile(r"\bprotein(?:s)?\b", re.IGNORECASE),
    "docking": re.compile(r"\bdocking\b", re.IGNORECASE),
    "TDCN": re.compile(r"\bTDCN\b"),
}


def _contains_cjk(text: str) -> bool:
    return bool(re.search(r"[\u3400-\u4DBF\u4E00-\u9FFF]", text))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _check(
    rows: list[dict[str, Any]],
    name: str,
    passed: bool,
    evidence: str,
    *,
    severity: str = "error",
) -> None:
    rows.append(
        {
            "check": name,
            "severity": severity,
            "status": "pass" if passed else "fail",
            "evidence": evidence,
        }
    )


def _existing_required(
    directory: Path, filenames: tuple[str, ...]
) -> tuple[bool, list[str]]:
    missing = [
        filename
        for filename in filenames
        if not (directory / filename).is_file()
    ]
    return not missing, missing


def _read_text_files(paths: list[Path]) -> str:
    return "\n".join(
        path.read_text(encoding="utf-8") for path in paths if path.is_file()
    )


def _manifest_files(
    directories: list[tuple[str, Path]],
    output_dir: Path,
) -> pd.DataFrame:
    rows = []
    output_resolved = output_dir.resolve()
    for category, directory in directories:
        if not directory.exists():
            continue
        for path in sorted(directory.rglob("*")):
            if not path.is_file():
                continue
            if "__pycache__" in path.parts or ".pytest_cache" in path.parts:
                continue
            if path.suffix in {".pyc", ".pyo"}:
                continue
            try:
                path.resolve().relative_to(output_resolved)
                continue
            except ValueError:
                pass
            rows.append(
                {
                    "category": category,
                    "path": str(Path(category) / path.relative_to(directory)),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
            )
    return pd.DataFrame(rows)


def audit_release(
    code_dir: Path,
    provenance_dir: Path,
    space_dir: Path,
    validation_dir: Path,
    analysis_dir: Path,
    figures_dir: Path,
    manuscript_dir: Path,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    checks: list[dict[str, Any]] = []

    required_code_files = (
        "pyproject.toml",
        "README.md",
        "docs/STUDY_PROTOCOL.md",
        "docs/REPRODUCIBLE_WORKFLOW.md",
        "src/taxane_reaction_grammar_space/generate.py",
        "src/taxane_reaction_grammar_space/analyze.py",
        "src/taxane_reaction_grammar_space/study_figures.py",
        "src/taxane_reaction_grammar_space/manuscript.py",
        "src/taxane_reaction_grammar_space/audit.py",
    )
    code_ok, missing_code = _existing_required(
        code_dir, required_code_files
    )
    _check(
        checks,
        "complete_reproducible_code_package",
        code_ok,
        "missing=" + ",".join(missing_code) if missing_code else "all present",
    )

    required_provenance_names = {
        "software_environment.tsv",
        "input_file_hashes.tsv",
        "environment_and_input_audit.json",
        "rule_library_provenance_summary.json",
        "test_results.txt",
    }
    observed_provenance_names = {
        path.name for path in provenance_dir.rglob("*") if path.is_file()
    }
    missing_provenance = sorted(
        required_provenance_names - observed_provenance_names
    )
    _check(
        checks,
        "complete_environment_input_and_rule_provenance",
        not missing_provenance,
        (
            "all present"
            if not missing_provenance
            else "missing=" + ",".join(missing_provenance)
        ),
    )
    environment_summary_candidates = list(
        provenance_dir.rglob("environment_and_input_audit.json")
    )
    if environment_summary_candidates:
        environment_summary = json.loads(
            environment_summary_candidates[0].read_text(encoding="utf-8")
        )
        _check(
            checks,
            "input_snapshot_hashes_match",
            environment_summary.get("status") == "pass"
            and not environment_summary.get("input_hash_mismatches"),
            json.dumps(environment_summary, sort_keys=True),
        )
    else:
        _check(
            checks,
            "input_snapshot_hashes_match",
            False,
            "missing environment_and_input_audit.json",
        )
    test_result_candidates = list(provenance_dir.rglob("test_results.txt"))
    test_result_text = (
        test_result_candidates[0].read_text(
            encoding="utf-8", errors="replace"
        )
        if test_result_candidates
        else ""
    )
    _check(
        checks,
        "automated_test_suite_passes",
        bool(re.search(r"\nOK\s*$", test_result_text)),
        (
            "unittest report ends in OK"
            if re.search(r"\nOK\s*$", test_result_text)
            else "missing or nonpassing test report"
        ),
    )

    space_ok, missing_space = _existing_required(
        space_dir, REQUIRED_SPACE_FILES
    )
    _check(
        checks,
        "complete_G0_G3_space_artifacts",
        space_ok,
        "missing=" + ",".join(missing_space) if missing_space else "all present",
    )

    analysis_ok, missing_analysis = _existing_required(
        analysis_dir, REQUIRED_ANALYSIS_FILES
    )
    _check(
        checks,
        "complete_analysis_artifacts",
        analysis_ok,
        (
            "missing=" + ",".join(missing_analysis)
            if missing_analysis
            else "all present"
        ),
    )

    manuscript_ok, missing_manuscript = _existing_required(
        manuscript_dir, REQUIRED_MANUSCRIPT_FILES
    )
    _check(
        checks,
        "complete_manuscript_artifacts",
        manuscript_ok,
        (
            "missing=" + ",".join(missing_manuscript)
            if missing_manuscript
            else "all present"
        ),
    )

    generation_path = analysis_dir / "generation_expansion_summary.tsv"
    if generation_path.is_file():
        generation = read_table(generation_path)
        generation_ids = set(
            pd.to_numeric(generation["generation"], errors="coerce")
            .dropna()
            .astype(int)
        )
        _check(
            checks,
            "analysis_contains_exact_G0_G3_layers",
            generation_ids == {0, 1, 2, 3},
            f"observed_generations={sorted(generation_ids)}",
        )
        layer_map = {
            int(row["generation"]): str(row["interpretation_layer"])
            for row in generation.to_dict("records")
        }
        _check(
            checks,
            "generation_interpretation_boundary",
            layer_map.get(0) == "known_taxane_seed_space"
            and layer_map.get(1) == "primary_near_seed_chemical_space"
            and layer_map.get(2) == "primary_near_seed_chemical_space"
            and layer_map.get(3) == "exploratory_frontier",
            json.dumps(layer_map, sort_keys=True),
        )
    else:
        _check(
            checks,
            "analysis_contains_exact_G0_G3_layers",
            False,
            f"missing={generation_path}",
        )
        _check(
            checks,
            "generation_interpretation_boundary",
            False,
            f"missing={generation_path}",
        )

    analysis_summary_path = (
        analysis_dir / "chemical_space_analysis_summary.json"
    )
    convergence_path = (
        analysis_dir / "convergence_and_route_multiplicity.tsv"
    )
    if analysis_summary_path.is_file() and convergence_path.is_file():
        analysis_summary = json.loads(
            analysis_summary_path.read_text(encoding="utf-8")
        )
        convergence_columns = set(
            pd.read_csv(convergence_path, sep="\t", nrows=0).columns
        )
        required_route_columns = {
            "unique_semantic_group_count",
            "all_unique_parent_count",
            "later_rediscovery_event_count",
            "structural_path_count",
            "semantic_edge_path_count",
            "raw_rule_event_path_count",
        }
        route_definition_ok = (
            analysis_summary.get("convergence_definition")
            == (
                "at_least_two_distinct_parent_structures_in_the_targets_"
                "first_observed_generation"
            )
            and analysis_summary.get("path_count_layers", {}).get("primary")
            == "distinct_source_target_structural_edges"
            and analysis_summary.get("path_event_inclusion")
            == (
                "all_events_recorded_in_the_targets_first_observed_"
                "generation_regardless_of_insertion_flag"
            )
            and required_route_columns.issubset(convergence_columns)
        )
        _check(
            checks,
            "route_multiplicity_is_template_redundancy_aware",
            route_definition_ok,
            json.dumps(
                {
                    "convergence_definition": analysis_summary.get(
                        "convergence_definition"
                    ),
                    "path_count_layers": analysis_summary.get(
                        "path_count_layers"
                    ),
                    "path_event_inclusion": analysis_summary.get(
                        "path_event_inclusion"
                    ),
                    "missing_columns": sorted(
                        required_route_columns - convergence_columns
                    ),
                },
                sort_keys=True,
            ),
        )
    else:
        _check(
            checks,
            "route_multiplicity_is_template_redundancy_aware",
            False,
            (
                f"missing={analysis_summary_path}"
                if not analysis_summary_path.is_file()
                else f"missing={convergence_path}"
            ),
        )

    validation_summary_path = (
        validation_dir / "generated_space_validation_summary.json"
    )
    if validation_summary_path.is_file():
        validation_summary = json.loads(
            validation_summary_path.read_text(encoding="utf-8")
        )
        _check(
            checks,
            "generated_space_validation_passes",
            validation_summary.get("validation_status") == "pass"
            and int(validation_summary.get("failed_error_checks", -1)) == 0,
            json.dumps(validation_summary, sort_keys=True),
        )
    else:
        _check(
            checks,
            "generated_space_validation_passes",
            False,
            f"missing={validation_summary_path}",
        )

    scaffold_columns: list[str] = []
    if analysis_dir.exists():
        for path in analysis_dir.glob("*.tsv"):
            try:
                columns = list(pd.read_csv(path, sep="\t", nrows=0).columns)
            except pd.errors.EmptyDataError:
                continue
            scaffold_columns.extend(
                f"{path.name}:{column}"
                for column in columns
                if "scaffold" in column.lower()
            )
    _check(
        checks,
        "analysis_is_scaffold_independent",
        not scaffold_columns,
        "none" if not scaffold_columns else ";".join(scaffold_columns),
    )

    source_manifest_path = figures_dir / "figure_source_data_manifest.tsv"
    if source_manifest_path.is_file():
        source_manifest = read_table(source_manifest_path)
        missing_sources = [
            value
            for value in source_manifest["source_data_file"].astype(str)
            if not (
                Path(value)
                if Path(value).is_absolute()
                else figures_dir / Path(value)
            ).is_file()
        ]
        figure_names = set(source_manifest["figure"].astype(str))
        expected_figures = {
            *(f"Fig{value}" for value in range(1, 6)),
            *(f"FigS{value}" for value in range(1, 9)),
        }
        _check(
            checks,
            "all_main_and_supplementary_panels_have_source_data",
            len(source_manifest) >= 43
            and expected_figures.issubset(figure_names)
            and not missing_sources,
            (
                f"panel_sources={len(source_manifest)};"
                f"figures={sorted(figure_names)};"
                f"missing_sources={missing_sources}"
            ),
        )
        figure_source_scaffold_columns = []
        for value in source_manifest["source_data_file"].astype(str):
            path = Path(value)
            if not path.is_absolute():
                path = figures_dir / path
            if not path.is_file():
                continue
            try:
                columns = list(pd.read_csv(path, sep="\t", nrows=0).columns)
            except pd.errors.EmptyDataError:
                continue
            figure_source_scaffold_columns.extend(
                f"{path.name}:{column}"
                for column in columns
                if "scaffold" in column.lower()
            )
        _check(
            checks,
            "figure_source_data_are_scaffold_independent",
            not figure_source_scaffold_columns,
            (
                "none"
                if not figure_source_scaffold_columns
                else ";".join(figure_source_scaffold_columns)
            ),
        )
    else:
        _check(
            checks,
            "all_main_and_supplementary_panels_have_source_data",
            False,
            f"missing={source_manifest_path}",
        )
        _check(
            checks,
            "figure_source_data_are_scaffold_independent",
            False,
            f"missing={source_manifest_path}",
        )

    expected_image_stems = [
        "Fig1_evidence_stratified_reaction_grammar",
        "Fig2_validation_and_evidence_sensitivity",
        "Fig3_iterative_taxane_space_expansion",
        "Fig4_reaction_grammar_transformation_landscape",
        "Fig5_convergence_and_latent_bridges",
        "FigS1_fingerprint_projection",
        "FigS2_physicochemical_distributions",
        "FigS3_rule_activation_across_generations",
        "FigS4_product_rejection_reasons",
        "FigS5_cycles_and_known_space_reconnections",
        "FigS6_extended_functional_state_transitions",
        "FigS7_reaction_edit_similarity_and_locality",
        "FigS8_evidence_layer_G1_sensitivity",
    ]
    missing_images = [
        f"{stem}.{suffix}"
        for stem in expected_image_stems
        for suffix in ("pdf", "svg", "png")
        if not (figures_dir / f"{stem}.{suffix}").is_file()
    ]
    _check(
        checks,
        "all_figures_have_pdf_svg_and_png",
        not missing_images,
        "all present" if not missing_images else ",".join(missing_images),
    )

    manuscript_paths = [
        manuscript_dir / filename for filename in REQUIRED_MANUSCRIPT_FILES
    ]
    manuscript_text = _read_text_files(manuscript_paths)
    narrative_manuscript_text = _read_text_files(
        [
            manuscript_dir / filename
            for filename in REQUIRED_MANUSCRIPT_FILES
            if filename != "REFERENCES.md"
        ]
    )
    _check(
        checks,
        "manuscript_has_no_unresolved_final_placeholders",
        "{{FINAL_" not in manuscript_text,
        "none" if "{{FINAL_" not in manuscript_text else "found {{FINAL_",
    )
    has_cjk = _contains_cjk(manuscript_text)
    _check(
        checks,
        "manuscript_is_English_only",
        not has_cjk,
        "no CJK characters" if not has_cjk else "CJK found",
    )
    forbidden_hits = [
        name
        for name, pattern in FORBIDDEN_MANUSCRIPT_PATTERNS.items()
        if pattern.search(narrative_manuscript_text)
    ]
    _check(
        checks,
        "manuscript_excludes_downstream_annotation_topics",
        not forbidden_hits,
        "none" if not forbidden_hits else ",".join(forbidden_hits),
    )

    check_frame = pd.DataFrame(checks)
    failed = check_frame[
        (check_frame["severity"] == "error")
        & (check_frame["status"] == "fail")
    ]
    checks_path = output_dir / "release_audit_checks.tsv"
    write_table(check_frame, checks_path)

    manifest = _manifest_files(
        [
            ("code", code_dir),
            ("00_provenance", provenance_dir),
            ("03_primary_G0_G3", space_dir),
            ("04_validation", validation_dir),
            ("06_analysis", analysis_dir),
            ("07_figures", figures_dir),
            ("08_manuscript", manuscript_dir),
        ],
        output_dir,
    )
    manifest_path = output_dir / "release_file_manifest.tsv"
    write_table(manifest, manifest_path)
    checksum_path = output_dir / "release_checksums.sha256"
    checksum_path.write_text(
        "".join(
            f"{row.sha256}  {row.path}\n"
            for row in manifest.itertuples(index=False)
        ),
        encoding="utf-8",
    )
    summary = {
        "release_status": "pass" if failed.empty else "fail",
        "checks": int(len(check_frame)),
        "passed_checks": int((check_frame["status"] == "pass").sum()),
        "failed_error_checks": int(len(failed)),
        "failed_error_names": list(failed["check"]),
        "manifested_files": int(len(manifest)),
        "manifested_bytes": int(manifest["size_bytes"].sum())
        if not manifest.empty
        else 0,
        "G3_interpretation": "exploratory_frontier",
        "outputs": {
            "checks": checks_path.name,
            "manifest": manifest_path.name,
            "checksums": checksum_path.name,
        },
    }
    summary_path = output_dir / "release_audit_summary.json"
    write_json(summary, summary_path)
    if not failed.empty:
        raise RuntimeError(
            "Release audit failed: " + ", ".join(failed["check"])
        )
    return {
        "checks": checks_path,
        "manifest": manifest_path,
        "checksums": checksum_path,
        "summary": summary_path,
    }
