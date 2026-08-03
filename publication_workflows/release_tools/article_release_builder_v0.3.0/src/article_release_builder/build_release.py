from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

from .data_exports import (
    TABLE_GROUP_TITLES,
    build_supplementary_workbook,
    export_bridge_tables,
    export_generation_csvs,
)
from .documents import convert_docx_to_pdf, markdown_to_docx
from .figures import (
    SUPPLEMENTARY_FIGURE_GROUPS,
    build_curated_supplementary_figures,
)
from .utils import (
    copy_file,
    copy_tree,
    ensure_empty_output,
    hardlink_tree,
    package_manifest,
    sha256_file,
    write_json,
    write_tsv,
)
from .validation import validate_completed_release


ARTICLE_TITLE = (
    "A provenance-resolved enzymatic reaction grammar delineates "
    "the generative chemical space of taxanes"
)

MAIN_FIGURE_RELEASES = {
    1: "figure1_redesign_v6_20260730",
    2: "figure2_integrated_v12_balanced_g2_20260731",
    3: "figure3_g0_g3_v7_20260730",
    4: "figure4_redesigned_v8_20260730",
}

MAIN_FIGURE_STEMS = {
    1: "Figure_1_Evidence_Stratified_Grammar_V6",
    2: "Figure_2_T1_Taxane_Chemical_Space_V12_Balanced_G2",
    3: "Figure_3_T1_Chemical_Regularities_V7",
    4: "Figure_4_Bridge_Hypotheses_to_Global_Topology_V8",
}

MAIN_CAPTIONS = {
    1: "FIGURE_1_CAPTION.md",
    2: "FIGURE_2_CAPTION_V12.md",
    3: "FIGURE_3_CAPTION_V7.md",
    4: "FIGURE_4_CAPTION_V8.md",
}

MAIN_AUDITS = {
    1: "FIGURE_1_NUMERICAL_AUDIT.tsv",
    2: "NUMERICAL_AUDIT_V12.tsv",
    3: "NUMERICAL_AUDIT_V7.tsv",
    4: "NUMERICAL_AUDIT_FIGURE_4_V8.tsv",
}

SUPPLEMENTARY_CAPTIONS = {
    "FIGURE_S1": (
        "**Supplementary Figure S1 | Reaction-source provenance and "
        "rule-build attrition.** **(A)** Contributions of the frozen source "
        "collections to normalized reaction records. Source counts are "
        "overlapping support counts where appropriate and do not imply "
        "equivalent evidence quality. **(B)** Tier-specific quality-control "
        "attrition and semantic compression from release rows to executable "
        "representatives."
    ),
    "FIGURE_S2": (
        "**Supplementary Figure S2 | Evidence-tier sensitivity of the "
        "one-step chemical space.** **(A)** Physicochemical shifts of the "
        "independently generated T1, T2, and T3 G1 sets relative to G0. "
        "**(B)** Distributions of maximum Morgan-Tanimoto similarity to G0. "
        "The tiers were evaluated independently and were not pooled."
    ),
    "FIGURE_S3": (
        "**Supplementary Figure S3 | Product quality control and executable-"
        "grammar activation.** **(A)** Complete generation-resolved rejection "
        "audit by product-quality category. **(B)** Rule activation and use "
        "of the frozen T1 grammar across G0-G3. Rejected applications and "
        "inactive rules remain available in the audit tables."
    ),
    "FIGURE_S4": (
        "**Supplementary Figure S4 | Physicochemical displacement and "
        "reaction-edit locality.** **(A)** Generation-specific distributions "
        "of molecular descriptors. **(B)** Source-product similarity and "
        "structure-derived edit locality among accepted directional events. "
        "Generation depth is grammar distance and not necessarily enzyme-step "
        "count."
    ),
    "FIGURE_S5": (
        "**Supplementary Figure S5 | Extended functional-state transition "
        "atlas.** Structure-derived changes in predefined, non-mutually "
        "exclusive functional states across accepted directional derivations. "
        "Transitions were recalculated from source and product structures and "
        "were not inferred from rule names."
    ),
    "FIGURE_S6": (
        "**Supplementary Figure S6 | Convergence, reverse cycles, and "
        "known-space reconnection.** Generation-resolved summaries distinguish "
        "first-observation parent convergence, immediate reverse cycles, and "
        "connectivity-level returns to G0. Multiplicity is topological support "
        "and not biological replication."
    ),
    "FIGURE_S7": (
        "**Supplementary Figure S7 | Molecular structures of high-support "
        "latent bridge candidates.** Displayed structures were selected "
        "deterministically from the frozen bridge analysis. They are "
        "reaction-grammar-accessible hypotheses and are not asserted as "
        "observed taxane metabolites."
    ),
    "FIGURE_S8": (
        "**Supplementary Figure S8 | Complete directed G0 reconnection "
        "network.** All 309 ordered known-taxane pairs supported by the 227 "
        "conservative G1/G2 bridge candidates are shown. The graph records "
        "directed grammar paths and does not establish organism-level "
        "biosynthetic routes."
    ),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-release", type=Path, required=True)
    parser.add_argument("--editorial-release", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--skip-raw-hardlinks",
        action="store_true",
        help="Do not create a self-contained hard-linked primary-release snapshot.",
    )
    return parser.parse_args(argv)


def require_paths(paths: list[Path]) -> None:
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required inputs:\n" + "\n".join(str(path) for path in missing)
        )


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8").strip() + "\n"


def normalize_embedded_heading(text: str, replacement: str) -> str:
    lines = text.splitlines()
    if lines and lines[0].startswith("# "):
        lines[0] = replacement
    return "\n".join(lines).strip() + "\n"


def normalize_caption(text: str) -> str:
    lines = text.strip().splitlines()
    if lines and lines[0].startswith("# "):
        heading = lines[0][2:].strip()
        lines[0] = f"**{heading}**"
    return "\n".join(lines).strip()


def copy_main_figures(
    editorial_release: Path,
    output: Path,
) -> tuple[dict[str, Path], list[Path], dict[str, str]]:
    image_map: dict[str, Path] = {}
    audit_paths: list[Path] = []
    caption_map: dict[str, str] = {}
    for number in range(1, 5):
        source_release = editorial_release / MAIN_FIGURE_RELEASES[number]
        stem = MAIN_FIGURE_STEMS[number]
        target_source = output / "source_data" / "main_figures" / f"Figure_{number}"
        copy_tree(source_release / "source_data", target_source / "panel_source_data")
        copy_tree(source_release / "workflow", target_source / "workflow")
        for metadata_name in (
            MAIN_CAPTIONS[number],
            MAIN_AUDITS[number],
            "README.md",
        ):
            candidate = source_release / metadata_name
            if candidate.exists():
                copy_file(candidate, target_source / metadata_name)
        for extension in ("pdf", "svg", "png"):
            source = source_release / "figures" / f"{stem}.{extension}"
            target = output / "figures" / "main" / f"Figure_{number}.{extension}"
            copy_file(source, target)
        key = f"FIGURE_{number}"
        image_map[key] = output / "figures" / "main" / f"Figure_{number}.png"
        audit_paths.append(source_release / MAIN_AUDITS[number])
        caption_map[key] = normalize_caption(
            read_text(source_release / MAIN_CAPTIONS[number])
        )
    return image_map, audit_paths, caption_map


def copy_supplementary_figures(
    editorial_release: Path,
    output: Path,
) -> dict[str, Path]:
    source_release = editorial_release / "manuscript_v4_figure_redesign"
    source_figures = source_release / "supplementary_figures"
    image_map = build_curated_supplementary_figures(
        source_figures,
        output / "figures" / "supplementary",
    )
    archive = output / "source_data" / "supplementary_figures"
    copy_tree(source_release / "source_data", archive / "panel_source_data")
    copy_tree(source_figures, archive / "original_15_figure_set")
    write_tsv(
        archive / "SUPPLEMENTARY_FIGURE_CONSOLIDATION.tsv",
        [
            ["new_figure", "source_figures", "editorial_action"],
            *[
                [
                    f"Figure S{new}",
                    ";".join(f"old Figure S{old}" for old in old_group),
                    "merged" if len(old_group) > 1 else "retained",
                ]
                for new, old_group in SUPPLEMENTARY_FIGURE_GROUPS.items()
            ],
            ["omitted", "old Figure S3", "duplicates main Figure 1E"],
            ["omitted", "old Figure S11", "decorative molecular gallery"],
            ["omitted", "old Figure S13", "overlaps final Figure S8"],
            ["omitted", "old Figure S14", "routes represented in main Figure 4"],
        ],
    )
    return image_map


def write_benchmark_table(primary_release: Path, target: Path) -> None:
    benchmark_root = primary_release / "05_sensitivity_and_benchmarks"
    sources = (
        benchmark_root / "external_leakage_control" / "benchmark_summary.json",
        benchmark_root / "domain_replay_calibration" / "benchmark_summary.json",
    )
    rows = []
    for source in sources:
        payload = json.loads(source.read_text(encoding="utf-8"))
        odds_ratio = payload.get("fisher_exact_odds_ratio")
        if isinstance(odds_ratio, float) and math.isnan(odds_ratio):
            odds_ratio = "NA"
        rows.append(
            [
                payload["mode"],
                payload["compiled_evaluated_rules"],
                payload["positive_reactions"],
                payload["positive_connectivity_recovered"],
                payload["positive_full_stereo_recovered"],
                payload["decoy_pairs"],
                payload["decoy_connectivity_matched"],
                odds_ratio,
                payload["fisher_exact_p_value"],
                "specificity control"
                if payload["exclude_taxol_derived"]
                else "internal domain-replay calibration",
            ]
        )
    write_tsv(
        target,
        [
            [
                "mode",
                "compiled_rules",
                "curated_reactions",
                "connectivity_recovered",
                "full_stereochemistry_recovered",
                "decoy_pairs",
                "decoy_connectivity_matches",
                "fisher_odds_ratio",
                "fisher_p_value",
                "interpretation",
            ],
            *rows,
        ],
    )


def prepare_supplementary_tables(
    editorial_release: Path,
    primary_release: Path,
    output: Path,
) -> dict[int, list[dict[str, object]]]:
    old_tables = (
        editorial_release / "manuscript_v4_figure_redesign" / "supplementary_tables"
    )
    archive = output / "source_data" / "supplementary_tables" / "original_S1-S20"
    copy_tree(old_tables, archive)
    target = output / "supplementary_tables" / "tsv"
    target.mkdir(parents=True, exist_ok=True)

    def component(
        group: int,
        suffix: str,
        old_number: int,
        title: str,
    ) -> dict[str, object]:
        destination = target / f"Table_S{group}{suffix}_{title.replace(' ', '_')}.tsv"
        copy_file(old_tables / f"Table_S{old_number}_V4.tsv", destination)
        return {
            "label": f"Table S{group}{suffix}",
            "title": title,
            "path": destination,
        }

    benchmark = target / "Table_S1E_Pathway_calibration.tsv"
    write_benchmark_table(primary_release, benchmark)
    bridge_candidates, directed_pairs = export_bridge_tables(
        primary_release / "06_analysis",
        target,
    )

    groups: dict[int, list[dict[str, object]]] = {
        1: [
            component(1, "A", 1, "Reaction_sources_and_provenance"),
            component(1, "B", 2, "Evidence_tiers_and_build_attrition"),
            component(1, "C", 3, "Independent_one_step_tier_outputs"),
            component(1, "D", 4, "G1_evidence_tier_set_overlap"),
            {
                "label": "Table S1E",
                "title": "Leakage-controlled and domain-informed pathway calibration",
                "path": benchmark,
            },
        ],
        2: [
            component(2, "A", 5, "Primary_T1_grammar_provenance"),
            component(2, "B", 6, "G0_rule_activation_audit"),
        ],
        3: [
            component(3, "", 7, "Generation_level_expansion_summary"),
        ],
        4: [
            component(4, "A", 8, "Molecular_state_table_manifest"),
            component(4, "B", 9, "Directed_derivation_event_manifest"),
            component(4, "C", 10, "Parent_rule_application_audit_manifest"),
            component(4, "D", 11, "Product_rejection_audit_manifest"),
        ],
        5: [
            component(5, "A", 12, "Quality_control_comparison"),
            component(5, "B", 13, "Parameter_robustness"),
            component(5, "C", 14, "Computational_performance"),
        ],
        6: [
            component(6, "", 15, "Descriptors_and_nearest_G0_similarity"),
        ],
        7: [
            component(7, "", 16, "Functional_state_transitions"),
        ],
        8: [
            component(8, "A", 17, "Grammar_usage_and_concentration"),
            component(8, "B", 18, "Reaction_edit_locality_and_formula_deltas"),
        ],
        9: [
            component(9, "", 19, "Convergence_and_route_multiplicity"),
        ],
        10: [
            {
                "label": "Table S10",
                "title": "All latent bridge candidates",
                "path": bridge_candidates,
            },
        ],
        11: [
            {
                "label": "Table S11",
                "title": "All directed G0-pair bridge records",
                "path": directed_pairs,
            },
        ],
    }
    write_tsv(
        output
        / "source_data"
        / "supplementary_tables"
        / "SUPPLEMENTARY_TABLE_CONSOLIDATION.tsv",
        [
            ["new_table", "title", "component", "authoritative_tsv"],
            *[
                [
                    f"Table S{number}",
                    TABLE_GROUP_TITLES[number],
                    component_row["label"],
                    Path(component_row["path"]).name,
                ]
                for number, components in groups.items()
                for component_row in components
            ],
            [
                "omitted",
                "Compact bridge summary",
                "old Table S20",
                "duplicated by complete Tables S10-S11",
            ],
        ],
    )
    return groups


def assemble_markdown(
    template: str,
    replacements: dict[str, str],
    image_relative_paths: dict[str, str],
    caption_map: dict[str, str] | None = None,
) -> str:
    caption_map = caption_map or {}
    text = template
    for key, value in replacements.items():
        text = text.replace(f"{{{{{key}}}}}", value.strip())
    for key, relative_path in image_relative_paths.items():
        replacement = f"![{key}]({relative_path})"
        caption = caption_map.get(key, "").strip()
        if caption:
            replacement += f"\n\n{caption}"
        text = text.replace(f"{{{{{key}}}}}", replacement)
    unresolved = [
        line for line in text.splitlines() if line.strip().startswith("{{")
    ]
    if unresolved:
        raise ValueError(f"Unresolved Markdown placeholders: {unresolved}")
    return text.strip() + "\n"


def write_release_readme(output: Path, raw_snapshot_status: dict[str, object]) -> None:
    text = f"""# Reproducible taxane reaction-grammar article release

This release assembles the frozen taxane reaction-grammar study into a
submission-oriented manuscript, figure set, curated Supplementary Information,
11 logically consolidated supplementary tables, generation-specific molecular
CSVs, source data, and reproduction code.

**{ARTICLE_TITLE}**

## Start here

- `manuscript/Main_Manuscript_with_Figures.docx`
- `manuscript/Main_Manuscript_with_Figures.pdf`
- `manuscript/Figure_Legends.docx`
- `supplementary_information/Supplementary_Information.docx`
- `supplementary_tables/Supplementary_Tables_S1-S11.xlsx`
- `data/generation_csv/G0_known_taxanes.csv`
- `data/generation_csv/G1_inferred_intermediates.csv`
- `data/generation_csv/G2_inferred_intermediates.csv`
- `data/generation_csv/G3_exploratory_intermediates.csv`

## Editorial consolidation

The four main figures preserve the frozen scientific results. Figure 2 uses
the V12 multiscale display: G0-G2 are node resolved, G3 is represented as a
complete density background, and the bounded G2 relaxation changes display
coordinates only. Main-figure captions appear immediately below their figures
in the manuscript. The 15-figure supplementary set was consolidated into eight
thematic figures. The 22-table set was consolidated into 11 logical tables;
component TSV files remain separate and checksum-bearing. Original
supplementary figures and tables are retained under `source_data` for complete
traceability.

## Interpretation boundary

G0 contains known taxane seeds. G1 and G2 define the primary near-seed
reaction-grammar space. G3 is exploratory and right-censored. Generated
structures are hypotheses of grammar accessibility rather than experimentally
observed metabolites.

## Raw-data snapshot

Raw snapshot status: `{json.dumps(raw_snapshot_status, sort_keys=True)}`.

Hard-linked files are ordinary readable files that share immutable filesystem
storage with the frozen primary release. Copying this package to another
filesystem materializes independent files.

## Submission metadata still required

Author list, affiliations, corresponding-author details, contributions,
funding, acknowledgements, competing interests, and permanent repository
identifiers remain explicit placeholders because those inputs were not
available.
"""
    (output / "README.md").write_text(text, encoding="utf-8")


def build_release(args: argparse.Namespace) -> dict[str, object]:
    primary = args.primary_release.resolve()
    editorial = args.editorial_release.resolve()
    output = args.output.resolve()
    package_root = Path(__file__).resolve().parents[2]
    templates = package_root / "templates"

    required = [
        primary / "03_primary_G0_G3" / "chemical_space_nodes.tsv",
        primary / "03_primary_G0_G3" / "derivation_events.tsv",
        primary / "06_analysis" / "latent_bridge_candidates.tsv",
        primary / "06_analysis" / "known_G0_pair_bridge_summary.tsv",
        primary
        / "05_sensitivity_and_benchmarks"
        / "external_leakage_control"
        / "benchmark_summary.json",
        primary
        / "05_sensitivity_and_benchmarks"
        / "domain_replay_calibration"
        / "benchmark_summary.json",
        templates / "MANUSCRIPT_CORE.md",
        templates / "SUPPLEMENTARY_INFORMATION.md",
        templates / "MATERIALS_AND_METHODS.md",
        templates / "REFERENCES.md",
    ]
    for number in range(1, 5):
        source_release = editorial / MAIN_FIGURE_RELEASES[number]
        stem = MAIN_FIGURE_STEMS[number]
        required.extend(
            [
                source_release / "figures" / f"{stem}.pdf",
                source_release / "figures" / f"{stem}.svg",
                source_release / "figures" / f"{stem}.png",
                source_release / MAIN_CAPTIONS[number],
                source_release / MAIN_AUDITS[number],
            ]
        )
    require_paths(required)
    ensure_empty_output(output)

    for directory in (
        "manuscript",
        "supplementary_information",
        "supplementary_tables",
        "figures/main",
        "figures/supplementary",
        "source_data",
        "data/generation_csv",
        "code",
        "manifests",
        "validation",
    ):
        (output / directory).mkdir(parents=True, exist_ok=True)

    main_image_map, main_audits, main_caption_map = copy_main_figures(
        editorial,
        output,
    )
    supplementary_image_map = copy_supplementary_figures(editorial, output)

    generation_audit = export_generation_csvs(
        primary / "03_primary_G0_G3" / "chemical_space_nodes.tsv",
        output / "data" / "generation_csv",
    )
    table_groups = prepare_supplementary_tables(editorial, primary, output)
    workbook_path = (
        output
        / "supplementary_tables"
        / "Supplementary_Tables_S1-S11.xlsx"
    )
    table_audit = build_supplementary_workbook(
        table_groups,
        workbook_path,
        generation_audit,
        main_audits,
    )

    methods = normalize_embedded_heading(
        read_text(templates / "MATERIALS_AND_METHODS.md"),
        "## Materials and methods",
    )
    references = normalize_embedded_heading(
        read_text(templates / "REFERENCES.md"),
        "## References",
    )
    manuscript_template = read_text(templates / "MANUSCRIPT_CORE.md")
    manuscript_docx_text = manuscript_template.replace(
        "{{METHODS}}",
        methods.strip(),
    ).replace(
        "{{REFERENCES}}",
        references.strip(),
    )
    manuscript_markdown = assemble_markdown(
        manuscript_template,
        {"METHODS": methods, "REFERENCES": references},
        {
            f"FIGURE_{number}": f"../figures/main/Figure_{number}.png"
            for number in range(1, 5)
        },
        main_caption_map,
    )
    manuscript_md = output / "manuscript" / "Main_Manuscript_with_Figures.md"
    manuscript_md.write_text(manuscript_markdown, encoding="utf-8")
    manuscript_docx = output / "manuscript" / "Main_Manuscript_with_Figures.docx"
    markdown_to_docx(
        manuscript_docx_text,
        manuscript_docx,
        main_image_map,
        ARTICLE_TITLE,
        main_caption_map,
    )
    manuscript_pdf = convert_docx_to_pdf(manuscript_docx, output / "manuscript")

    figure_legends_text = "# Figure legends\n\n" + "\n\n".join(
        main_caption_map[f"FIGURE_{number}"] for number in range(1, 5)
    )
    figure_legends_md = output / "manuscript" / "Figure_Legends.md"
    figure_legends_md.write_text(figure_legends_text + "\n", encoding="utf-8")
    figure_legends_docx = output / "manuscript" / "Figure_Legends.docx"
    markdown_to_docx(
        figure_legends_text,
        figure_legends_docx,
        {},
        f"Figure legends: {ARTICLE_TITLE}",
    )
    figure_legends_pdf = convert_docx_to_pdf(
        figure_legends_docx,
        output / "manuscript",
    )

    supplementary_template = read_text(templates / "SUPPLEMENTARY_INFORMATION.md")
    supplementary_markdown = assemble_markdown(
        supplementary_template,
        {},
        {
            f"FIGURE_S{number}": (
                f"../figures/supplementary/Figure_S{number}.png"
            )
            for number in range(1, 9)
        },
        SUPPLEMENTARY_CAPTIONS,
    )
    supplementary_md = (
        output / "supplementary_information" / "Supplementary_Information.md"
    )
    supplementary_md.write_text(supplementary_markdown, encoding="utf-8")
    supplementary_docx = (
        output / "supplementary_information" / "Supplementary_Information.docx"
    )
    markdown_to_docx(
        supplementary_template,
        supplementary_docx,
        supplementary_image_map,
        f"Supplementary Information: {ARTICLE_TITLE}",
        SUPPLEMENTARY_CAPTIONS,
    )
    supplementary_pdf = convert_docx_to_pdf(
        supplementary_docx,
        output / "supplementary_information",
    )

    copy_tree(package_root, output / "code" / "article_release_builder_v0.3.0")
    copy_tree(primary / "code", output / "code" / "primary_G0_G3_workflow")
    for number in range(1, 5):
        copy_tree(
            editorial / MAIN_FIGURE_RELEASES[number] / "workflow",
            output / "code" / f"figure_{number}_workflow",
        )
    copy_tree(
        editorial / "manuscript_v4_figure_redesign" / "workflow",
        output / "code" / "supplementary_figure_workflow",
    )

    raw_status: dict[str, object]
    if args.skip_raw_hardlinks:
        raw_status = {"mode": "skipped_by_request"}
    else:
        raw_status = {
            "mode": "hardlink_snapshot",
            **hardlink_tree(primary, output / "data" / "raw_primary_release"),
        }
    write_release_readme(output, raw_status)

    copy_file(
        primary / "09_release_audit" / "release_file_manifest.tsv",
        output / "manifests" / "RAW_PRIMARY_RELEASE_FILE_MANIFEST.tsv",
    )
    copy_file(
        primary / "09_release_audit" / "release_checksums.sha256",
        output / "manifests" / "RAW_PRIMARY_RELEASE_CHECKSUMS.sha256",
    )

    validation = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "primary_release": str(primary),
        "editorial_release": str(editorial),
        "output": str(output),
        "generation_csvs": generation_audit,
        "supplementary_tables": table_audit,
        "main_figure_audits": [str(path) for path in main_audits],
        "main_manuscript_docx": {
            "path": str(manuscript_docx.relative_to(output)),
            "size_bytes": manuscript_docx.stat().st_size,
            "sha256": sha256_file(manuscript_docx),
        },
        "main_manuscript_pdf": {
            "path": str(manuscript_pdf.relative_to(output)),
            "size_bytes": manuscript_pdf.stat().st_size,
            "sha256": sha256_file(manuscript_pdf),
        },
        "figure_legends_pdf": {
            "path": str(figure_legends_pdf.relative_to(output)),
            "size_bytes": figure_legends_pdf.stat().st_size,
            "sha256": sha256_file(figure_legends_pdf),
        },
        "supplementary_docx": {
            "path": str(supplementary_docx.relative_to(output)),
            "size_bytes": supplementary_docx.stat().st_size,
            "sha256": sha256_file(supplementary_docx),
        },
        "supplementary_pdf": {
            "path": str(supplementary_pdf.relative_to(output)),
            "size_bytes": supplementary_pdf.stat().st_size,
            "sha256": sha256_file(supplementary_pdf),
        },
        "raw_snapshot": raw_status,
        "status": "PASS",
    }
    final_audit = validate_completed_release(
        output,
        primary,
        main_audits,
        table_groups,
    )
    validation["final_audit"] = final_audit
    write_json(output / "validation" / "BUILD_VALIDATION.json", validation)

    package_rows = package_manifest(
        output,
        excluded_prefixes=("data/raw_primary_release/",),
    )
    write_tsv(
        output / "manifests" / "ARTICLE_PACKAGE_FILE_MANIFEST.tsv",
        [
            ["relative_path", "size_bytes", "sha256"],
            *[
                [row["relative_path"], row["size_bytes"], row["sha256"]]
                for row in package_rows
            ],
        ],
    )
    return validation


def main(argv: list[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        validation = build_release(args)
        print(json.dumps(validation, indent=2))
        return 0
    except Exception as exc:
        print(f"BUILD FAILED: {exc}", file=sys.stderr)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
