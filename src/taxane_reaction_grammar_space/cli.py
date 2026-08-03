from __future__ import annotations

import argparse
from pathlib import Path

from .audit import audit_release
from .analyze import analyze_chemical_space
from .benchmark import run_leakage_controlled_benchmark
from .environment import record_environment
from .figures import render_publication_figures
from .generate import generate_chemical_space
from .manuscript import render_manuscript
from .provenance import summarize_rule_library_provenance
from .rules import prepare_generative_grammar, prepare_taxane_domain_grammar
from .screen import screen_grammar_against_seeds
from .sensitivity import compare_g1_sensitivity_spaces
from .study_figures import render_study_figures
from .select import (
    assemble_open_grammar,
    augment_open_grammar_with_domain,
    select_taxane_activated_grammar,
)
from .validate import validate_generated_space


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="taxane-grammar-space",
        description="Build and analyze a reaction-grammar-accessible taxane chemical space.",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare-rules",
        help="Create a strict, semantically nonredundant T1 generative grammar.",
    )
    prepare.add_argument("--rules", required=True, type=Path)
    prepare.add_argument("--output-dir", required=True, type=Path)
    prepare.add_argument("--chunk-size", type=int, default=2000)
    prepare.add_argument("--representatives-per-group", type=int, default=3)
    prepare.add_argument("--max-reactant-atoms", type=int, default=48)
    prepare.add_argument("--max-rules", type=int, default=None)
    prepare.add_argument(
        "--release-tier",
        choices=["T1", "T2", "T3"],
        default="T1",
    )
    prepare.add_argument(
        "--require-single-center",
        action=argparse.BooleanOptionalAction,
        default=True,
    )

    domain = subparsers.add_parser(
        "prepare-domain-rules",
        help="Harmonize reviewed taxane-domain grammar for executable use.",
    )
    domain.add_argument("--rules", required=True, type=Path)
    domain.add_argument("--output-dir", required=True, type=Path)

    screen = subparsers.add_parser(
        "screen-grammar",
        help="Retain grammar productions compatible with at least one G0 taxane.",
    )
    screen.add_argument("--grammar", required=True, type=Path)
    screen.add_argument("--nodes", required=True, type=Path)
    screen.add_argument("--output-dir", required=True, type=Path)
    screen.add_argument("--validate-prefilter-rules", type=int, default=100)
    screen.add_argument("--max-site-matches-per-seed", type=int, default=512)
    screen.add_argument("--matched-seed-id-sample-size", type=int, default=25)

    generate = subparsers.add_parser(
        "generate-space",
        help="Derive an audited G1-G3 chemical space from the G0 taxane library.",
    )
    generate.add_argument("--grammar", required=True, type=Path)
    generate.add_argument("--nodes", required=True, type=Path)
    generate.add_argument("--output-dir", required=True, type=Path)
    generate.add_argument("--max-generation", type=int, default=3)
    generate.add_argument("--max-products-per-parent-rule", type=int, default=256)
    generate.add_argument("--min-source-atom-retention", type=float, default=0.65)
    generate.add_argument("--max-abs-formal-charge", type=int, default=2)
    generate.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=False,
    )

    select = subparsers.add_parser(
        "select-grammar",
        help="Build primary and extended taxane-activated grammar tiers.",
    )
    select.add_argument("--activated", required=True, type=Path)
    select.add_argument("--output-dir", required=True, type=Path)
    select.add_argument("--representatives-per-group", type=int, default=3)
    select.add_argument("--max-heavy-atom-gain", type=int, default=24)
    select.add_argument("--max-product-pattern-growth", type=int, default=32)

    assemble = subparsers.add_parser(
        "assemble-open-grammar",
        help="Combine G0-domain and global representatives for dynamic activation.",
    )
    assemble.add_argument("--global-selected", required=True, type=Path)
    assemble.add_argument("--g0-selected", required=True, type=Path)
    assemble.add_argument("--output-dir", required=True, type=Path)
    assemble.add_argument("--representatives-per-group", type=int, default=4)

    augment = subparsers.add_parser(
        "augment-domain-grammar",
        help="Add reviewed family-domain rules to an external open grammar.",
    )
    augment.add_argument("--open-grammar", required=True, type=Path)
    augment.add_argument("--domain-grammar", required=True, type=Path)
    augment.add_argument("--output-dir", required=True, type=Path)

    benchmark = subparsers.add_parser(
        "benchmark",
        help="Run an exact-pair leakage-controlled Taxol pathway recovery benchmark.",
    )
    benchmark.add_argument("--activated-grammar", required=True, type=Path)
    benchmark.add_argument("--taxol-pathway", required=True, type=Path)
    benchmark.add_argument("--nodes", required=True, type=Path)
    benchmark.add_argument("--output-dir", required=True, type=Path)
    benchmark.add_argument("--decoys-per-positive", type=int, default=20)
    benchmark.add_argument("--max-products-per-rule", type=int, default=256)
    benchmark.add_argument(
        "--exclude-taxol-derived",
        action=argparse.BooleanOptionalAction,
        default=True,
    )

    provenance = subparsers.add_parser(
        "summarize-provenance",
        help="Summarize source-to-rule provenance and evidence retention.",
    )
    provenance.add_argument("--build-root", required=True, type=Path)
    provenance.add_argument("--final-grammar", required=True, type=Path)
    provenance.add_argument("--output-dir", required=True, type=Path)
    provenance.add_argument("--prepared-t1-summary", type=Path, default=None)
    provenance.add_argument("--prepared-t2-summary", type=Path, default=None)
    provenance.add_argument("--prepared-t3-summary", type=Path, default=None)

    sensitivity = subparsers.add_parser(
        "compare-g1-sensitivity",
        help="Compare full-stereochemistry G1 spaces across exclusive evidence tiers.",
    )
    sensitivity.add_argument("--primary-nodes", required=True, type=Path)
    sensitivity.add_argument("--primary-summary", required=True, type=Path)
    sensitivity.add_argument("--t2-nodes", required=True, type=Path)
    sensitivity.add_argument("--t2-summary", required=True, type=Path)
    sensitivity.add_argument("--t3-nodes", required=True, type=Path)
    sensitivity.add_argument("--t3-summary", required=True, type=Path)
    sensitivity.add_argument("--output-dir", required=True, type=Path)

    analyze = subparsers.add_parser(
        "analyze-space",
        help="Compute structure-based G0-G3 chemical-space analyses.",
    )
    analyze.add_argument("--nodes", required=True, type=Path)
    analyze.add_argument("--events", required=True, type=Path)
    analyze.add_argument("--application-audit", required=True, type=Path)
    analyze.add_argument("--rejections", required=True, type=Path)
    analyze.add_argument("--output-dir", required=True, type=Path)
    analyze.add_argument("--projection-max-nodes-per-generation", type=int, default=20000)
    analyze.add_argument(
        "--similarity-max-nodes-per-generation", type=int, default=50_000
    )
    analyze.add_argument("--random-seed", type=int, default=1729)

    validate = subparsers.add_parser(
        "validate-space",
        help="Audit generated-space identity and generation invariants.",
    )
    validate.add_argument("--database", required=True, type=Path)
    validate.add_argument("--output-dir", required=True, type=Path)

    figures = subparsers.add_parser(
        "render-figures",
        help="Render publication figures and a chart provenance map.",
    )
    figures.add_argument("--analysis-dir", required=True, type=Path)
    figures.add_argument("--output-dir", required=True, type=Path)
    figures.add_argument("--grammar-summary", type=Path, default=None)
    figures.add_argument("--screen-summary", type=Path, default=None)
    figures.add_argument("--selection-summary", type=Path, default=None)
    figures.add_argument("--selected-grammar", type=Path, default=None)
    figures.add_argument("--benchmark-dir", type=Path, default=None)

    study_figures = subparsers.add_parser(
        "render-study-figures",
        help="Render publication figures with panel-level source data.",
    )
    study_figures.add_argument("--analysis-dir", required=True, type=Path)
    study_figures.add_argument("--provenance-dir", required=True, type=Path)
    study_figures.add_argument("--sensitivity-dir", required=True, type=Path)
    study_figures.add_argument(
        "--external-benchmark-dir", required=True, type=Path
    )
    study_figures.add_argument(
        "--domain-benchmark-dir", required=True, type=Path
    )
    study_figures.add_argument(
        "--t2-analysis-dir", required=True, type=Path
    )
    study_figures.add_argument(
        "--t3-analysis-dir", required=True, type=Path
    )
    study_figures.add_argument("--output-dir", required=True, type=Path)

    manuscript = subparsers.add_parser(
        "render-manuscript",
        help="Populate the English manuscript from validated G0-G3 analyses.",
    )
    manuscript.add_argument("--template", required=True, type=Path)
    manuscript.add_argument("--analysis-dir", required=True, type=Path)
    manuscript.add_argument("--output-dir", required=True, type=Path)

    audit = subparsers.add_parser(
        "audit-release",
        help="Verify the complete G0-G3 study release and write checksums.",
    )
    audit.add_argument("--code-dir", required=True, type=Path)
    audit.add_argument("--provenance-dir", required=True, type=Path)
    audit.add_argument("--space-dir", required=True, type=Path)
    audit.add_argument("--validation-dir", required=True, type=Path)
    audit.add_argument("--analysis-dir", required=True, type=Path)
    audit.add_argument("--figures-dir", required=True, type=Path)
    audit.add_argument("--manuscript-dir", required=True, type=Path)
    audit.add_argument("--output-dir", required=True, type=Path)

    environment = subparsers.add_parser(
        "record-environment",
        help="Record software versions and verify configured input hashes.",
    )
    environment.add_argument("--output-dir", required=True, type=Path)
    environment.add_argument("--study-config", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "prepare-rules":
        paths = prepare_generative_grammar(
            args.rules,
            args.output_dir,
            chunk_size=args.chunk_size,
            representatives_per_group=args.representatives_per_group,
            max_reactant_atoms=args.max_reactant_atoms,
            require_single_center=args.require_single_center,
            max_rules=args.max_rules,
            release_tier=args.release_tier,
        )
        print("Generative grammar preparation complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "prepare-domain-rules":
        paths = prepare_taxane_domain_grammar(args.rules, args.output_dir)
        print("Taxane-domain grammar harmonization complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "screen-grammar":
        paths = screen_grammar_against_seeds(
            args.grammar,
            args.nodes,
            args.output_dir,
            validate_prefilter_rules=args.validate_prefilter_rules,
            max_site_matches_per_seed=args.max_site_matches_per_seed,
            matched_seed_id_sample_size=args.matched_seed_id_sample_size,
        )
        print("G0 grammar compatibility screen complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "generate-space":
        paths = generate_chemical_space(
            args.grammar,
            args.nodes,
            args.output_dir,
            max_generation=args.max_generation,
            max_products_per_parent_rule=args.max_products_per_parent_rule,
            min_source_atom_retention=args.min_source_atom_retention,
            max_abs_formal_charge=args.max_abs_formal_charge,
            resume=args.resume,
        )
        print("Taxane reaction-grammar chemical-space derivation complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "select-grammar":
        paths = select_taxane_activated_grammar(
            args.activated,
            args.output_dir,
            representatives_per_group=args.representatives_per_group,
            max_heavy_atom_gain=args.max_heavy_atom_gain,
            max_product_pattern_growth=args.max_product_pattern_growth,
        )
        print("Taxane-activated grammar selection complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "assemble-open-grammar":
        paths = assemble_open_grammar(
            args.global_selected,
            args.g0_selected,
            args.output_dir,
            representatives_per_group=args.representatives_per_group,
        )
        print("Domain-seeded open grammar assembly complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "augment-domain-grammar":
        paths = augment_open_grammar_with_domain(
            args.open_grammar,
            args.domain_grammar,
            args.output_dir,
        )
        print("Reviewed domain grammar augmentation complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "benchmark":
        paths = run_leakage_controlled_benchmark(
            args.activated_grammar,
            args.taxol_pathway,
            args.nodes,
            args.output_dir,
            decoys_per_positive=args.decoys_per_positive,
            max_products_per_rule=args.max_products_per_rule,
            exclude_taxol_derived=args.exclude_taxol_derived,
        )
        print("Leakage-controlled pathway benchmark complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "summarize-provenance":
        paths = summarize_rule_library_provenance(
            args.build_root,
            args.final_grammar,
            args.output_dir,
            prepared_t1_summary=args.prepared_t1_summary,
            prepared_t2_summary=args.prepared_t2_summary,
            prepared_t3_summary=args.prepared_t3_summary,
        )
        print("Rule-library provenance summary complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "compare-g1-sensitivity":
        paths = compare_g1_sensitivity_spaces(
            args.primary_nodes,
            args.primary_summary,
            args.t2_nodes,
            args.t2_summary,
            args.t3_nodes,
            args.t3_summary,
            args.output_dir,
        )
        print("G1 evidence-layer sensitivity comparison complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "analyze-space":
        paths = analyze_chemical_space(
            args.nodes,
            args.events,
            args.application_audit,
            args.rejections,
            args.output_dir,
            projection_max_nodes_per_generation=(
                args.projection_max_nodes_per_generation
            ),
            similarity_max_nodes_per_generation=(
                args.similarity_max_nodes_per_generation
            ),
            random_seed=args.random_seed,
        )
        print("Structure-based chemical-space analysis complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "validate-space":
        paths = validate_generated_space(args.database, args.output_dir)
        print("Generated-space validation complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "render-figures":
        paths = render_publication_figures(
            args.analysis_dir,
            args.output_dir,
            grammar_summary_path=args.grammar_summary,
            screen_summary_path=args.screen_summary,
            selection_summary_path=args.selection_summary,
            selected_grammar_path=args.selected_grammar,
            benchmark_dir=args.benchmark_dir,
        )
        print("Publication figure rendering complete.")
        for name, outputs in paths.items():
            print(f"{name}: {', '.join(str(path) for path in outputs)}")
    elif args.command == "render-study-figures":
        result = render_study_figures(
            args.analysis_dir,
            args.provenance_dir,
            args.sensitivity_dir,
            args.external_benchmark_dir,
            args.domain_benchmark_dir,
            args.t2_analysis_dir,
            args.t3_analysis_dir,
            args.output_dir,
        )
        print("Publication figures and panel source data complete.")
        print(f"source_manifest: {result['source_manifest']}")
        print(f"summary: {result['summary']}")
    elif args.command == "render-manuscript":
        paths = render_manuscript(
            args.template,
            args.analysis_dir,
            args.output_dir,
        )
        print("Final English manuscript rendering complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "audit-release":
        paths = audit_release(
            args.code_dir,
            args.provenance_dir,
            args.space_dir,
            args.validation_dir,
            args.analysis_dir,
            args.figures_dir,
            args.manuscript_dir,
            args.output_dir,
        )
        print("Complete G0-G3 release audit passed.")
        for name, path in paths.items():
            print(f"{name}: {path}")
    elif args.command == "record-environment":
        paths = record_environment(
            args.output_dir,
            study_config=args.study_config,
        )
        print("Software environment and input snapshot audit complete.")
        for name, path in paths.items():
            print(f"{name}: {path}")


if __name__ == "__main__":
    main()
