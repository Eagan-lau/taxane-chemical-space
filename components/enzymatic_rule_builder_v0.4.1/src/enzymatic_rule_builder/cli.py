from __future__ import annotations

import argparse
import json
from pathlib import Path

from .inventory import write_inventory
from .participant_registry import build_participant_registry_from_manifest
from .pipeline import run_build_all, run_build_from_raw
from .taxol_multisubstrate import augment_multisubstrate_taxol_release
from . import __version__


def _quiet_rdkit_logs() -> None:
    """Suppress repetitive RDKit parser diagnostics without hiding Python exceptions."""
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass


def _split_layers(text: str) -> list[str]:
    return [x.strip() for x in str(text).replace(";", ",").split(",") if x.strip()]


def _add_anchor_generalization_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--generalize-exact-anchors",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run a separate RXNMapper/RDChiral exact-anchor generalization flow and write reaction_smarts_library.anchor_derived.tsv; does not mix anchor-derived SMARTS into core.",
    )
    parser.add_argument(
        "--anchor-generalization-layers",
        default="T1_Bio_Core,T2_Bio_Extended,T3_Chem_like",
        help="Comma/semicolon separated evidence layers eligible for the separate anchor_generalization flow.",
    )
    parser.add_argument(
        "--anchor-generalization-mode",
        choices=["strict", "permissive", "both"],
        default="strict",
        help="strict keeps only replay-validated exact-anchor SMARTS; permissive keeps all extractable SMARTS with QC warnings; both writes both releases.",
    )
    parser.add_argument("--anchor-generalization-batch-size", type=int, default=16)
    parser.add_argument("--anchor-generalization-min-mapper-confidence", type=float, default=0.50)
    parser.add_argument("--anchor-generalization-max-records", type=int, default=None)
    parser.add_argument("--anchor-generalization-max-reaction-chars", type=int, default=2000)
    parser.add_argument("--anchor-generalization-max-mapper-tokens", type=int, default=500)
    parser.add_argument("--anchor-generalization-mapper-timeout-seconds", type=int, default=120)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzymatic-rules",
        description="Build a database-derived directional reaction SMARTS rule library",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build-rules", help="Build a general transformation rule library from a source manifest")
    build.add_argument("--manifest", required=True, help="YAML manifest describing reaction/template datasets")
    build.add_argument("--output-root", required=True, help="Output directory")
    build.add_argument("--taxol-pathway", default=None, help="Optional override path for taxol_pathway_csv dataset")
    build.add_argument("--cofactor-yaml", default=None, help="Optional participant/cofactor registry YAML; preferably database-derived")
    build.add_argument("--derive-participant-registry", action=argparse.BooleanOptionalAction, default=True, help="Derive a participant-role registry from input reactions when --cofactor-yaml is not supplied")
    build.add_argument("--participant-registry-min-frequency", type=int, default=3)
    build.add_argument("--participant-registry-max-heavy-atoms", type=int, default=80)
    build.add_argument("--scoring-yaml", default=None, help="Optional scoring weights YAML")
    build.add_argument("--allow-exact-pairs-as-predictive-rules", action="store_true", help="Dangerous: promote exact pairs to predictive rules; off by default")
    build.add_argument("--abstract-exact-reactions", action=argparse.BooleanOptionalAction, default=None, help="Use RXNMapper/RDChiral to generalize selected exact reactions into predictive SMARTS. Default: auto-enabled when curated Taxol exact anchors are present.")
    build.add_argument("--exact-abstraction-layers", default="T1_Bio_Core,T2_Bio_Extended", help="Comma/semicolon separated evidence layers eligible for exact-reaction abstraction")
    build.add_argument("--exact-abstraction-batch-size", type=int, default=16)
    build.add_argument("--exact-abstraction-min-mapper-confidence", type=float, default=0.50)
    build.add_argument("--exact-abstraction-max-records", type=int, default=None)
    _add_anchor_generalization_args(build)
    build.add_argument("--max-decoys", type=int, default=500)
    build.add_argument("--skip-benchmark", action="store_true", help="Skip recall/decoy benchmark files after release generation; writes empty benchmark tables with a skipped summary")
    build.add_argument("--require-smarts-rules", action="store_true", help="Fail if no predictive reaction SMARTS rules are produced in the all tier")
    build.add_argument("--require-core-rules", action="store_true", help="Fail if the strict-core SMARTS release used for downstream network construction is empty")
    build.add_argument("--transfer-family-consensus", action=argparse.BooleanOptionalAction, default=True, help="Generate replay-validated donor-transfer family SMARTS from role-aware main substrate-product projections")
    build.add_argument("--transfer-family-mcs-timeout", type=int, default=10)
    build.add_argument("--data-driven-consensus", action=argparse.BooleanOptionalAction, default=True, help="Cluster all QC-passing SMARTS across databases and promote high-support consensus representative rules")
    build.add_argument("--consensus-min-evidence-rows", type=int, default=3)
    build.add_argument("--consensus-min-source-database-count", type=int, default=1)
    build.add_argument("--consensus-min-template-count", type=int, default=2)

    raw = sub.add_parser(
        "build-from-raw",
        help="Discover raw external database files and build a v0.4.1-native SMARTS rule library without network_ready inputs",
    )
    raw.add_argument("--external-db-root", required=True, help="Raw external database root directory")
    raw.add_argument("--output-root", required=True, help="Output directory")
    raw.add_argument("--taxol-pathway", default=None, help="Optional Taxol pathway anchor CSV")
    raw.add_argument("--cofactor-yaml", default=None, help="Optional participant/cofactor registry YAML; preferably database-derived")
    raw.add_argument("--derive-participant-registry", action=argparse.BooleanOptionalAction, default=True, help="Derive a participant-role registry from input reactions when --cofactor-yaml is not supplied")
    raw.add_argument("--participant-registry-min-frequency", type=int, default=3)
    raw.add_argument("--participant-registry-max-heavy-atoms", type=int, default=80)
    raw.add_argument("--scoring-yaml", default=None, help="Optional scoring weights YAML")
    raw.add_argument("--include-taxol-anchors", action=argparse.BooleanOptionalAction, default=None, help="Include curated Taxol exact anchors. Default: auto-enabled when --taxol-pathway is supplied")
    raw.add_argument("--include-uspto", action=argparse.BooleanOptionalAction, default=True, help="Include BioNavi USPTO_NPL as T3 exploratory evidence")
    raw.add_argument("--include-annotation-only", action=argparse.BooleanOptionalAction, default=True, help="Include KEGG/MetaNetX annotation-only evidence")
    raw.add_argument("--use-processed-bionavi", action=argparse.BooleanOptionalAction, default=True, help="Use prevalidated BioNavi processed reaction-SMILES lines when present")
    raw.add_argument("--max-retrorules", type=int, default=None)
    raw.add_argument("--max-rhea", type=int, default=None)
    raw.add_argument("--max-bionavi-per-file", type=int, default=None)
    raw.add_argument("--max-uspto", type=int, default=None)
    raw.add_argument("--max-kegg", type=int, default=None)
    raw.add_argument("--max-metanetx", type=int, default=None)
    raw.add_argument("--allow-exact-pairs-as-predictive-rules", action="store_true", help="Dangerous: promote exact pairs to predictive rules; off by default")
    raw.add_argument("--abstract-exact-reactions", action=argparse.BooleanOptionalAction, default=None, help="Use RXNMapper/RDChiral to generalize selected exact reactions into predictive SMARTS. Default: auto-enabled when Taxol anchors are included.")
    raw.add_argument("--exact-abstraction-layers", default="T1_Bio_Core,T2_Bio_Extended", help="Comma/semicolon separated evidence layers eligible for exact-reaction abstraction")
    raw.add_argument("--exact-abstraction-batch-size", type=int, default=16)
    raw.add_argument("--exact-abstraction-min-mapper-confidence", type=float, default=0.50)
    raw.add_argument("--exact-abstraction-max-records", type=int, default=None)
    _add_anchor_generalization_args(raw)
    raw.add_argument("--max-decoys", type=int, default=500)
    raw.add_argument("--skip-benchmark", action="store_true", help="Skip recall/decoy benchmark files after release generation; writes empty benchmark tables with a skipped summary")
    raw.add_argument("--require-smarts-rules", action="store_true", help="Fail if no predictive reaction SMARTS rules are produced in the all tier")
    raw.add_argument("--require-core-rules", action="store_true", help="Fail if the strict-core SMARTS release used for downstream network construction is empty")
    raw.add_argument("--transfer-family-consensus", action=argparse.BooleanOptionalAction, default=True, help="Generate replay-validated donor-transfer family SMARTS from role-aware main substrate-product projections")
    raw.add_argument("--transfer-family-mcs-timeout", type=int, default=10)
    raw.add_argument("--data-driven-consensus", action=argparse.BooleanOptionalAction, default=True, help="Cluster all QC-passing SMARTS across databases and promote high-support consensus representative rules")
    raw.add_argument("--consensus-min-evidence-rows", type=int, default=3)
    raw.add_argument("--consensus-min-source-database-count", type=int, default=1)
    raw.add_argument("--consensus-min-template-count", type=int, default=2)

    preg = sub.add_parser("build-participant-registry", help="Build a database-derived participant/cofactor registry from a manifest")
    preg.add_argument("--manifest", required=True, help="YAML manifest describing reaction/template datasets")
    preg.add_argument("--output-dir", required=True, help="Output directory")
    preg.add_argument("--min-frequency", type=int, default=3)
    preg.add_argument("--max-heavy-atoms", type=int, default=80)

    inv = sub.add_parser("inventory", help="Write checksum inventory for input files")
    inv.add_argument("--paths", nargs="+", required=True)
    inv.add_argument("--output-dir", required=True)

    aug = sub.add_parser(
        "augment-taxol-multisubstrate",
        help="Add explicit multi-substrate rule annotations and Taxol inferred external participant variants to an existing release",
    )
    aug.add_argument("--input-library", required=True, help="Existing reaction_smarts_library.*.tsv release")
    aug.add_argument("--rule-qc", required=True, help="Rule compile QC table with n_reactants/n_products/status")
    aug.add_argument("--taxol-pathway", required=True, help="Taxol known pathway CSV with Enzyme/Substrate/Product/EC")
    aug.add_argument("--output-dir", required=True, help="Output directory for augmented release files")
    aug.add_argument("--chunk-size", type=int, default=50000)

    return parser


def main(argv: list[str] | None = None) -> int:
    _quiet_rdkit_logs()
    args = build_parser().parse_args(argv)
    if args.command == "build-rules":
        summary = run_build_all(
            manifest=args.manifest,
            output_root=args.output_root,
            taxol_pathway=args.taxol_pathway,
            cofactor_yaml=args.cofactor_yaml,
            derive_participant_registry=args.derive_participant_registry,
            participant_registry_min_frequency=args.participant_registry_min_frequency,
            participant_registry_max_heavy_atoms=args.participant_registry_max_heavy_atoms,
            scoring_yaml=args.scoring_yaml,
            allow_exact_pairs_as_predictive_rules=args.allow_exact_pairs_as_predictive_rules,
            abstract_exact_reactions=args.abstract_exact_reactions,
            exact_abstraction_layers=_split_layers(args.exact_abstraction_layers),
            exact_abstraction_batch_size=args.exact_abstraction_batch_size,
            exact_abstraction_min_mapper_confidence=args.exact_abstraction_min_mapper_confidence,
            exact_abstraction_max_records=args.exact_abstraction_max_records,
            generalize_exact_anchors=args.generalize_exact_anchors,
            anchor_generalization_mode=args.anchor_generalization_mode,
            anchor_generalization_layers=_split_layers(args.anchor_generalization_layers),
            anchor_generalization_batch_size=args.anchor_generalization_batch_size,
            anchor_generalization_min_mapper_confidence=args.anchor_generalization_min_mapper_confidence,
            anchor_generalization_max_records=args.anchor_generalization_max_records,
            anchor_generalization_max_reaction_chars=args.anchor_generalization_max_reaction_chars,
            anchor_generalization_max_mapper_tokens=args.anchor_generalization_max_mapper_tokens,
            anchor_generalization_mapper_timeout_seconds=args.anchor_generalization_mapper_timeout_seconds,
            max_decoys=args.max_decoys,
            skip_benchmark=args.skip_benchmark,
            require_smarts_rules=args.require_smarts_rules,
            require_core_rules=args.require_core_rules,
            transfer_family_consensus=args.transfer_family_consensus,
            transfer_family_mcs_timeout=args.transfer_family_mcs_timeout,
            data_driven_consensus=args.data_driven_consensus,
            consensus_min_evidence_rows=args.consensus_min_evidence_rows,
            consensus_min_source_database_count=args.consensus_min_source_database_count,
            consensus_min_template_count=args.consensus_min_template_count,
        )
        print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))
        return 0
    if args.command == "build-from-raw":
        summary = run_build_from_raw(
            external_db_root=args.external_db_root,
            output_root=args.output_root,
            taxol_pathway=args.taxol_pathway,
            cofactor_yaml=args.cofactor_yaml,
            derive_participant_registry=args.derive_participant_registry,
            participant_registry_min_frequency=args.participant_registry_min_frequency,
            participant_registry_max_heavy_atoms=args.participant_registry_max_heavy_atoms,
            scoring_yaml=args.scoring_yaml,
            include_taxol_anchors=args.include_taxol_anchors,
            include_uspto=args.include_uspto,
            include_annotation_only=args.include_annotation_only,
            use_processed_bionavi=args.use_processed_bionavi,
            max_retrorules=args.max_retrorules,
            max_rhea=args.max_rhea,
            max_bionavi_per_file=args.max_bionavi_per_file,
            max_uspto=args.max_uspto,
            max_kegg=args.max_kegg,
            max_metanetx=args.max_metanetx,
            allow_exact_pairs_as_predictive_rules=args.allow_exact_pairs_as_predictive_rules,
            abstract_exact_reactions=args.abstract_exact_reactions,
            exact_abstraction_layers=_split_layers(args.exact_abstraction_layers),
            exact_abstraction_batch_size=args.exact_abstraction_batch_size,
            exact_abstraction_min_mapper_confidence=args.exact_abstraction_min_mapper_confidence,
            exact_abstraction_max_records=args.exact_abstraction_max_records,
            generalize_exact_anchors=args.generalize_exact_anchors,
            anchor_generalization_mode=args.anchor_generalization_mode,
            anchor_generalization_layers=_split_layers(args.anchor_generalization_layers),
            anchor_generalization_batch_size=args.anchor_generalization_batch_size,
            anchor_generalization_min_mapper_confidence=args.anchor_generalization_min_mapper_confidence,
            anchor_generalization_max_records=args.anchor_generalization_max_records,
            anchor_generalization_max_reaction_chars=args.anchor_generalization_max_reaction_chars,
            anchor_generalization_max_mapper_tokens=args.anchor_generalization_max_mapper_tokens,
            anchor_generalization_mapper_timeout_seconds=args.anchor_generalization_mapper_timeout_seconds,
            max_decoys=args.max_decoys,
            skip_benchmark=args.skip_benchmark,
            require_smarts_rules=args.require_smarts_rules,
            require_core_rules=args.require_core_rules,
            transfer_family_consensus=args.transfer_family_consensus,
            transfer_family_mcs_timeout=args.transfer_family_mcs_timeout,
            data_driven_consensus=args.data_driven_consensus,
            consensus_min_evidence_rows=args.consensus_min_evidence_rows,
            consensus_min_source_database_count=args.consensus_min_source_database_count,
            consensus_min_template_count=args.consensus_min_template_count,
        )
        print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))
        return 0
    if args.command == "build-participant-registry":
        summary = build_participant_registry_from_manifest(
            args.manifest,
            args.output_dir,
            min_occurrence=args.min_frequency,
            max_heavy_atoms=args.max_heavy_atoms,
        )
        print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))
        return 0
    if args.command == "inventory":
        df = write_inventory(args.paths, args.output_dir)
        print(json.dumps({"files": len(df), "output_dir": str(Path(args.output_dir).resolve())}, indent=2))
        return 0
    if args.command == "augment-taxol-multisubstrate":
        summary = augment_multisubstrate_taxol_release(
            input_library=args.input_library,
            rule_qc=args.rule_qc,
            taxol_pathway=args.taxol_pathway,
            output_dir=args.output_dir,
            chunk_size=args.chunk_size,
        )
        print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
