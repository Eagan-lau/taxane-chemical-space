from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from .abstraction import abstract_exact_reactions_to_smarts
from .anchor_generalization import generalize_exact_anchor_templates
from .anchors import build_anchor_edges
from .benchmark import decoy_pair_benchmark, known_pathway_recall
from .consensus import build_data_driven_consensus_rules
from .inventory import write_inventory
from .normalization import add_main_pairs
from .participant_registry import build_participant_registry_from_sources, registry_to_yaml_payload
from .raw import build_raw_source_manifest
from .release import write_smarts_releases
from .rules import build_rule_library
from .scoring import DEFAULT_WEIGHTS
from .sources import family_evidence_paths, load_datasets, load_family_evidence
from .templates import build_templates, deduplicate_templates, qc_templates
from .transfer_family import build_transfer_family_consensus_templates
from .utils import ensure_dir, read_yaml, write_json, write_yaml


def package_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_config(name: str) -> Path:
    candidates = [package_root() / "configs" / name, Path(__file__).resolve().parents[3] / "configs" / name, Path.cwd() / "configs" / name]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"Cannot find default config {name}; searched {candidates}")


def _bool_mask(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([False] * len(df), index=df.index)
    return df[col].astype(str).str.lower().isin(["true", "1", "yes"])


def run_build_all(
    *,
    manifest: str | Path,
    output_root: str | Path,
    taxol_pathway: str | Path | None = None,
    cofactor_yaml: str | Path | None = None,
    derive_participant_registry: bool = True,
    participant_registry_min_frequency: int = 3,
    participant_registry_max_heavy_atoms: int = 80,
    scoring_yaml: str | Path | None = None,
    allow_exact_pairs_as_predictive_rules: bool = False,
    abstract_exact_reactions: bool | None = None,
    exact_abstraction_layers: list[str] | None = None,
    exact_abstraction_batch_size: int = 16,
    exact_abstraction_min_mapper_confidence: float = 0.50,
    exact_abstraction_max_records: int | None = None,
    generalize_exact_anchors: bool = False,
    anchor_generalization_mode: str = "strict",
    anchor_generalization_layers: list[str] | None = None,
    anchor_generalization_batch_size: int = 16,
    anchor_generalization_min_mapper_confidence: float = 0.50,
    anchor_generalization_max_records: int | None = None,
    anchor_generalization_max_reaction_chars: int = 2000,
    anchor_generalization_max_mapper_tokens: int = 500,
    anchor_generalization_mapper_timeout_seconds: int = 120,
    max_decoys: int = 500,
    skip_benchmark: bool = False,
    require_smarts_rules: bool = False,
    require_core_rules: bool = False,
    transfer_family_consensus: bool = True,
    transfer_family_mcs_timeout: int = 10,
    data_driven_consensus: bool = True,
    consensus_min_evidence_rows: int = 3,
    consensus_min_source_database_count: int = 1,
    consensus_min_template_count: int = 2,
) -> dict:
    out = ensure_dir(output_root)
    dirs = {
        "inventory": ensure_dir(out / "00_inventory"),
        "sources": ensure_dir(out / "01_sources"),
        "templates": ensure_dir(out / "02_templates"),
        "rules": ensure_dir(out / "03_rules"),
        "release": ensure_dir(out / "04_release"),
        "benchmark": ensure_dir(out / "05_benchmark"),
    }
    shutil.copyfile(manifest, dirs["inventory"] / "build_manifest.resolved.yaml")

    scoring = DEFAULT_WEIGHTS
    if scoring_yaml:
        scoring = read_yaml(scoring_yaml)

    source_df, used_paths = load_datasets(manifest, override_taxol_path=taxol_pathway)
    source_df.to_csv(dirs["sources"] / "source_reactions.normalized.tsv", sep="	", index=False)

    participant_registry_path = str(cofactor_yaml) if cofactor_yaml else ""
    participant_registry_summary = {"derived": False, "path": participant_registry_path, "records": 0}
    if not participant_registry_path and derive_participant_registry:
        reg_df = build_participant_registry_from_sources(
            source_df,
            min_occurrence=participant_registry_min_frequency,
            max_heavy_atoms=participant_registry_max_heavy_atoms,
        )
        registry_tsv = dirs["sources"] / "participant_role_registry.tsv"
        registry_yaml = dirs["sources"] / "participant_role_registry.yaml"
        reg_df.to_csv(registry_tsv, sep="\t", index=False)
        write_yaml(registry_yaml, registry_to_yaml_payload(reg_df))
        participant_registry_path = str(registry_yaml)
        participant_registry_summary = {
            "derived": True,
            "path": participant_registry_path,
            "tsv_path": str(registry_tsv),
            "records": int(len(reg_df)),
            "min_occurrence": int(participant_registry_min_frequency),
            "max_heavy_atoms": int(participant_registry_max_heavy_atoms),
        }

    inventory_paths = list(used_paths) + family_evidence_paths(manifest)
    if participant_registry_path:
        inventory_paths.append(str(participant_registry_path))
    if scoring_yaml:
        inventory_paths.append(str(scoring_yaml))
    inventory_paths.append(str(manifest))
    if inventory_paths:
        write_inventory(inventory_paths, dirs["inventory"])

    family_evidence = load_family_evidence(manifest)
    family_evidence.to_csv(dirs["sources"] / "family_evidence.normalized.tsv", sep="\t", index=False)

    main_df = add_main_pairs(source_df, participant_registry_path or None)
    if abstract_exact_reactions is None:
        # Production default: curated Taxol/T1 exact reactions are high-quality
        # sources and should be attempted for exact-to-SMARTS abstraction when
        # present. Missing optional mapper dependencies are handled gracefully by
        # the abstraction module and reported in the summary.
        has_curated_taxol_exact = bool(source_df.get("curated_taxol_anchor", pd.Series(dtype=str)).astype(str).str.lower().isin(["true", "1", "yes"]).any()) if len(source_df) else False
        abstract_exact_reactions = has_curated_taxol_exact and not generalize_exact_anchors
    exact_abstraction_summary = {
        "enabled": bool(abstract_exact_reactions),
        "candidate_records": 0,
        "attempted_records": 0,
        "successful_templates": 0,
        "failed_records": 0,
    }
    if abstract_exact_reactions:
        main_df, _abstraction_report, exact_abstraction_summary = abstract_exact_reactions_to_smarts(
            main_df,
            dirs["sources"],
            enabled_layers=exact_abstraction_layers,
            batch_size=exact_abstraction_batch_size,
            min_mapper_confidence=exact_abstraction_min_mapper_confidence,
            max_records=exact_abstraction_max_records,
        )
        exact_abstraction_summary["enabled"] = True
    main_df.to_csv(dirs["sources"] / "source_reactions.main_pair.tsv", sep="\t", index=False)

    raw_templates = build_templates(main_df, allow_exact_pairs_as_predictive_rules=allow_exact_pairs_as_predictive_rules)
    transfer_family_consensus_summary = {
        "enabled": bool(transfer_family_consensus),
        "candidate_rows": 0,
        "released_template_rows": 0,
        "unique_reaction_smarts": 0,
    }
    if transfer_family_consensus:
        transfer_templates, _transfer_report, transfer_family_consensus_summary = build_transfer_family_consensus_templates(
            main_df,
            output_dir=dirs["rules"],
            mcs_timeout=transfer_family_mcs_timeout,
        )
        if len(transfer_templates):
            raw_templates = pd.concat([raw_templates, transfer_templates], ignore_index=True)
    raw_templates.to_csv(dirs["templates"] / "templates.raw.tsv", sep="\t", index=False)
    qc = qc_templates(raw_templates)
    qc.to_csv(dirs["templates"] / "templates.qc.tsv", sep="\t", index=False)
    dedup = deduplicate_templates(qc)
    dedup.to_csv(dirs["templates"] / "templates.deduplicated.tsv", sep="\t", index=False)

    rules = build_rule_library(dedup, family_evidence, evidence_weights=scoring)
    base_rule_count = int(len(rules))
    data_driven_consensus_summary = {
        "enabled": bool(data_driven_consensus),
        "candidate_groups": 0,
        "promoted_consensus_rules": 0,
    }
    if data_driven_consensus:
        consensus_rules, _consensus_report, data_driven_consensus_summary = build_data_driven_consensus_rules(
            rules,
            output_dir=dirs["rules"],
            min_evidence_rows=consensus_min_evidence_rows,
            min_source_database_count=consensus_min_source_database_count,
            min_template_count=consensus_min_template_count,
        )
        if len(consensus_rules):
            rules = pd.concat([rules, consensus_rules], ignore_index=True)
    # Audit table: includes exact anchors, generalized templates, and promoted consensus rules. Not for direct network construction.
    rules.to_csv(dirs["rules"] / "general_transformation_rules.audit.tsv", sep="\t", index=False)
    # Backward-compatible audit alias.
    rules.to_csv(dirs["rules"] / "general_transformation_rules.all.tsv", sep="\t", index=False)

    # Final predictive release: reaction SMARTS only.
    smarts_release_summary = write_smarts_releases(rules, dirs["release"])

    anchors = build_anchor_edges(dedup)
    anchors.to_csv(dirs["release"] / "curated_exact_anchor_edges.tsv", sep="\t", index=False)

    anchor_generalization_summary = {
        "enabled": bool(generalize_exact_anchors),
        "candidate_anchor_records": 0,
        "attempted_records": 0,
        "successful_templates": 0,
        "failed_records": 0,
        "anchor_derived_released_smarts": 0,
    }
    if generalize_exact_anchors:
        requested_mode = str(anchor_generalization_mode or "strict").strip().lower()
        if requested_mode not in {"strict", "permissive", "both"}:
            raise ValueError(f"Unknown anchor_generalization_mode: {anchor_generalization_mode!r}")
        modes = ["strict", "permissive"] if requested_mode == "both" else [requested_mode]
        mode_summaries = {}
        for mode in modes:
            mode_summaries[mode] = generalize_exact_anchor_templates(
                dedup,
                family_evidence,
                dirs["release"],
                enabled_layers=anchor_generalization_layers,
                batch_size=anchor_generalization_batch_size,
                min_mapper_confidence=anchor_generalization_min_mapper_confidence,
                max_records=anchor_generalization_max_records,
                max_reaction_chars=anchor_generalization_max_reaction_chars,
                max_mapper_tokens=anchor_generalization_max_mapper_tokens,
                mapper_timeout_seconds=anchor_generalization_mapper_timeout_seconds,
                release_mode=mode,
            )
        if len(mode_summaries) == 1:
            anchor_generalization_summary = next(iter(mode_summaries.values()))
            anchor_generalization_summary["mode_summaries"] = mode_summaries
        else:
            anchor_generalization_summary = {
                "enabled": True,
                "release_mode": requested_mode,
                "mode_summaries": mode_summaries,
                "candidate_anchor_records": max(int(s.get("candidate_anchor_records", 0)) for s in mode_summaries.values()),
                "selected_anchor_records": max(int(s.get("selected_anchor_records", 0)) for s in mode_summaries.values()),
                "attempted_records": sum(int(s.get("attempted_records", 0)) for s in mode_summaries.values()),
                "successful_templates": sum(int(s.get("successful_templates", 0)) for s in mode_summaries.values()),
                "failed_records": sum(int(s.get("failed_records", 0)) for s in mode_summaries.values()),
                "strict_anchor_derived_released_smarts": int(mode_summaries.get("strict", {}).get("anchor_derived_released_smarts", 0)),
                "permissive_anchor_derived_released_smarts": int(mode_summaries.get("permissive", {}).get("anchor_derived_released_smarts", 0)),
            }

    benchmark_summary = {"enabled": not bool(skip_benchmark), "skipped": bool(skip_benchmark)}
    if skip_benchmark:
        recall = pd.DataFrame(columns=[
            "known_reaction_id", "enzyme_name", "known_ec", "canonical_known_reaction_smiles",
            "recovered_by_external_generalized_rule", "best_rule_id", "best_rule_confidence",
            "all_hit_rule_ids", "leakage_control",
        ])
        decoys = pd.DataFrame(columns=[
            "decoy_id", "substrate_reaction_id", "product_reaction_id",
            "matched_rule_count", "max_rule_confidence", "matched_rule_ids",
        ])
        benchmark_summary["reason"] = "skipped_by_user_request"
    else:
        recall = known_pathway_recall(rules, source_df)
        decoys = decoy_pair_benchmark(rules, source_df, max_decoys=max_decoys)
    recall.to_csv(dirs["benchmark"] / "known_taxol_pathway_recall.external_generalized_only.tsv", sep="\t", index=False)
    decoys.to_csv(dirs["benchmark"] / "decoy_pair_benchmark.tsv", sep="\t", index=False)

    release_warnings: list[str] = []
    if smarts_release_summary.get("reaction_smarts_rules_all", 0) > 0 and smarts_release_summary.get("reaction_smarts_rules_core", 0) == 0:
        release_warnings.append(
            "Predictive SMARTS rules were produced, but the strict core SMARTS release is empty; "
            "use --require-core-rules for production builds whose downstream input is reaction_smarts_library.T1_core.tsv, "
            "or inspect T2_extended/T3_exploratory tiers explicitly."
        )

    summary = {
        "source_records": int(len(source_df)),
        "family_evidence_records": int(len(family_evidence)),
        "main_pair_records": int(len(main_df)),
        "raw_templates_or_anchors": int(len(raw_templates)),
        "deduplicated_templates_or_anchors": int(len(dedup)),
        "base_audit_records_before_data_driven_consensus": base_rule_count,
        "audit_records_all": int(len(rules)),
        "predictive_generalized_rules": int(((rules["template_scope"] == "generalized_template") & (_bool_mask(rules, "predictive_rule_use"))).sum()) if len(rules) else 0,
        "exact_anchor_rules": int((rules["template_scope"] == "exact_anchor").sum()) if len(rules) else 0,
        "anchor_edges": int(len(anchors)),
        "known_taxol_reactions": int(len(recall)),
        "known_taxol_recovered_by_external_generalized_rule": int(recall["recovered_by_external_generalized_rule"].sum()) if len(recall) else 0,
        "benchmark": benchmark_summary,
        "exact_pairs_promoted_to_predictive_rules": bool(allow_exact_pairs_as_predictive_rules),
        "exact_reaction_abstraction": exact_abstraction_summary,
        "anchor_generalization": anchor_generalization_summary,
        "data_driven_consensus": data_driven_consensus_summary,
        "transfer_family_consensus": transfer_family_consensus_summary,
        "participant_role_registry": participant_registry_summary,
        "rule_library_scope": "database_derived_general_one_step_enzymatic_transformation_smarts_library",
        "primary_downstream_network_input": "04_release/reaction_smarts_library.T1_core.tsv",
        "release_warnings": release_warnings,
        "core_release_empty_warning": bool(smarts_release_summary.get("reaction_smarts_rules_all", 0) > 0 and smarts_release_summary.get("reaction_smarts_rules_core", 0) == 0),
        **smarts_release_summary,
    }
    write_json(out / "build_summary.json", summary)
    if require_smarts_rules and smarts_release_summary.get("reaction_smarts_rules_all", 0) == 0:
        raise ValueError(
            "No validated predictive reaction SMARTS rules were produced in the all tier. "
            "Provide generalized reaction_smarts templates from MolLink/RetroRules/Rhea-derived sources before network construction."
        )
    if require_core_rules and smarts_release_summary.get("reaction_smarts_rules_core", 0) == 0:
        raise ValueError(
            "No strict-core predictive reaction SMARTS rules were produced. "
            "The primary downstream network input reaction_smarts_library.T1_core.tsv is empty; inspect build_summary.json and T2/T3 releases or relax scoring/evidence thresholds."
        )
    return summary


def run_build_from_raw(
    *,
    external_db_root: str | Path,
    output_root: str | Path,
    taxol_pathway: str | Path | None = None,
    cofactor_yaml: str | Path | None = None,
    derive_participant_registry: bool = True,
    participant_registry_min_frequency: int = 3,
    participant_registry_max_heavy_atoms: int = 80,
    scoring_yaml: str | Path | None = None,
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
    allow_exact_pairs_as_predictive_rules: bool = False,
    abstract_exact_reactions: bool | None = None,
    exact_abstraction_layers: list[str] | None = None,
    exact_abstraction_batch_size: int = 16,
    exact_abstraction_min_mapper_confidence: float = 0.50,
    exact_abstraction_max_records: int | None = None,
    generalize_exact_anchors: bool = False,
    anchor_generalization_mode: str = "strict",
    anchor_generalization_layers: list[str] | None = None,
    anchor_generalization_batch_size: int = 16,
    anchor_generalization_min_mapper_confidence: float = 0.50,
    anchor_generalization_max_records: int | None = None,
    anchor_generalization_max_reaction_chars: int = 2000,
    anchor_generalization_max_mapper_tokens: int = 500,
    anchor_generalization_mapper_timeout_seconds: int = 120,
    max_decoys: int = 500,
    skip_benchmark: bool = False,
    require_smarts_rules: bool = False,
    require_core_rules: bool = False,
    transfer_family_consensus: bool = True,
    transfer_family_mcs_timeout: int = 10,
    data_driven_consensus: bool = True,
    consensus_min_evidence_rows: int = 3,
    consensus_min_source_database_count: int = 1,
    consensus_min_template_count: int = 2,
) -> dict:
    out = ensure_dir(output_root)
    manifest_dir = ensure_dir(out / "00_raw_manifest")
    resolved_include_taxol_anchors = bool(taxol_pathway) if include_taxol_anchors is None else bool(include_taxol_anchors)
    raw_manifest, raw_discovery = build_raw_source_manifest(
        external_db_root=external_db_root,
        output_dir=manifest_dir,
        taxol_pathway=taxol_pathway,
        include_taxol_anchors=resolved_include_taxol_anchors,
        include_uspto=include_uspto,
        include_annotation_only=include_annotation_only,
        use_processed_bionavi=use_processed_bionavi,
        max_retrorules=max_retrorules,
        max_rhea=max_rhea,
        max_bionavi_per_file=max_bionavi_per_file,
        max_uspto=max_uspto,
        max_kegg=max_kegg,
        max_metanetx=max_metanetx,
    )
    if abstract_exact_reactions is None:
        abstract_exact_reactions = bool(resolved_include_taxol_anchors) and not generalize_exact_anchors
    summary = run_build_all(
        manifest=raw_manifest,
        output_root=out,
        taxol_pathway=None,
        cofactor_yaml=cofactor_yaml,
        derive_participant_registry=derive_participant_registry,
        participant_registry_min_frequency=participant_registry_min_frequency,
        participant_registry_max_heavy_atoms=participant_registry_max_heavy_atoms,
        scoring_yaml=scoring_yaml,
        allow_exact_pairs_as_predictive_rules=allow_exact_pairs_as_predictive_rules,
        abstract_exact_reactions=abstract_exact_reactions,
        exact_abstraction_layers=exact_abstraction_layers,
        exact_abstraction_batch_size=exact_abstraction_batch_size,
        exact_abstraction_min_mapper_confidence=exact_abstraction_min_mapper_confidence,
        exact_abstraction_max_records=exact_abstraction_max_records,
        generalize_exact_anchors=generalize_exact_anchors,
        anchor_generalization_mode=anchor_generalization_mode,
        anchor_generalization_layers=anchor_generalization_layers,
        anchor_generalization_batch_size=anchor_generalization_batch_size,
        anchor_generalization_min_mapper_confidence=anchor_generalization_min_mapper_confidence,
        anchor_generalization_max_records=anchor_generalization_max_records,
        anchor_generalization_max_reaction_chars=anchor_generalization_max_reaction_chars,
        anchor_generalization_max_mapper_tokens=anchor_generalization_max_mapper_tokens,
        anchor_generalization_mapper_timeout_seconds=anchor_generalization_mapper_timeout_seconds,
        max_decoys=max_decoys,
        skip_benchmark=skip_benchmark,
        require_smarts_rules=require_smarts_rules,
        require_core_rules=require_core_rules,
        transfer_family_consensus=transfer_family_consensus,
        transfer_family_mcs_timeout=transfer_family_mcs_timeout,
        data_driven_consensus=data_driven_consensus,
        consensus_min_evidence_rows=consensus_min_evidence_rows,
        consensus_min_source_database_count=consensus_min_source_database_count,
        consensus_min_template_count=consensus_min_template_count,
    )
    combined = {
        **summary,
        "raw_build": {
            "requested_external_db_root": raw_discovery.get("requested_external_db_root", str(Path(external_db_root).resolve())),
            "resolved_external_db_root": raw_discovery.get("resolved_external_db_root", raw_discovery.get("external_db_root", "")),
            "external_db_root": raw_discovery.get("resolved_external_db_root", raw_discovery.get("external_db_root", "")),
            "root_resolution_mode": raw_discovery.get("root_resolution_mode", ""),
            "manifest": str(raw_manifest.resolve()),
            "discovery_summary": raw_discovery,
            "legacy_network_ready_inputs": "not_used",
            "include_taxol_anchors_resolved": bool(resolved_include_taxol_anchors),
        },
    }
    write_json(out / "build_from_raw_summary.json", combined)
    return combined
