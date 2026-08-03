from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .io import ensure_dir, read_table, write_json, write_table


SOURCE_REACTIONS = Path("01_sources/source_reactions.normalized.tsv")
MAIN_PAIRS = Path("01_sources/source_reactions.main_pair.tsv")
RELEASE_FILES = {
    "T1": Path("04_release/reaction_smarts_library.T1_only.tsv"),
    "T2": Path("04_release/reaction_smarts_library.T2_only.tsv"),
    "T3": Path("04_release/reaction_smarts_library.T3_only.tsv"),
}


def _split_sources(value: object) -> list[str]:
    sources = {
        token.strip()
        for token in str(value or "").split(";")
        if token.strip() and token.strip().lower() != "nan"
    }
    return sorted(sources)


def _chunk_counts(
    path: Path,
    column: str,
    *,
    chunksize: int = 100_000,
) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    total = 0
    for frame in pd.read_csv(
        path,
        sep="\t",
        usecols=[column],
        dtype=str,
        chunksize=chunksize,
    ):
        values = frame[column].fillna("").replace("", "unattributed")
        counts.update({str(key): int(value) for key, value in values.value_counts().items()})
        total += len(frame)
    return counts, total


def _release_source_support(
    path: Path,
    *,
    chunksize: int = 50_000,
) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    total = 0
    for frame in pd.read_csv(
        path,
        sep="\t",
        usecols=["template_sources"],
        dtype=str,
        chunksize=chunksize,
    ):
        for value in frame["template_sources"].fillna(""):
            sources = _split_sources(value)
            if not sources:
                counts["unattributed"] += 1
            else:
                counts.update(sources)
        total += len(frame)
    return counts, total


def _source_table(build_root: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    source_counts, source_total = _chunk_counts(
        build_root / SOURCE_REACTIONS, "source_database"
    )
    pair_counts, pair_total = _chunk_counts(
        build_root / MAIN_PAIRS, "source_database"
    )
    tier_counts: dict[str, Counter[str]] = {}
    tier_totals: dict[str, int] = {}
    for tier, relative_path in RELEASE_FILES.items():
        tier_counts[tier], tier_totals[tier] = _release_source_support(
            build_root / relative_path
        )

    sources = sorted(
        set(source_counts)
        | set(pair_counts)
        | {source for counts in tier_counts.values() for source in counts}
    )
    rows = []
    for source in sources:
        rows.append(
            {
                "source_database": source,
                "normalized_source_reaction_records": source_counts[source],
                "main_pair_records": pair_counts[source],
                "T1_rule_rows_supported": tier_counts["T1"][source],
                "T2_rule_rows_supported": tier_counts["T2"][source],
                "T3_rule_rows_supported": tier_counts["T3"][source],
                "release_support_counting_scope": (
                    "overlapping database support; a multi-source rule contributes "
                    "one count to every supporting database"
                ),
            }
        )
    totals = {
        "normalized_source_reaction_records": source_total,
        "main_pair_records": pair_total,
        **{f"{tier}_exclusive_rule_rows": tier_totals[tier] for tier in RELEASE_FILES},
    }
    return pd.DataFrame(rows), totals


def _stage_table(build_summary: dict[str, Any]) -> pd.DataFrame:
    exact = build_summary.get("exact_reaction_abstraction", {})
    anchor = build_summary.get("anchor_generalization", {})
    strict_anchor = anchor.get("mode_summaries", {}).get("strict", {})
    permissive_anchor = anchor.get("mode_summaries", {}).get("permissive", {})
    consensus = build_summary.get("data_driven_consensus", {})
    rows = [
        (
            "normalized source reactions",
            build_summary.get("source_records", 0),
            "Database records after parser-level normalization.",
        ),
        (
            "main-substrate/main-product projections",
            build_summary.get("main_pair_records", 0),
            "Direction-normalized main molecular pairs after participant-role handling.",
        ),
        (
            "raw templates or exact anchors",
            build_summary.get("raw_templates_or_anchors", 0),
            "All template and exact-pair records entering rule QC.",
        ),
        (
            "deduplicated templates or anchors",
            build_summary.get("deduplicated_templates_or_anchors", 0),
            "Template/anchor records after build-level deduplication.",
        ),
        (
            "exact reactions considered for abstraction",
            exact.get("attempted_records", 0),
            "Exact reactions submitted to atom mapping and reaction-center abstraction.",
        ),
        (
            "successful exact-reaction abstractions",
            exact.get("successful_templates", 0),
            "Exact reactions yielding an extracted template before release-level selection.",
        ),
        (
            "strict replay-valid anchor-derived SMARTS",
            strict_anchor.get("anchor_derived_released_smarts", 0),
            "Anchor-derived rules passing strict mapper-confidence and replay gates.",
        ),
        (
            "permissive annotated anchor-derived SMARTS",
            permissive_anchor.get("anchor_derived_released_smarts", 0),
            "Extractable anchor-derived rules retained with explicit QC annotations.",
        ),
        (
            "role-aware transfer-family candidate rows",
            build_summary.get("transfer_family_consensus", {}).get(
                "released_template_rows", 0
            ),
            "Main-pair projected transfer templates passing donor/acceptor and replay QC.",
        ),
        (
            "data-driven consensus rules",
            consensus.get("promoted_consensus_rules", 0),
            "Recurrent generalized SMARTS promoted from evidence-cluster consensus.",
        ),
        (
            "predictive generalized SMARTS rules",
            build_summary.get("predictive_generalized_rules", 0),
            "All released generalized SMARTS across exclusive evidence tiers.",
        ),
        (
            "T1 evidence tier",
            build_summary.get("reaction_smarts_library_T1_only", 0),
            "Exclusive T1 evidence-supported generalized rule rows.",
        ),
        (
            "T2 evidence tier",
            build_summary.get("reaction_smarts_library_T2_only", 0),
            "Exclusive T2 extended-evidence generalized rule rows.",
        ),
        (
            "T3 evidence tier",
            build_summary.get("reaction_smarts_library_T3_only", 0),
            "Exclusive T3 exploratory generalized rule rows.",
        ),
    ]
    return pd.DataFrame(
        rows,
        columns=["build_stage", "record_count", "operational_definition"],
    )


def _prepared_tier_table(
    summary_paths: Iterable[tuple[str, Path | None]],
) -> pd.DataFrame:
    rows = []
    for tier, path in summary_paths:
        if path is None or not path.exists():
            continue
        summary = json.loads(path.read_text(encoding="utf-8"))
        rows.append(
            {
                "tier": tier.upper(),
                "input_release_rows": int(summary.get("rows_seen", 0)),
                "grammar_qc_eligible_rows": int(summary.get("eligible_rules", 0)),
                "semantic_groups": int(summary.get("semantic_groups", 0)),
                "initial_executable_representatives": int(
                    summary.get("selected_grammar_rules", 0)
                ),
                "single_center_required": bool(
                    summary.get("require_single_center", False)
                ),
                "summary_path": str(path),
            }
        )
    return pd.DataFrame(rows).sort_values("tier") if rows else pd.DataFrame()


def _final_grammar_tables(
    grammar_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    grammar = read_table(grammar_path).fillna("")
    forbidden = [column for column in grammar.columns if "scaffold" in column.lower()]
    if forbidden:
        raise ValueError(
            "Final grammar contains forbidden scaffold fields: " + ", ".join(forbidden)
        )

    source_counts: Counter[str] = Counter()
    for value in grammar.get("template_sources", pd.Series(dtype=str)):
        sources = _split_sources(value)
        source_counts.update(sources or ["unattributed"])
    source_frame = pd.DataFrame(
        [
            {
                "source_database": source,
                "final_grammar_rule_rows_supported": count,
                "counting_scope": (
                    "overlapping source support; multi-source rules count once per source"
                ),
            }
            for source, count in source_counts.most_common()
        ]
    )
    reaction_frame = (
        grammar.get("reaction_type", pd.Series(dtype=str))
        .replace("", "unassigned")
        .value_counts()
        .rename_axis("reaction_type")
        .reset_index(name="final_grammar_rule_count")
    )
    preferred_columns = [
        "final_grammar_rule_id",
        "reaction_smarts",
        "reaction_smarts_hash",
        "semantic_group_id",
        "reaction_type",
        "biochemical_step_granularity",
        "normalized_direction",
        "template_sources",
        "evidence_layer_best",
        "grammar_provenance_scope",
        "open_grammar_scope",
        "expected_element_delta",
        "structural_element_delta",
        "grammar_selection_score",
        "g0_match_count",
        "g0_site_match_count",
    ]
    provenance = grammar[
        [column for column in preferred_columns if column in grammar.columns]
    ].copy()
    return provenance, source_frame, reaction_frame


def summarize_rule_library_provenance(
    build_root: Path,
    final_grammar_path: Path,
    output_dir: Path,
    *,
    prepared_t1_summary: Path | None = None,
    prepared_t2_summary: Path | None = None,
    prepared_t3_summary: Path | None = None,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    build_summary_path = build_root / "build_summary.json"
    build_summary = json.loads(build_summary_path.read_text(encoding="utf-8"))

    source_frame, release_totals = _source_table(build_root)
    stage_frame = _stage_table(build_summary)
    tier_frame = _prepared_tier_table(
        [
            ("T1", prepared_t1_summary),
            ("T2", prepared_t2_summary),
            ("T3", prepared_t3_summary),
        ]
    )
    final_provenance, final_sources, final_reactions = _final_grammar_tables(
        final_grammar_path
    )

    paths = {
        "build_stages": output_dir / "rule_library_build_stage_counts.tsv",
        "source_contributions": output_dir
        / "source_database_rule_contributions.tsv",
        "prepared_tiers": output_dir / "prepared_grammar_tier_composition.tsv",
        "final_rule_provenance": output_dir / "final_primary_rule_provenance.tsv",
        "final_source_support": output_dir
        / "final_primary_grammar_source_support.tsv",
        "final_reaction_types": output_dir
        / "final_primary_grammar_reaction_types.tsv",
        "summary": output_dir / "rule_library_provenance_summary.json",
    }
    write_table(stage_frame, paths["build_stages"])
    write_table(source_frame, paths["source_contributions"])
    write_table(tier_frame, paths["prepared_tiers"])
    write_table(final_provenance, paths["final_rule_provenance"])
    write_table(final_sources, paths["final_source_support"])
    write_table(final_reactions, paths["final_reaction_types"])

    summary = {
        "mode": "rule_library_provenance",
        "build_root": str(build_root),
        "build_summary": str(build_summary_path),
        "final_grammar": str(final_grammar_path),
        "release_totals": release_totals,
        "final_primary_rule_rows": int(len(final_provenance)),
        "final_primary_reaction_types": int(len(final_reactions)),
        "final_primary_supporting_databases": int(len(final_sources)),
        "source_support_counts_are_overlapping": True,
        "exclusive_release_tiers": ["T1", "T2", "T3"],
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
