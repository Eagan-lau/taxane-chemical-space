from __future__ import annotations

import json
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .chemistry import require_rdkit
from .generate import _fingerprint
from .io import ensure_dir, read_table, write_json, write_table


FUNCTIONAL_GROUP_SMARTS = {
    "free_hydroxyl": "[OX2H1;!$([OX2H1][CX3](=O))]",
    "ester": "[CX3](=[OX1])[OX2][#6]",
    "carboxylic_acid_or_carboxylate": "[CX3](=[OX1])[OX1-,OX2H1]",
    "ketone_or_aldehyde": "[CX3](=[OX1])([#6,#1])[#6,#1]",
    "ether": "[OD2]([#6])[#6]",
    "amide": "[NX3][CX3](=[OX1])",
    "epoxide": "[OX2r3]1[CX4r3][CX4r3]1",
    "alkene": "[CX3]=[CX3]",
    "aromatic_ring_atom": "[a]",
    "phosphate": "[PX4](=[OX1])([OX1-,OX2])([OX1-,OX2])[OX1-,OX2]",
}


def _interpretation_layer(generation: int) -> str:
    if generation == 0:
        return "known_taxane_seed_space"
    if generation in (1, 2):
        return "primary_near_seed_chemical_space"
    return "exploratory_frontier"


def _compile_functional_queries() -> dict[str, Any]:
    Chem, *_rest = require_rdkit()
    queries = {
        name: Chem.MolFromSmarts(smarts)
        for name, smarts in FUNCTIONAL_GROUP_SMARTS.items()
    }
    invalid = [name for name, query in queries.items() if query is None]
    if invalid:
        raise RuntimeError(f"Invalid functional-group SMARTS: {invalid}")
    return queries


def _functional_state(mol, queries: dict[str, Any]) -> dict[str, int]:
    return {
        name: len(mol.GetSubstructMatches(query, uniquify=True))
        for name, query in queries.items()
    }


def _functional_transition_label(
    source_state: dict[str, int], target_state: dict[str, int]
) -> str:
    tokens = []
    for name in FUNCTIONAL_GROUP_SMARTS:
        delta = target_state[name] - source_state[name]
        if delta:
            tokens.append(f"{name}:{delta:+d}")
    return ";".join(tokens) if tokens else "no_counted_functional_state_change"


def _generation_summary(
    nodes: pd.DataFrame,
    events: pd.DataFrame,
    applications: pd.DataFrame,
    rejections: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    max_generation = int(nodes["generation_first"].max())
    for generation in range(max_generation + 1):
        generation_nodes = nodes[nodes["generation_first"] == generation]
        generation_events = (
            events[events["generation"] == generation]
            if generation > 0
            else events.iloc[0:0]
        )
        generation_applications = (
            applications[applications["generation"] == generation]
            if generation > 0
            else applications.iloc[0:0]
        )
        generation_rejections = (
            rejections[rejections["generation"] == generation]
            if generation > 0
            else rejections.iloc[0:0]
        )
        rows.append(
            {
                "generation": generation,
                "unique_nodes_first_observed": int(len(generation_nodes)),
                "derivation_events": int(len(generation_events)),
                "unique_parent_nodes": int(
                    generation_events["source_space_id"].nunique()
                ),
                "unique_target_nodes": int(
                    generation_events["target_space_id"].nunique()
                ),
                "unique_rules_used": int(
                    generation_events["grammar_rule_id"].nunique()
                ),
                "matched_parent_rule_pairs": int(len(generation_applications)),
                "raw_product_tuples": int(
                    pd.to_numeric(
                        generation_applications.get(
                            "raw_product_tuple_count", pd.Series(dtype=float)
                        ),
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),
                "rejected_product_events": int(len(generation_rejections)),
                "known_G0_full_recovery_events": int(
                    pd.to_numeric(
                        generation_events.get(
                            "known_g0_full_match", pd.Series(dtype=float)
                        ),
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),
                "known_G0_connectivity_recovery_events": int(
                    pd.to_numeric(
                        generation_events.get(
                            "known_g0_connectivity_match", pd.Series(dtype=float)
                        ),
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),
                "immediate_reverse_cycle_events": int(
                    pd.to_numeric(
                        generation_events.get(
                            "immediate_reverse_cycle", pd.Series(dtype=float)
                        ),
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),
            }
        )
    result = pd.DataFrame(rows)
    result["cumulative_unique_nodes"] = result[
        "unique_nodes_first_observed"
    ].cumsum()
    result["new_node_growth_ratio_vs_previous_generation"] = (
        result["unique_nodes_first_observed"]
        / result["unique_nodes_first_observed"].shift(1).replace(0, np.nan)
    )
    result["unique_node_yield_per_raw_product"] = (
        result["unique_nodes_first_observed"]
        / result["raw_product_tuples"].replace(0, np.nan)
    )
    result["accepted_event_fraction_of_raw_products"] = (
        result["derivation_events"]
        / result["raw_product_tuples"].replace(0, np.nan)
    )
    result["known_connectivity_reconnection_fraction"] = (
        result["known_G0_connectivity_recovery_events"]
        / result["derivation_events"].replace(0, np.nan)
    )
    result["interpretation_layer"] = result["generation"].map(
        _interpretation_layer
    )
    return result


def _descriptor_summary(nodes: pd.DataFrame) -> pd.DataFrame:
    descriptor_columns = [
        "exact_mass",
        "heavy_atom_count",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "ring_count",
        "fraction_csp3",
        "formal_charge",
    ]
    rows = []
    for generation, group in nodes.groupby("generation_first", sort=True):
        for descriptor in descriptor_columns:
            values = pd.to_numeric(group[descriptor], errors="coerce").dropna()
            rows.append(
                {
                    "generation": int(generation),
                    "interpretation_layer": _interpretation_layer(int(generation)),
                    "descriptor": descriptor,
                    "n": int(len(values)),
                    "mean": float(values.mean()),
                    "standard_deviation": float(values.std(ddof=1)),
                    "median": float(values.median()),
                    "q1": float(values.quantile(0.25)),
                    "q3": float(values.quantile(0.75)),
                    "minimum": float(values.min()),
                    "maximum": float(values.max()),
                }
            )
    return pd.DataFrame(rows)


def _nearest_g0_similarity(
    nodes: pd.DataFrame,
    *,
    max_nodes_per_generation: int | None = None,
    random_seed: int = 1729,
) -> pd.DataFrame:
    Chem, *_rest = require_rdkit()
    from rdkit import DataStructs

    g0 = nodes[nodes["generation_first"] == 0].copy()
    g0_mols = [Chem.MolFromSmiles(smiles) for smiles in g0["smiles"]]
    g0_fps = [_fingerprint(mol) for mol in g0_mols]
    g0_ids = list(g0["space_id"])
    selected = []
    rng = np.random.default_rng(random_seed)
    for generation, group in nodes.groupby("generation_first", sort=True):
        population_size = len(group)
        if (
            max_nodes_per_generation is not None
            and len(group) > max_nodes_per_generation
        ):
            indices = rng.choice(
                group.index.to_numpy(),
                size=max_nodes_per_generation,
                replace=False,
            )
            group = group.loc[np.sort(indices)]
        group = group.copy()
        group["_generation_population_size"] = population_size
        group["_generation_sample_size"] = len(group)
        selected.append(group)
    sampled = pd.concat(selected, ignore_index=True)
    rows = []
    for index, record in enumerate(sampled.to_dict("records"), start=1):
        mol = Chem.MolFromSmiles(str(record["smiles"]))
        fp = _fingerprint(mol)
        similarities = DataStructs.BulkTanimotoSimilarity(fp, g0_fps)
        best_index = int(np.argmax(similarities))
        rows.append(
            {
                "space_id": record["space_id"],
                "generation": int(record["generation_first"]),
                "interpretation_layer": _interpretation_layer(
                    int(record["generation_first"])
                ),
                "nearest_G0_space_id": g0_ids[best_index],
                "nearest_G0_tanimoto": round(float(similarities[best_index]), 6),
                "generation_population_size": int(
                    record["_generation_population_size"]
                ),
                "generation_sample_size": int(record["_generation_sample_size"]),
                "sampling_fraction": round(
                    float(record["_generation_sample_size"])
                    / float(record["_generation_population_size"]),
                    8,
                ),
                "sampling_policy": (
                    "all_nodes"
                    if int(record["_generation_sample_size"])
                    == int(record["_generation_population_size"])
                    else "deterministic_stratified_random_sample_within_generation"
                ),
            }
        )
        if index % 10_000 == 0:
            print(
                f"[analyze] nearest-G0 similarities={index}/{len(sampled)}",
                flush=True,
            )
    return pd.DataFrame(rows)


def _chemical_space_projection(
    nodes: pd.DataFrame,
    *,
    max_nodes_per_generation: int = 20_000,
    random_seed: int = 1729,
) -> pd.DataFrame:
    Chem, *_rest = require_rdkit()
    from scipy import sparse
    from sklearn.decomposition import TruncatedSVD

    rng = np.random.default_rng(random_seed)
    samples = []
    for generation, group in nodes.groupby("generation_first", sort=True):
        population_size = len(group)
        if len(group) > max_nodes_per_generation:
            indices = rng.choice(
                group.index.to_numpy(),
                size=max_nodes_per_generation,
                replace=False,
            )
            group = group.loc[np.sort(indices)]
        group = group.copy()
        group["_generation_population_size"] = population_size
        group["_generation_sample_size"] = len(group)
        samples.append(group)
    sample = pd.concat(samples, ignore_index=True)
    row_indices: list[int] = []
    column_indices: list[int] = []
    for row_index, smiles in enumerate(sample["smiles"]):
        mol = Chem.MolFromSmiles(str(smiles))
        fp = _fingerprint(mol)
        on_bits = list(fp.GetOnBits())
        row_indices.extend([row_index] * len(on_bits))
        column_indices.extend(on_bits)
    matrix = sparse.csr_matrix(
        (
            np.ones(len(row_indices), dtype=np.float32),
            (row_indices, column_indices),
        ),
        shape=(len(sample), 2048),
    )
    model = TruncatedSVD(n_components=2, random_state=random_seed)
    coordinates = model.fit_transform(matrix)
    result = sample[
        [
            "space_id",
            "generation_first",
            "smiles",
            "full_inchikey",
            "connectivity_key",
            "formula",
        ]
    ].copy()
    result = result.rename(columns={"generation_first": "generation"})
    result["interpretation_layer"] = result["generation"].map(
        _interpretation_layer
    )
    result["generation_population_size"] = sample[
        "_generation_population_size"
    ].astype(int)
    result["generation_sample_size"] = sample[
        "_generation_sample_size"
    ].astype(int)
    result["sampling_fraction"] = (
        result["generation_sample_size"] / result["generation_population_size"]
    )
    result["sampling_policy"] = np.where(
        result["generation_sample_size"] == result["generation_population_size"],
        "all_nodes",
        "deterministic_stratified_random_sample_within_generation",
    )
    result["molecular_fp_axis_1"] = coordinates[:, 0]
    result["molecular_fp_axis_2"] = coordinates[:, 1]
    result["axis_1_explained_variance_ratio"] = float(
        model.explained_variance_ratio_[0]
    )
    result["axis_2_explained_variance_ratio"] = float(
        model.explained_variance_ratio_[1]
    )
    result["projection_method"] = "Morgan_radius2_2048bit_TruncatedSVD"
    return result


def _functional_states_and_transitions(
    nodes: pd.DataFrame, events: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    Chem, *_rest = require_rdkit()
    queries = _compile_functional_queries()
    state_columns = list(FUNCTIONAL_GROUP_SMARTS)
    state_matrix = np.empty(
        (len(nodes), len(state_columns)), dtype=np.int16
    )
    node_view = nodes[
        ["space_id", "generation_first", "smiles"]
    ].reset_index(drop=True)
    for row_index, record in enumerate(
        node_view.itertuples(index=False, name=None)
    ):
        _space_id, _generation, smiles = record
        mol = Chem.MolFromSmiles(str(smiles))
        state = _functional_state(mol, queries)
        state_matrix[row_index, :] = [
            state[name] for name in state_columns
        ]
        if (row_index + 1) % 100_000 == 0:
            print(
                "[analyze] functional states="
                f"{row_index + 1}/{len(node_view)}",
                flush=True,
            )
    states = pd.DataFrame(state_matrix, columns=state_columns)
    states.insert(
        0,
        "interpretation_layer",
        node_view["generation_first"].map(
            lambda value: _interpretation_layer(int(value))
        ),
    )
    states.insert(
        0,
        "generation",
        node_view["generation_first"].to_numpy(dtype=np.int16),
    )
    states.insert(0, "space_id", node_view["space_id"].astype(str).to_numpy())
    state_indices = pd.Series(
        np.arange(len(states), dtype=np.int64),
        index=states["space_id"].astype(str),
    )
    source_indices = (
        events["source_space_id"].astype(str).map(state_indices).to_numpy()
    )
    target_indices = (
        events["target_space_id"].astype(str).map(state_indices).to_numpy()
    )
    if np.isnan(source_indices).any() or np.isnan(target_indices).any():
        raise ValueError("Functional-state analysis found unknown event endpoints")
    source_indices = source_indices.astype(np.int64)
    target_indices = target_indices.astype(np.int64)
    delta_matrix = (
        state_matrix[target_indices] - state_matrix[source_indices]
    ).astype(np.int16)
    unique_patterns, transition_codes = np.unique(
        delta_matrix, axis=0, return_inverse=True
    )
    definition_rows = []
    label_by_code: dict[int, str] = {}
    for code, pattern in enumerate(unique_patterns):
        source_state = {name: 0 for name in FUNCTIONAL_GROUP_SMARTS}
        target_state = {
            name: int(value)
            for name, value in zip(FUNCTIONAL_GROUP_SMARTS, pattern)
        }
        label = _functional_transition_label(source_state, target_state)
        label_by_code[code] = label
        definition_rows.append(
            {
                "functional_transition_code": code,
                "functional_state_transition": label,
                **{
                    f"delta_{name}": int(value)
                    for name, value in zip(FUNCTIONAL_GROUP_SMARTS, pattern)
                },
            }
        )
    definitions = pd.DataFrame(definition_rows)
    event_view = events[
        [
            "event_id",
            "generation",
            "source_space_id",
            "target_space_id",
            "grammar_rule_id",
            "semantic_group_id",
            "reaction_type",
        ]
    ].copy()
    event_view["functional_transition_code"] = transition_codes
    summary = (
        event_view.groupby(
            ["generation", "functional_transition_code"], dropna=False
        )
        .agg(
            derivation_event_count=("event_id", "size"),
            unique_source_count=("source_space_id", "nunique"),
            unique_target_count=("target_space_id", "nunique"),
            unique_rule_count=("grammar_rule_id", "nunique"),
        )
        .reset_index()
        .sort_values(
            ["generation", "derivation_event_count"],
            ascending=[True, False],
        )
    )
    summary["functional_state_transition"] = summary[
        "functional_transition_code"
    ].map(label_by_code)
    summary["interpretation_layer"] = summary["generation"].map(
        lambda value: _interpretation_layer(int(value))
    )
    rule_summary = (
        event_view.groupby(
            [
                "generation",
                "grammar_rule_id",
                "semantic_group_id",
                "reaction_type",
                "functional_transition_code",
            ],
            dropna=False,
        )
        .agg(
            derivation_event_count=("event_id", "size"),
            unique_source_count=("source_space_id", "nunique"),
            unique_target_count=("target_space_id", "nunique"),
        )
        .reset_index()
        .sort_values(
            ["generation", "derivation_event_count"],
            ascending=[True, False],
        )
    )
    rule_summary["functional_state_transition"] = rule_summary[
        "functional_transition_code"
    ].map(label_by_code)
    rule_summary["interpretation_layer"] = rule_summary["generation"].map(
        lambda value: _interpretation_layer(int(value))
    )
    return states, definitions, rule_summary, summary


def _grammar_usage(events: pd.DataFrame) -> pd.DataFrame:
    return (
        events.groupby(
            [
                "grammar_rule_id",
                "smarts_rule_id",
                "semantic_group_id",
                "reaction_type",
                "evidence_layer",
                "expected_element_delta",
            ],
            dropna=False,
        )
        .agg(
            derivation_event_count=("event_id", "size"),
            generation_min=("generation", "min"),
            generation_max=("generation", "max"),
            unique_source_count=("source_space_id", "nunique"),
            unique_target_count=("target_space_id", "nunique"),
            new_target_count=("target_is_new", "sum"),
            median_source_product_tanimoto=("source_product_tanimoto", "median"),
            median_source_atom_retention=("source_atom_retention", "median"),
        )
        .reset_index()
        .sort_values("derivation_event_count", ascending=False)
    )


def _gini(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[values >= 0]
    if not len(values) or float(values.sum()) == 0:
        return 0.0
    ordered = np.sort(values)
    cumulative = np.cumsum(ordered)
    n = len(ordered)
    return float(
        (n + 1 - 2 * float(cumulative.sum()) / float(cumulative[-1])) / n
    )


def _grammar_use_concentration(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    generation_sets = [
        ("all_generations", events),
        *[
            (f"G{int(generation)}", group)
            for generation, group in events.groupby("generation", sort=True)
        ],
    ]
    for generation_label, group in generation_sets:
        counts = group["grammar_rule_id"].value_counts().to_numpy(dtype=float)
        semantic_counts = (
            group["semantic_group_id"].value_counts().to_numpy(dtype=float)
        )
        total = float(counts.sum())
        probabilities = counts / total if total else np.asarray([], dtype=float)
        entropy = float(
            -(probabilities * np.log(probabilities)).sum()
        ) if len(probabilities) else 0.0
        normalized_entropy = (
            entropy / math.log(len(probabilities))
            if len(probabilities) > 1
            else 0.0
        )
        ordered = np.sort(counts)[::-1]
        rows.append(
            {
                "generation_scope": generation_label,
                "interpretation_layer": (
                    "all_layers"
                    if generation_label == "all_generations"
                    else _interpretation_layer(int(generation_label[1:]))
                ),
                "derivation_events": int(total),
                "rules_used": int(len(counts)),
                "semantic_groups_used": int(len(semantic_counts)),
                "shannon_entropy": entropy,
                "normalized_shannon_entropy": normalized_entropy,
                "effective_rule_number_exp_shannon": float(math.exp(entropy)),
                "herfindahl_hirschman_index": float(
                    np.square(probabilities).sum()
                )
                if len(probabilities)
                else 0.0,
                "gini_coefficient": _gini(counts),
                "top1_rule_event_fraction": float(ordered[:1].sum() / total)
                if total
                else 0.0,
                "top5_rule_event_fraction": float(ordered[:5].sum() / total)
                if total
                else 0.0,
                "top10_rule_event_fraction": float(ordered[:10].sum() / total)
                if total
                else 0.0,
            }
        )
    return pd.DataFrame(rows)


def _convergence_and_paths(
    nodes: pd.DataFrame, events: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_event_aggregation = (
        events.groupby("target_space_id", dropna=False)
        .agg(
            all_incoming_derivation_events=("event_id", "size"),
            all_unique_parent_count=("source_space_id", "nunique"),
            all_unique_rule_count=("grammar_rule_id", "nunique"),
            all_unique_semantic_group_count=(
                "semantic_group_id",
                "nunique",
            ),
            all_earliest_event_generation=("generation", "min"),
        )
        .reset_index()
    )
    first_observation_events = events.loc[
        events["generation"] == events["target_generation_first"]
    ]
    event_aggregation = (
        first_observation_events.groupby("target_space_id", dropna=False)
        .agg(
            incoming_derivation_events=("event_id", "size"),
            unique_parent_count=("source_space_id", "nunique"),
            unique_rule_count=("grammar_rule_id", "nunique"),
            unique_semantic_group_count=("semantic_group_id", "nunique"),
            earliest_event_generation=("generation", "min"),
        )
        .reset_index()
    )
    node_ids = nodes["space_id"].astype(str).to_numpy()
    node_index = pd.Series(np.arange(len(nodes), dtype=np.int64), index=node_ids)
    g0_mask = nodes["generation_first"].to_numpy(dtype=int) == 0

    raw_rule_event_path_counts = np.zeros(len(nodes), dtype=np.float64)
    semantic_edge_path_counts = np.zeros(len(nodes), dtype=np.float64)
    structural_path_counts = np.zeros(len(nodes), dtype=np.float64)
    for path_counts in (
        raw_rule_event_path_counts,
        semantic_edge_path_counts,
        structural_path_counts,
    ):
        path_counts[g0_mask] = 1.0

    for generation in sorted(
        int(value)
        for value in nodes["generation_first"].unique()
        if int(value) > 0
    ):
        generation_events = events.loc[
            (events["generation"] == generation)
            & (events["target_generation_first"] == generation),
            [
                "source_space_id",
                "target_space_id",
                "semantic_group_id",
            ],
        ]
        pair_weights = (
            generation_events.groupby(
                ["target_space_id", "source_space_id"],
                sort=False,
                observed=True,
            )
            .agg(
                raw_rule_event_multiplicity=(
                    "semantic_group_id",
                    "size",
                ),
                semantic_edge_multiplicity=(
                    "semantic_group_id",
                    "nunique",
                ),
            )
            .reset_index()
        )
        source_indices = (
            pair_weights["source_space_id"]
            .astype(str)
            .map(node_index)
            .to_numpy(dtype=np.int64)
        )
        target_ids = pair_weights["target_space_id"].astype(str).to_numpy()
        weighted = pd.DataFrame(
            {
                "target_space_id": target_ids,
                "structural": structural_path_counts[source_indices],
                "semantic": (
                    semantic_edge_path_counts[source_indices]
                    * pair_weights[
                        "semantic_edge_multiplicity"
                    ].to_numpy(dtype=np.float64)
                ),
                "raw": (
                    raw_rule_event_path_counts[source_indices]
                    * pair_weights[
                        "raw_rule_event_multiplicity"
                    ].to_numpy(dtype=np.float64)
                ),
            }
        )
        target_counts = weighted.groupby(
            "target_space_id", sort=False
        )[["structural", "semantic", "raw"]].sum()
        target_indices = (
            target_counts.index.to_series()
            .map(node_index)
            .to_numpy(dtype=np.int64)
        )
        structural_path_counts[target_indices] = target_counts[
            "structural"
        ].to_numpy(dtype=float)
        semantic_edge_path_counts[target_indices] = target_counts[
            "semantic"
        ].to_numpy(dtype=float)
        raw_rule_event_path_counts[target_indices] = target_counts[
            "raw"
        ].to_numpy(dtype=float)
        print(
            "[analyze] path layers="
            f"G{generation};pairs={len(pair_weights)}",
            flush=True,
        )
    path_frame = pd.DataFrame(
        {
            "space_id": node_ids,
            "structural_path_count": structural_path_counts,
            "log10_structural_path_count": np.where(
                structural_path_counts > 0,
                np.log10(structural_path_counts),
                np.nan,
            ),
            "semantic_edge_path_count": semantic_edge_path_counts,
            "log10_semantic_edge_path_count": np.where(
                semantic_edge_path_counts > 0,
                np.log10(semantic_edge_path_counts),
                np.nan,
            ),
            "raw_rule_event_path_count": raw_rule_event_path_counts,
            "log10_raw_rule_event_path_count": np.where(
                raw_rule_event_path_counts > 0,
                np.log10(raw_rule_event_path_counts),
                np.nan,
            ),
        }
    )
    convergence = nodes[
        ["space_id", "generation_first", "smiles", "formula", "full_inchikey"]
    ].merge(event_aggregation, left_on="space_id", right_on="target_space_id", how="left")
    convergence = convergence.drop(
        columns=["target_space_id"], errors="ignore"
    )
    convergence = convergence.merge(
        all_event_aggregation,
        left_on="space_id",
        right_on="target_space_id",
        how="left",
    )
    convergence = convergence.merge(
        path_frame, on="space_id", how="left", suffixes=("", "_path")
    )
    convergence = convergence.drop(columns=["target_space_id"], errors="ignore")
    convergence = convergence.rename(columns={"generation_first": "generation"})
    for column in (
        "incoming_derivation_events",
        "unique_parent_count",
        "unique_rule_count",
        "unique_semantic_group_count",
        "all_incoming_derivation_events",
        "all_unique_parent_count",
        "all_unique_rule_count",
        "all_unique_semantic_group_count",
    ):
        convergence[column] = convergence[column].fillna(0).astype(int)
    convergence["later_rediscovery_event_count"] = (
        convergence["all_incoming_derivation_events"]
        - convergence["incoming_derivation_events"]
    )
    convergence["multiple_parent_routes"] = (
        convergence["unique_parent_count"] >= 2
    )
    convergence["multiple_semantic_route_support"] = (
        convergence["unique_semantic_group_count"] >= 2
    )
    convergence["multiple_rule_support"] = (
        convergence["unique_rule_count"] >= 2
    )
    convergence["is_convergent"] = convergence["multiple_parent_routes"]
    convergence["convergence_definition"] = (
        "at_least_two_distinct_parent_structures_in_the_targets_"
        "first_observed_generation"
    )
    convergence["interpretation_layer"] = convergence["generation"].map(
        lambda value: _interpretation_layer(int(value))
    )

    generation_by_id = pd.Series(
        nodes["generation_first"].to_numpy(dtype=int),
        index=node_ids,
    )
    edge_frame = events[
        ["source_space_id", "target_space_id"]
    ].drop_duplicates()
    edge_frame["source_generation"] = (
        edge_frame["source_space_id"].astype(str).map(generation_by_id)
    )
    edge_frame["target_generation"] = (
        edge_frame["target_space_id"].astype(str).map(generation_by_id)
    )
    forward_edges = edge_frame[
        edge_frame["target_generation"]
        == edge_frame["source_generation"] + 1
    ].copy()
    recovery_edges = edge_frame[
        (edge_frame["source_generation"] > 0)
        & (edge_frame["target_generation"] == 0)
    ].copy()
    ordered_g0_ids = sorted(
        nodes.loc[nodes["generation_first"] == 0, "space_id"].astype(str)
    )
    g0_bit = {
        space_id: 1 << index for index, space_id in enumerate(ordered_g0_ids)
    }
    ancestors = [g0_bit.get(space_id, 0) for space_id in node_ids]
    descendants = [g0_bit.get(space_id, 0) for space_id in node_ids]
    forward_edges = forward_edges.sort_values(
        ["source_generation", "source_space_id", "target_space_id"]
    )
    forward_source_indices = (
        forward_edges["source_space_id"].astype(str).map(node_index).to_numpy(dtype=np.int64)
    )
    forward_target_indices = (
        forward_edges["target_space_id"].astype(str).map(node_index).to_numpy(dtype=np.int64)
    )
    for source_index, target_index in zip(
        forward_source_indices, forward_target_indices
    ):
        ancestors[target_index] |= ancestors[source_index]
    recovery_source_indices = (
        recovery_edges["source_space_id"].astype(str).map(node_index).to_numpy(dtype=np.int64)
    )
    recovery_target_indices = (
        recovery_edges["target_space_id"].astype(str).map(node_index).to_numpy(dtype=np.int64)
    )
    for source_index, target_index in zip(
        recovery_source_indices, recovery_target_indices
    ):
        descendants[source_index] |= descendants[target_index]
    for source_index, target_index in zip(
        reversed(forward_source_indices), reversed(forward_target_indices)
    ):
        descendants[source_index] |= descendants[target_index]

    def bit_ids(mask: int, limit: int = 20) -> tuple[list[str], str, str]:
        indices: list[int] = []
        remaining = mask
        while remaining:
            least_significant_bit = remaining & -remaining
            indices.append(least_significant_bit.bit_length() - 1)
            remaining ^= least_significant_bit
        ids = [ordered_g0_ids[index] for index in indices]
        sample = ";".join(ids[:limit])
        digest = (
            __import__("hashlib")
            .sha256(";".join(ids).encode("utf-8"))
            .hexdigest()[:20]
        )
        return ids, sample, digest

    bridge_rows: list[tuple[Any, ...]] = []
    pair_bridge_counts: Counter[tuple[int, str, str]] = Counter()
    node_generations = nodes["generation_first"].to_numpy(dtype=np.int16)
    generated_indices = np.flatnonzero(
        node_generations > 0
    )
    for progress_index, node_position in enumerate(generated_indices, start=1):
        node_id = node_ids[node_position]
        ancestor_mask = ancestors[node_position]
        descendant_mask = descendants[node_position]
        ancestor_count = int(ancestor_mask.bit_count())
        descendant_count = int(descendant_mask.bit_count())
        shared_count = int(
            (ancestor_mask & descendant_mask).bit_count()
        )
        cross_pair_count = ancestor_count * descendant_count - shared_count
        generation_value = int(node_generations[node_position])
        if cross_pair_count:
            ancestor_ids, ancestor_sample, ancestor_hash = bit_ids(
                ancestor_mask
            )
            descendant_ids, descendant_sample, descendant_hash = bit_ids(
                descendant_mask
            )
            for ancestor_id in ancestor_ids:
                for descendant_id in descendant_ids:
                    if ancestor_id != descendant_id:
                        pair_bridge_counts[
                            (generation_value, ancestor_id, descendant_id)
                        ] += 1
        else:
            ancestor_sample = ""
            descendant_sample = ""
            ancestor_hash = ""
            descendant_hash = ""
        bridge_rows.append(
            (
                node_id,
                generation_value,
                _interpretation_layer(generation_value),
                ancestor_count,
                descendant_count,
                cross_pair_count,
                ancestor_sample,
                descendant_sample,
                ancestor_hash,
                descendant_hash,
                bool(cross_pair_count),
                "generation_increasing_novel_edges_plus_direct_G0_recovery_edges",
            )
        )
        if progress_index % 250_000 == 0:
            print(
                "[analyze] bridge lineage="
                f"{progress_index}/{len(generated_indices)}",
                flush=True,
            )
    bridges = pd.DataFrame(
        bridge_rows,
        columns=[
            "space_id",
            "generation",
            "interpretation_layer",
            "known_G0_ancestor_count",
            "known_G0_descendant_count",
            "distinct_G0_pair_bridge_count",
            "known_G0_ancestor_ids_sample",
            "known_G0_descendant_ids_sample",
            "known_G0_ancestor_set_hash",
            "known_G0_descendant_set_hash",
            "latent_bridge_candidate",
            "bridge_graph_policy",
        ],
    )
    pair_rows = [
        {
            "generation": generation,
            "interpretation_layer": _interpretation_layer(generation),
            "known_G0_source_space_id": source,
            "known_G0_target_space_id": target,
            "supporting_latent_intermediate_count": count,
            "pair_direction": "source_to_target",
            "bridge_graph_policy": (
                "generation_increasing_novel_edges_plus_direct_G0_recovery_edges"
            ),
        }
        for (generation, source, target), count in pair_bridge_counts.most_common()
    ]
    pair_columns = [
        "generation",
        "interpretation_layer",
        "known_G0_source_space_id",
        "known_G0_target_space_id",
        "supporting_latent_intermediate_count",
        "pair_direction",
        "bridge_graph_policy",
    ]
    return convergence, bridges, pd.DataFrame(pair_rows, columns=pair_columns)


def _edit_landscape(events: pd.DataFrame) -> pd.DataFrame:
    return (
        events.groupby(
            [
                "generation",
                "observed_element_delta",
                "expected_element_delta",
                "observed_changed_source_atoms",
            ],
            dropna=False,
        )
        .agg(
            derivation_event_count=("event_id", "size"),
            unique_rule_count=("grammar_rule_id", "nunique"),
            unique_source_count=("source_space_id", "nunique"),
            unique_target_count=("target_space_id", "nunique"),
            median_source_atom_retention=("source_atom_retention", "median"),
            median_source_product_tanimoto=("source_product_tanimoto", "median"),
        )
        .reset_index()
        .sort_values(
            ["generation", "derivation_event_count"],
            ascending=[True, False],
        )
    )


def analyze_chemical_space(
    nodes_path: Path,
    events_path: Path,
    application_audit_path: Path,
    rejections_path: Path,
    output_dir: Path,
    *,
    projection_max_nodes_per_generation: int = 20_000,
    similarity_max_nodes_per_generation: int | None = 50_000,
    random_seed: int = 1729,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    nodes = read_table(nodes_path).fillna("")
    events = read_table(events_path).fillna("")
    applications = read_table(application_audit_path).fillna("")
    rejections = read_table(rejections_path).fillna("")
    numeric_node_columns = [
        "generation_first",
        "exact_mass",
        "heavy_atom_count",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "ring_count",
        "fraction_csp3",
        "formal_charge",
    ]
    for column in numeric_node_columns:
        nodes[column] = pd.to_numeric(nodes[column], errors="coerce")
    numeric_event_columns = [
        "event_id",
        "generation",
        "target_is_new",
        "target_generation_first",
        "source_atom_retention",
        "observed_changed_source_atoms",
        "source_product_tanimoto",
        "known_g0_full_match",
        "known_g0_connectivity_match",
        "immediate_reverse_cycle",
    ]
    for column in numeric_event_columns:
        events[column] = pd.to_numeric(events[column], errors="coerce").fillna(0)
    if not applications.empty:
        applications["generation"] = pd.to_numeric(
            applications["generation"], errors="coerce"
        ).fillna(0)
    if not rejections.empty:
        rejections["generation"] = pd.to_numeric(
            rejections["generation"], errors="coerce"
        ).fillna(0)

    generation_summary = _generation_summary(
        nodes, events, applications, rejections
    )
    descriptor_summary = _descriptor_summary(nodes)
    nearest_similarity = _nearest_g0_similarity(
        nodes,
        max_nodes_per_generation=similarity_max_nodes_per_generation,
        random_seed=random_seed,
    )
    projection = _chemical_space_projection(
        nodes,
        max_nodes_per_generation=projection_max_nodes_per_generation,
        random_seed=random_seed,
    )
    (
        states,
        transition_definitions,
        transition_rule_summary,
        transition_summary,
    ) = _functional_states_and_transitions(nodes, events)
    grammar_usage = _grammar_usage(events)
    grammar_concentration = _grammar_use_concentration(events)
    convergence, bridges, bridge_pairs = _convergence_and_paths(nodes, events)
    edit_landscape = _edit_landscape(events)
    rejection_summary = (
        rejections.groupby(["generation", "rejection_reason"], dropna=False)
        .size()
        .reset_index(name="rejected_product_count")
        .sort_values(
            ["generation", "rejected_product_count"], ascending=[True, False]
        )
        if not rejections.empty
        else pd.DataFrame(
            columns=["generation", "rejection_reason", "rejected_product_count"]
        )
    )
    paths = {
        "generation_summary": output_dir / "generation_expansion_summary.tsv",
        "descriptor_summary": output_dir / "physicochemical_descriptor_summary.tsv",
        "nearest_similarity": output_dir / "nearest_G0_similarity.tsv",
        "projection": output_dir / "chemical_space_fingerprint_projection.tsv",
        "functional_states": output_dir / "functional_state_counts.tsv",
        "functional_transition_definitions": output_dir
        / "functional_state_transition_definitions.tsv",
        "functional_transition_rule_summary": output_dir
        / "functional_state_transition_rule_summary.tsv",
        "functional_transition_summary": output_dir
        / "functional_state_transition_summary.tsv",
        "grammar_usage": output_dir / "reaction_grammar_usage.tsv",
        "grammar_concentration": output_dir
        / "reaction_grammar_use_concentration.tsv",
        "convergence": output_dir / "convergence_and_route_multiplicity.tsv",
        "bridges": output_dir / "latent_bridge_candidates.tsv",
        "bridge_pairs": output_dir / "known_G0_pair_bridge_summary.tsv",
        "edit_landscape": output_dir / "reaction_edit_landscape.tsv",
        "rejection_summary": output_dir / "generation_rejection_summary.tsv",
        "summary": output_dir / "chemical_space_analysis_summary.json",
    }
    table_map = {
        "generation_summary": generation_summary,
        "descriptor_summary": descriptor_summary,
        "nearest_similarity": nearest_similarity,
        "projection": projection,
        "functional_states": states,
        "functional_transition_definitions": transition_definitions,
        "functional_transition_rule_summary": transition_rule_summary,
        "functional_transition_summary": transition_summary,
        "grammar_usage": grammar_usage,
        "grammar_concentration": grammar_concentration,
        "convergence": convergence,
        "bridges": bridges,
        "bridge_pairs": bridge_pairs,
        "edit_landscape": edit_landscape,
        "rejection_summary": rejection_summary,
    }
    for key, frame in table_map.items():
        write_table(frame, paths[key])
    summary = {
        "mode": "taxane_reaction_grammar_space_analysis",
        "input_nodes": int(len(nodes)),
        "input_derivation_events": int(len(events)),
        "generations": sorted(int(value) for value in nodes["generation_first"].unique()),
        "rules_used": int(events["grammar_rule_id"].nunique()),
        "semantic_groups_used": int(events["semantic_group_id"].nunique()),
        "convergent_generated_nodes": int(
            convergence.loc[
                convergence["generation"] > 0, "is_convergent"
            ].astype(bool).sum()
        ),
        "convergence_definition": (
            "at_least_two_distinct_parent_structures_in_the_targets_"
            "first_observed_generation"
        ),
        "path_count_layers": {
            "primary": "distinct_source_target_structural_edges",
            "semantic_audit": (
                "distinct_source_target_semantic_group_edges"
            ),
            "raw_audit": "all_parent_rule_product_events",
        },
        "path_event_inclusion": (
            "all_events_recorded_in_the_targets_first_observed_generation_"
            "regardless_of_insertion_flag"
        ),
        "latent_bridge_candidates": int(
            bridges["latent_bridge_candidate"].astype(bool).sum()
        ),
        "projection_sample_size": int(len(projection)),
        "nearest_similarity_rows": int(len(nearest_similarity)),
        "projection_sampling_policy": (
            "deterministic_stratified_random_sample_within_generation"
        ),
        "projection_max_nodes_per_generation": int(
            projection_max_nodes_per_generation
        ),
        "nearest_similarity_sampling_policy": (
            "all_if_below_cap_else_deterministic_stratified_random_sample"
        ),
        "nearest_similarity_max_nodes_per_generation": (
            None
            if similarity_max_nodes_per_generation is None
            else int(similarity_max_nodes_per_generation)
        ),
        "generation_interpretation": {
            "G0": "known_taxane_seed_space",
            "G1_G2": "primary_near_seed_chemical_space",
            "G3": "exploratory_frontier_not_a_metabolite_catalogue",
        },
        "random_seed": random_seed,
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
