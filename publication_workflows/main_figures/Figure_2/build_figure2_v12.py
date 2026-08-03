#!/usr/bin/env python3
"""Render integrated Figure 2 V12 from frozen taxane chemical-space results.

This presentation-only workflow does not enumerate structures, rebuild rules,
or recalculate scientific metrics. It changes only the visual encoding of the
released G0-G3 chemical space: G0-G2 are node resolved and G3 remains a
complete component-normalized density layer with density-derived stippling.
All visual layers share one monotonic density-equalizing coordinate transform.
The G2 layer then undergoes bounded deterministic display relaxation to reduce
local overplotting and occupy nearby low-density regions. G2 membership,
parent links, and all G0/G1 graph topology remain unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sqlite3
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Patch, Rectangle
from scipy.ndimage import (
    binary_closing,
    binary_dilation,
    binary_fill_holes,
    gaussian_filter,
    map_coordinates,
)

import build_figure2_figure3_v6 as base


DEFAULT_WORK = Path(".")

G1_DISPLAY_LABELS = {
    "hydroxylation_or_oxygenation": "hydroxylation product",
    "acetylation_or_deacetylation_like_acyl_transfer": (
        "O-acetyl-transfer product"
    ),
    "deacetylation_or_acetyl_ester_hydrolysis": "deacetylation product",
}

G1_MARKERS = {
    "G1_00000695": "a",
    "G1_00003186": "b",
    "G1_00005175": "c",
}

GENERATION_COLORS = {
    "G0": base.COLORS["red"],
    "G1": base.COLORS["blue"],
    "G2": "#8a949c",
    "G3": base.COLORS["orange"],
}

G0_G1_NODE_AREA = 1.55
G2_NODE_AREA = G0_G1_NODE_AREA * 0.25
G2_DISPLAY_COLOR = "#8a949c"
G3_STIPPLE_COUNT = 50000
LAYOUT_EQUALIZATION_STRENGTH = 0.42
LAYOUT_TARGET_MARGIN_FRACTION = 0.065
G2_RELAXATION_GRID_SIZE = 180
G2_RELAXATION_ITERATIONS = 14
G2_RELAXATION_SMOOTHING_SIGMA = 2.35
G2_RELAXATION_STEP_FRACTION = 0.017
G2_RELAXATION_MAX_DISPLACEMENT_FRACTION = 0.042
G2_RELAXATION_BOUNDARY_MARGIN_FRACTION = 0.025

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_inputs(work: Path) -> dict[str, object]:
    data = base.load_inputs(work)
    projected = pd.read_csv(data["paths"]["projected"], sep="\t")
    network = data["network"]
    nodes = network[network["record_type"].eq("node")].copy()
    nodes["generation"] = pd.to_numeric(nodes["generation"]).astype(int)
    g1_nodes = nodes[nodes["generation"].eq(1)][
        [
            "space_id",
            "layout_x",
            "layout_y",
            "layout_component_index",
        ]
    ].rename(columns={"space_id": "target_space_id"})
    g1_examples = projected.merge(
        g1_nodes,
        on="target_space_id",
        how="left",
        validate="one_to_one",
    )
    g1_examples["display_label"] = g1_examples["reaction_type"].map(
        G1_DISPLAY_LABELS
    )
    g1_examples["network_marker"] = g1_examples["target_space_id"].map(
        G1_MARKERS
    )
    g1_examples["display_order"] = g1_examples["network_marker"].map(
        {"a": 0, "b": 1, "c": 2}
    )
    g1_examples = g1_examples.sort_values("display_order")
    data["g1_examples"] = g1_examples
    data["g2_points"] = recover_g2_layout(data["network"], work)
    data["g3_stipple"] = density_stipple_points(data["density"])
    data.update(
        build_equalized_display_layout(
            network=data["network"],
            density=data["density"],
            g2_points=data["g2_points"],
            g3_stipple=data["g3_stipple"],
            representative=data["representative"],
            g1_examples=data["g1_examples"],
        )
    )
    balanced_g2, relaxation_summary, relaxation_iterations = (
        relax_g2_display_layout(
            data["display_g2_points"],
            data["display_network"],
            np.asarray(data["display_density"]["bounds"], dtype=float),
        )
    )
    data["display_g2_points"] = balanced_g2
    data["g2_relaxation_summary"] = relaxation_summary
    data["g2_relaxation_iterations"] = relaxation_iterations
    data["paths"]["frozen_database"] = (
        work
        / "inputs/G0_G3_primary_release/03_primary_G0_G3/"
        "taxane_reaction_grammar_space.sqlite"
    )
    return data


def density_stipple_points(
    density,
    *,
    point_count: int = G3_STIPPLE_COUNT,
    seed: int = 20260730,
) -> pd.DataFrame:
    """Sample reproducible display marks from the frozen G3 density field.

    The marks are a pointillist rendering of the density background and are
    explicitly not interpreted as individual molecular structures.
    """
    intensity = np.asarray(
        density["g3_component_normalized_intensity"],
        dtype=np.float64,
    )
    bounds = np.asarray(density["bounds"], dtype=float)
    positive = np.isfinite(intensity) & (intensity > 0)
    if not positive.any():
        raise RuntimeError("Frozen G3 density contains no positive pixels")

    # A shallow display gamma preserves density support while making sparse
    # internal regions visible after publication-scale downsampling.
    weights = np.zeros_like(intensity, dtype=np.float64)
    weights[positive] = np.power(intensity[positive], 0.18)
    flat_weights = weights.ravel()
    flat_weights /= flat_weights.sum()

    rng = np.random.default_rng(seed)
    flat_indices = rng.choice(
        flat_weights.size,
        size=point_count,
        replace=True,
        p=flat_weights,
    )
    rows, columns = np.unravel_index(flat_indices, intensity.shape)
    jitter_x = rng.random(point_count)
    jitter_y = rng.random(point_count)
    pixel_width = (bounds[1] - bounds[0]) / intensity.shape[1]
    pixel_height = (bounds[3] - bounds[2]) / intensity.shape[0]
    layout_x = bounds[0] + (columns + jitter_x) * pixel_width
    layout_y = bounds[2] + (rows + jitter_y) * pixel_height

    return pd.DataFrame(
        {
            "display_mark_id": [
                f"G3_DENSITY_STIPPLE_{index:05d}"
                for index in range(1, point_count + 1)
            ],
            "layout_x": layout_x,
            "layout_y": layout_y,
            "source_density_intensity": intensity[rows, columns],
            "visual_role": "G3_density_stipple",
            "represents_individual_structure": False,
            "sampling_seed": seed,
            "display_gamma": 0.18,
        }
    )


def build_axis_warp(
    values: np.ndarray,
    lower_bound: float,
    upper_bound: float,
    *,
    axis: str,
    strength: float = LAYOUT_EQUALIZATION_STRENGTH,
    margin_fraction: float = LAYOUT_TARGET_MARGIN_FRACTION,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """Build a softened empirical-CDF coordinate warp for one display axis."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError(f"No finite coordinates available for axis {axis}")

    quantiles = np.linspace(0.0, 1.0, 257)
    original = np.quantile(values, quantiles)
    span = upper_bound - lower_bound
    target_lower = lower_bound + margin_fraction * span
    target_upper = upper_bound - margin_fraction * span
    target = target_lower + quantiles * (target_upper - target_lower)
    display = (1.0 - strength) * original + strength * target

    original_knots = np.concatenate(
        ([lower_bound], original, [upper_bound])
    )
    display_knots = np.concatenate(
        ([lower_bound], display, [upper_bound])
    )
    unique_original, unique_indices = np.unique(
        original_knots,
        return_index=True,
    )
    unique_display = display_knots[unique_indices]
    unique_display = np.maximum.accumulate(unique_display)
    unique_display[-1] = upper_bound

    table = pd.DataFrame(
        {
            "axis": axis,
            "quantile": quantiles,
            "original_coordinate": original,
            "equalized_target_coordinate": target,
            "display_coordinate": display,
            "equalization_strength": strength,
            "target_margin_fraction": margin_fraction,
        }
    )
    return (
        {
            "original_knots": unique_original,
            "display_knots": unique_display,
        },
        table,
    )


def apply_axis_warp(
    values: np.ndarray,
    warp: dict[str, np.ndarray],
) -> np.ndarray:
    return np.interp(
        np.asarray(values, dtype=float),
        warp["original_knots"],
        warp["display_knots"],
    )


def warp_frame_coordinates(
    frame: pd.DataFrame,
    x_warp: dict[str, np.ndarray],
    y_warp: dict[str, np.ndarray],
    *,
    preserve_layout_coordinates: bool = False,
) -> pd.DataFrame:
    output = frame.copy()
    if preserve_layout_coordinates:
        for column in ("layout_x", "layout_y"):
            if column in output:
                output[f"original_{column}"] = output[column]

    for column in ("x", "layout_x", "source_x", "target_x"):
        if column not in output:
            continue
        numeric = pd.to_numeric(output[column], errors="coerce")
        mask = numeric.notna()
        output.loc[mask, column] = apply_axis_warp(
            numeric.loc[mask].to_numpy(dtype=float),
            x_warp,
        )
    for column in ("y", "layout_y", "source_y", "target_y"):
        if column not in output:
            continue
        numeric = pd.to_numeric(output[column], errors="coerce")
        mask = numeric.notna()
        output.loc[mask, column] = apply_axis_warp(
            numeric.loc[mask].to_numpy(dtype=float),
            y_warp,
        )
    return output


def warp_density_layer(
    density: np.ndarray,
    bounds: np.ndarray,
    x_warp: dict[str, np.ndarray],
    y_warp: dict[str, np.ndarray],
) -> np.ndarray:
    """Resample a frozen density layer into the shared display coordinates."""
    density = np.asarray(density, dtype=np.float32)
    height, width = density.shape
    display_x = np.linspace(bounds[0], bounds[1], width)
    display_y = np.linspace(bounds[2], bounds[3], height)
    source_x = np.interp(
        display_x,
        x_warp["display_knots"],
        x_warp["original_knots"],
    )
    source_y = np.interp(
        display_y,
        y_warp["display_knots"],
        y_warp["original_knots"],
    )
    source_columns = (
        (source_x - bounds[0])
        / max(bounds[1] - bounds[0], 1e-12)
        * (width - 1)
    ).astype(np.float32)
    source_rows = (
        (source_y - bounds[2])
        / max(bounds[3] - bounds[2], 1e-12)
        * (height - 1)
    ).astype(np.float32)
    row_grid, column_grid = np.meshgrid(
        source_rows,
        source_columns,
        indexing="ij",
    )
    warped = map_coordinates(
        density,
        [row_grid, column_grid],
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    )
    return np.clip(warped, 0.0, 1.0).astype(np.float32)


def marginal_histogram_cv(
    values: np.ndarray,
    *,
    bins: int = 80,
) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    counts, _ = np.histogram(values, bins=bins)
    return float(counts.std() / max(counts.mean(), 1e-12))


def build_equalized_display_layout(
    *,
    network: pd.DataFrame,
    density,
    g2_points: pd.DataFrame,
    g3_stipple: pd.DataFrame,
    representative: pd.DataFrame,
    g1_examples: pd.DataFrame,
) -> dict[str, object]:
    """Apply one display-only coordinate warp to every plotted layer."""
    parts = split_network(network)
    node_x = np.concatenate(
        [
            parts["g0"]["x"].to_numpy(dtype=float),
            parts["g1"]["layout_x"].to_numpy(dtype=float),
            g2_points["layout_x"].to_numpy(dtype=float),
        ]
    )
    node_y = np.concatenate(
        [
            parts["g0"]["y"].to_numpy(dtype=float),
            parts["g1"]["layout_y"].to_numpy(dtype=float),
            g2_points["layout_y"].to_numpy(dtype=float),
        ]
    )
    bounds = np.asarray(density["bounds"], dtype=float)
    x_warp, x_table = build_axis_warp(
        node_x,
        float(bounds[0]),
        float(bounds[1]),
        axis="x",
    )
    y_warp, y_table = build_axis_warp(
        node_y,
        float(bounds[2]),
        float(bounds[3]),
        axis="y",
    )

    display_node_x = apply_axis_warp(node_x, x_warp)
    display_node_y = apply_axis_warp(node_y, y_warp)
    uniformity = pd.DataFrame(
        [
            {
                "axis": "x",
                "marginal_histogram_cv_before": marginal_histogram_cv(node_x),
                "marginal_histogram_cv_after": marginal_histogram_cv(
                    display_node_x
                ),
            },
            {
                "axis": "y",
                "marginal_histogram_cv_before": marginal_histogram_cv(node_y),
                "marginal_histogram_cv_after": marginal_histogram_cv(
                    display_node_y
                ),
            },
        ]
    )
    uniformity["relative_cv_reduction"] = (
        uniformity["marginal_histogram_cv_before"]
        - uniformity["marginal_histogram_cv_after"]
    ) / uniformity["marginal_histogram_cv_before"]

    display_density = {
        "bounds": bounds.copy(),
        "g3_component_normalized_intensity": warp_density_layer(
            density["g3_component_normalized_intensity"],
            bounds,
            x_warp,
            y_warp,
        ),
    }
    display_network = warp_frame_coordinates(network, x_warp, y_warp)
    display_network["display_layout"] = "soft_density_equalized_v12"
    display_g2 = warp_frame_coordinates(
        g2_points,
        x_warp,
        y_warp,
        preserve_layout_coordinates=True,
    )
    display_g3_stipple = warp_frame_coordinates(
        g3_stipple,
        x_warp,
        y_warp,
        preserve_layout_coordinates=True,
    )
    display_representative = warp_frame_coordinates(
        representative,
        x_warp,
        y_warp,
    )
    display_g1_examples = warp_frame_coordinates(
        g1_examples,
        x_warp,
        y_warp,
    )
    return {
        "display_network": display_network,
        "display_density": display_density,
        "display_g2_points": display_g2,
        "display_g3_stipple": display_g3_stipple,
        "display_representative": display_representative,
        "display_g1_examples": display_g1_examples,
        "layout_equalization_knots": pd.concat(
            [x_table, y_table],
            ignore_index=True,
        ),
        "layout_uniformity_summary": uniformity,
    }


def grid_occupancy_metrics(
    x_values: np.ndarray,
    y_values: np.ndarray,
    *,
    bins: int = 120,
    support_mask: np.ndarray | None = None,
) -> tuple[dict[str, float], np.ndarray]:
    """Summarize two-dimensional occupancy in a fixed normalized frame."""
    histogram, _, _ = np.histogram2d(
        np.clip(np.asarray(x_values, dtype=float), 0.0, 1.0),
        np.clip(np.asarray(y_values, dtype=float), 0.0, 1.0),
        bins=bins,
        range=((0.0, 1.0), (0.0, 1.0)),
    )
    occupied = histogram > 0
    if support_mask is None:
        support_mask = binary_fill_holes(
            binary_closing(
                binary_dilation(occupied, iterations=1),
                iterations=4,
            )
        )
    support_values = histogram[support_mask]
    mean = float(support_values.mean()) if len(support_values) else 0.0
    metrics = {
        "support_cell_count": int(support_mask.sum()),
        "occupied_support_cell_count": int(
            (occupied & support_mask).sum()
        ),
        "occupied_support_fraction": float(
            (occupied & support_mask).sum() / max(support_mask.sum(), 1)
        ),
        "empty_support_fraction": float(
            ((~occupied) & support_mask).sum()
            / max(support_mask.sum(), 1)
        ),
        "support_count_cv": float(
            support_values.std() / max(mean, 1e-12)
        ),
    }
    return metrics, support_mask


def relax_g2_display_layout(
    g2_points: pd.DataFrame,
    network: pd.DataFrame,
    bounds: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Spread dense G2 display clouds into nearby low-occupancy regions.

    The operation is deterministic, bounded, and applied only to the largest
    G0 component. It changes no node identity, parent assignment, generation,
    or explicit G0/G1 edge. Coordinates remain presentation-only.
    """
    output = g2_points.copy()
    output["pre_relaxation_layout_x"] = output["layout_x"]
    output["pre_relaxation_layout_y"] = output["layout_y"]

    parts = split_network(network)
    core_component = int(
        parts["g0"].groupby("component_index").size().idxmax()
    )
    core_mask = output["layout_component_index"].eq(core_component).to_numpy()
    if not core_mask.any():
        raise RuntimeError("No G2 points belong to the largest G0 component")

    x_span = float(bounds[1] - bounds[0])
    y_span = float(bounds[3] - bounds[2])
    x_all = (
        output["layout_x"].to_numpy(dtype=float) - float(bounds[0])
    ) / x_span
    y_all = (
        output["layout_y"].to_numpy(dtype=float) - float(bounds[2])
    ) / y_span
    x = x_all[core_mask].copy()
    y = y_all[core_mask].copy()
    original_x = x.copy()
    original_y = y.copy()

    before_metrics, support_mask = grid_occupancy_metrics(x, y)
    iteration_rows = []
    grid_size = G2_RELAXATION_GRID_SIZE
    margin = G2_RELAXATION_BOUNDARY_MARGIN_FRACTION

    for iteration in range(1, G2_RELAXATION_ITERATIONS + 1):
        histogram, _, _ = np.histogram2d(
            x,
            y,
            bins=grid_size,
            range=((0.0, 1.0), (0.0, 1.0)),
        )
        smoothed = gaussian_filter(
            histogram.astype(np.float64),
            sigma=G2_RELAXATION_SMOOTHING_SIGMA,
            mode="nearest",
        )
        positive = smoothed[smoothed > 0]
        reference = float(np.median(positive)) if len(positive) else 1.0
        potential = np.log1p(smoothed / max(reference, 1e-12))
        gradient_x, gradient_y = np.gradient(potential)
        coordinates = np.vstack(
            [x * (grid_size - 1), y * (grid_size - 1)]
        )
        sampled_gradient_x = map_coordinates(
            gradient_x,
            coordinates,
            order=1,
            mode="nearest",
            prefilter=False,
        )
        sampled_gradient_y = map_coordinates(
            gradient_y,
            coordinates,
            order=1,
            mode="nearest",
            prefilter=False,
        )
        magnitude = np.hypot(sampled_gradient_x, sampled_gradient_y)
        direction_x = sampled_gradient_x / np.maximum(magnitude, 1e-12)
        direction_y = sampled_gradient_y / np.maximum(magnitude, 1e-12)
        step = G2_RELAXATION_STEP_FRACTION * np.tanh(2.5 * magnitude)
        candidate_x = x - step * direction_x
        candidate_y = y - step * direction_y

        delta_x = candidate_x - original_x
        delta_y = candidate_y - original_y
        displacement = np.hypot(delta_x, delta_y)
        scale = np.minimum(
            1.0,
            G2_RELAXATION_MAX_DISPLACEMENT_FRACTION
            / np.maximum(displacement, 1e-12),
        )
        x = np.clip(
            original_x + delta_x * scale,
            margin,
            1.0 - margin,
        )
        y = np.clip(
            original_y + delta_y * scale,
            margin,
            1.0 - margin,
        )
        iteration_metrics, _ = grid_occupancy_metrics(
            x,
            y,
            support_mask=support_mask,
        )
        iteration_rows.append(
            {
                "iteration": iteration,
                **iteration_metrics,
                "median_normalized_displacement": float(
                    np.median(np.hypot(x - original_x, y - original_y))
                ),
                "maximum_normalized_displacement": float(
                    np.max(np.hypot(x - original_x, y - original_y))
                ),
            }
        )

    x_all[core_mask] = x
    y_all[core_mask] = y
    output["layout_x"] = float(bounds[0]) + x_all * x_span
    output["layout_y"] = float(bounds[2]) + y_all * y_span
    output["g2_relaxation_applied"] = core_mask
    output["relaxation_displacement"] = np.hypot(
        output["layout_x"] - output["pre_relaxation_layout_x"],
        output["layout_y"] - output["pre_relaxation_layout_y"],
    )
    output["relaxation_displacement_fraction_of_diagonal"] = (
        output["relaxation_displacement"] / np.hypot(x_span, y_span)
    )
    output["layout_method"] = output["layout_method"].astype(str).where(
        ~output["g2_relaxation_applied"],
        output["layout_method"].astype(str)
        + "+bounded_density_relaxation_v12",
    )

    after_metrics, _ = grid_occupancy_metrics(
        x,
        y,
        support_mask=support_mask,
    )
    displacement = output.loc[
        output["g2_relaxation_applied"], "relaxation_displacement"
    ].to_numpy(dtype=float)
    summary = pd.DataFrame(
        [
            {
                "metric": "core_G2_node_count",
                "before": int(core_mask.sum()),
                "after": int(core_mask.sum()),
                "change": 0,
            },
            {
                "metric": "occupied_support_fraction",
                "before": before_metrics["occupied_support_fraction"],
                "after": after_metrics["occupied_support_fraction"],
                "change": (
                    after_metrics["occupied_support_fraction"]
                    - before_metrics["occupied_support_fraction"]
                ),
            },
            {
                "metric": "empty_support_fraction",
                "before": before_metrics["empty_support_fraction"],
                "after": after_metrics["empty_support_fraction"],
                "change": (
                    after_metrics["empty_support_fraction"]
                    - before_metrics["empty_support_fraction"]
                ),
            },
            {
                "metric": "support_count_cv",
                "before": before_metrics["support_count_cv"],
                "after": after_metrics["support_count_cv"],
                "change": (
                    after_metrics["support_count_cv"]
                    - before_metrics["support_count_cv"]
                ),
            },
            {
                "metric": "median_display_displacement",
                "before": 0.0,
                "after": float(np.median(displacement)),
                "change": float(np.median(displacement)),
            },
            {
                "metric": "maximum_display_displacement",
                "before": 0.0,
                "after": float(np.max(displacement)),
                "change": float(np.max(displacement)),
            },
        ]
    )
    return output, summary, pd.DataFrame(iteration_rows)


def space_id_suffix(values: pd.Series) -> np.ndarray:
    """Extract the deterministic numeric suffix used by the frozen layout."""
    return (
        values.str.split("_", n=1)
        .str[-1]
        .astype(np.int64)
        .to_numpy()
    )


def pseudo_unit(
    values: np.ndarray,
    multiplier: int,
    increment: int,
) -> np.ndarray:
    raw = (
        values.astype(np.uint64) * np.uint64(multiplier)
        + np.uint64(increment)
    ) & np.uint64(0xFFFFFFFF)
    return raw.astype(np.float64) / float(2**32)


def recover_g2_layout(
    network: pd.DataFrame,
    work: Path,
) -> pd.DataFrame:
    """Recover the exact display coordinates used for the frozen G2 density.

    This is a rendering-only operation over frozen first-observation
    derivations. It neither reapplies reaction rules nor enumerates products.
    """
    database = (
        work
        / "inputs/G0_G3_primary_release/03_primary_G0_G3/"
        "taxane_reaction_grammar_space.sqlite"
    )
    if not database.is_file():
        raise FileNotFoundError(f"Frozen chemical-space database missing: {database}")

    nodes = network[network["record_type"].eq("node")].copy()
    nodes["generation"] = pd.to_numeric(nodes["generation"]).astype(int)
    g1 = nodes[nodes["generation"].eq(1)].copy()
    for column in ("layout_x", "layout_y", "layout_component_index"):
        g1[column] = pd.to_numeric(g1[column], errors="raise")

    connection = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
    try:
        g2_edges = pd.read_sql_query(
            """
            SELECT source_space_id, target_space_id
            FROM derivation_events
            WHERE generation = 2 AND target_is_new = 1
            ORDER BY target_space_id
            """,
            connection,
        )
    finally:
        connection.close()

    if len(g2_edges) != 223823:
        raise RuntimeError(
            "Unexpected frozen G2 first-observation edge count: "
            f"{len(g2_edges):,}"
        )
    if g2_edges["target_space_id"].duplicated().any():
        raise RuntimeError("Frozen G2 first-observation targets are not unique")

    g1_lookup = g1.set_index("space_id")
    source_x = g2_edges["source_space_id"].map(g1_lookup["layout_x"])
    source_y = g2_edges["source_space_id"].map(g1_lookup["layout_y"])
    components = g2_edges["source_space_id"].map(
        g1_lookup["layout_component_index"]
    )
    if source_x.isna().any() or source_y.isna().any() or components.isna().any():
        raise RuntimeError("G2 parent coordinate or component mapping is incomplete")

    numeric_ids = space_id_suffix(g2_edges["target_space_id"])
    generation = 2
    angle_unit = pseudo_unit(
        numeric_ids,
        2654435761 + generation * 101,
        2246822519 + generation * 17,
    )
    radius_unit = pseudo_unit(
        numeric_ids,
        3266489917 + generation * 131,
        668265263 + generation * 29,
    )
    angle = 2 * np.pi * angle_unit
    radius = 3.2 + (9.2 - 3.2) * np.sqrt(radius_unit)
    source_x_array = source_x.to_numpy(dtype=float)
    source_y_array = source_y.to_numpy(dtype=float)

    return pd.DataFrame(
        {
            "space_id": g2_edges["target_space_id"],
            "parent_space_id": g2_edges["source_space_id"],
            "layout_x": source_x_array + radius * np.cos(angle),
            "layout_y": source_y_array + radius * np.sin(angle),
            "layout_component_index": components.to_numpy(dtype=np.int32),
            "layout_method": "frozen_v4_2_deterministic_child_coordinates",
            "generation": "G2",
        }
    )


def split_network(network: pd.DataFrame) -> dict[str, pd.DataFrame]:
    return base.split_network(network)


def line_segments(
    frame: pd.DataFrame,
) -> list[list[tuple[float, float]]]:
    return [
        [(row.source_x, row.source_y), (row.target_x, row.target_y)]
        for row in frame.itertuples(index=False)
    ]


def draw_network_axis(
    ax,
    network: pd.DataFrame,
    density,
    g2_points: pd.DataFrame,
    g3_stipple: pd.DataFrame,
    selected_g0: pd.DataFrame,
    selected_g1: pd.DataFrame,
) -> None:
    parts = split_network(network)
    g0 = parts["g0"]
    g1 = parts["g1"]
    g0_edges = parts["g0_edges"]
    g1_edges = parts["g1_edges"]
    bounds = np.asarray(density["bounds"], dtype=float)
    extent = (bounds[0], bounds[1], bounds[2], bounds[3])

    ax.set_facecolor(base.COLORS["white"])
    ax.imshow(
        base.intensity_rgba(
            density["g3_component_normalized_intensity"],
            base.COLORS["g3"],
            0.22,
        ),
        extent=extent,
        origin="lower",
        interpolation="bilinear",
        zorder=0,
    )
    ax.scatter(
        g3_stipple["layout_x"],
        g3_stipple["layout_y"],
        s=0.13,
        color=base.COLORS["g3"],
        alpha=0.16,
        edgecolor="none",
        rasterized=True,
        zorder=1,
    )

    core_component = int(g0.groupby("component_index").size().idxmax())
    core_g0 = g0[g0["component_index"].eq(core_component)]
    embedded_g0 = g0[~g0["component_index"].eq(core_component)]
    core_g1 = g1[g1["layout_component_index"].eq(core_component)]
    embedded_g1 = g1[~g1["layout_component_index"].eq(core_component)]
    core_g2 = g2_points[
        g2_points["layout_component_index"].eq(core_component)
    ]
    embedded_g2 = g2_points[
        ~g2_points["layout_component_index"].eq(core_component)
    ]
    core_g0_edges = g0_edges[
        g0_edges["layout_component_index"].eq(core_component)
    ]
    embedded_g0_edges = g0_edges[
        ~g0_edges["layout_component_index"].eq(core_component)
    ]
    core_g1_edges = g1_edges[
        g1_edges["layout_component_index"].eq(core_component)
    ]
    embedded_g1_edges = g1_edges[
        ~g1_edges["layout_component_index"].eq(core_component)
    ]

    ax.scatter(
        embedded_g2["layout_x"],
        embedded_g2["layout_y"],
        s=G2_NODE_AREA,
        color=G2_DISPLAY_COLOR,
        alpha=0.14,
        edgecolor="none",
        rasterized=True,
        zorder=2,
    )
    ax.scatter(
        core_g2["layout_x"],
        core_g2["layout_y"],
        s=G2_NODE_AREA,
        color=G2_DISPLAY_COLOR,
        alpha=0.27,
        edgecolor="none",
        rasterized=True,
        zorder=2,
    )

    edge_layers = [
        (embedded_g1_edges, base.COLORS["blue"], 0.055, 0.006, 3),
        (embedded_g0_edges, base.COLORS["ink"], 0.12, 0.032, 4),
        (core_g1_edges, base.COLORS["blue"], 0.085, 0.022, 3),
        (core_g0_edges, base.COLORS["ink"], 0.23, 0.145, 5),
    ]
    for frame, color, width, alpha, zorder in edge_layers:
        ax.add_collection(
            LineCollection(
                line_segments(frame),
                colors=color,
                linewidths=width,
                alpha=alpha,
                zorder=zorder,
            )
        )

    marker_size = G0_G1_NODE_AREA
    ax.scatter(
        embedded_g1["layout_x"],
        embedded_g1["layout_y"],
        s=marker_size,
        color=base.COLORS["blue"],
        alpha=0.20,
        edgecolor="none",
        zorder=4,
    )
    ax.scatter(
        embedded_g0["x"],
        embedded_g0["y"],
        s=marker_size,
        color="#d92f3d",
        alpha=0.78,
        edgecolor="white",
        linewidth=0.055,
        zorder=7,
    )
    ax.scatter(
        core_g1["layout_x"],
        core_g1["layout_y"],
        s=marker_size,
        color=base.COLORS["blue"],
        alpha=0.73,
        edgecolor="none",
        zorder=5,
    )
    ax.scatter(
        core_g0["x"],
        core_g0["y"],
        s=marker_size,
        color="#d92f3d",
        alpha=1.0,
        edgecolor="white",
        linewidth=0.10,
        zorder=8,
    )

    for number, row in enumerate(
        selected_g0.itertuples(index=False), start=1
    ):
        ax.scatter(
            [row.x],
            [row.y],
            s=42,
            facecolor=base.COLORS["white"],
            edgecolor=base.COLORS["ink"],
            linewidth=0.70,
            alpha=0.90,
            zorder=12,
        )
        ax.text(
            row.x,
            row.y,
            str(number),
            ha="center",
            va="center",
            fontsize=4.3,
            fontweight="bold",
            color=base.COLORS["ink"],
            zorder=13,
        )

    for row in selected_g1.itertuples(index=False):
        ax.scatter(
            [float(row.layout_x)],
            [float(row.layout_y)],
            s=38,
            facecolor=base.COLORS["white"],
            edgecolor=base.COLORS["blue"],
            linewidth=0.75,
            alpha=0.94,
            zorder=12,
        )
        ax.text(
            float(row.layout_x),
            float(row.layout_y),
            str(row.network_marker),
            ha="center",
            va="center",
            fontsize=4.3,
            fontweight="bold",
            color=base.COLORS["blue"],
            zorder=13,
        )

    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[3], bounds[2])
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#cfd6dc")
        spine.set_linewidth(0.65)


def draw_generation_chart(ax, summary: pd.Series) -> None:
    generations = ["G0", "G1", "G2", "G3"]
    values = [int(summary.loc[generation]) for generation in generations]
    positions = np.arange(len(generations))
    bars = ax.bar(
        positions,
        values,
        width=0.62,
        color=[GENERATION_COLORS[value] for value in generations],
        edgecolor=base.COLORS["ink"],
        linewidth=0.35,
        zorder=2,
    )
    ax.set_yscale("log")
    ax.set_ylim(1e2, 1e7)
    ax.set_xticks(positions, generations)
    ax.set_ylabel("Unique structures\n(log scale)", labelpad=1)
    ax.yaxis.set_label_coords(-0.10, 0.5)
    ax.set_title("G0-G3 space expansion", pad=4)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value * 1.20,
            f"{value:,}",
            ha="center",
            va="bottom",
            fontsize=5.3,
            fontweight="bold",
            color=base.COLORS["ink"],
        )
    ax.grid(axis="y", color=base.COLORS["grid"], linewidth=0.55, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=5.7)


def draw_g1_card(ax, row) -> None:
    ax.set_axis_off()
    ax.add_patch(
        Rectangle(
            (0.01, 0.01),
            0.98,
            0.98,
            transform=ax.transAxes,
            facecolor=base.COLORS["white"],
            edgecolor="#ccd4da",
            linewidth=0.75,
        )
    )
    ax.text(
        0.055,
        0.91,
        str(row.network_marker),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        fontweight="bold",
        color=base.COLORS["white"],
        bbox={
            "boxstyle": "circle,pad=0.26",
            "facecolor": base.COLORS["blue"],
            "edgecolor": "none",
        },
    )
    ax.text(
        0.17,
        0.91,
        str(row.display_label),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        fontweight="bold",
        color=base.COLORS["ink"],
    )
    structure_ax = ax.inset_axes([0.06, 0.24, 0.88, 0.59])
    base.draw_molecule(structure_ax, row.target_smiles)
    ax.text(
        0.50,
        0.14,
        base.formula_label(row.target_formula),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.8,
        color=base.COLORS["ink"],
    )
    ax.text(
        0.50,
        0.055,
        str(row.target_space_id),
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=4.7,
        color=base.COLORS["gray"],
    )


def build_figure(data: dict[str, object], output: Path) -> list[Path]:
    network = data["display_network"]
    density = data["display_density"]
    g2_points = data["display_g2_points"]
    g3_stipple = data["display_g3_stipple"]
    selected_g0 = data["display_representative"]
    selected_g1 = data["display_g1_examples"]
    summary = (
        network[network["record_type"].eq("generation_summary")]
        .set_index("generation")["structure_count"]
        .astype(int)
    )
    edge_count = int(network["record_type"].eq("edge").sum())

    fig = plt.figure(figsize=(15.0, 10.4))
    grid = fig.add_gridspec(
        3,
        3,
        width_ratios=[2.25, 10.5, 2.25],
        height_ratios=[1.0, 1.0, 0.38],
        hspace=0.10,
        wspace=0.05,
    )
    left_axes = [fig.add_subplot(grid[index, 0]) for index in range(2)]
    network_ax = fig.add_subplot(grid[:2, 1])
    right_axes = [fig.add_subplot(grid[index, 2]) for index in range(2)]
    card_axes = [left_axes[0], left_axes[1], right_axes[0], right_axes[1]]

    draw_network_axis(
        network_ax,
        network,
        density,
        g2_points,
        g3_stipple,
        selected_g0,
        selected_g1,
    )
    for index, (card_ax, row) in enumerate(
        zip(card_axes, selected_g0.itertuples(index=False)), start=1
    ):
        base.draw_molecule_card(card_ax, row, index)
        card_anchor = (0.99, 0.50) if index <= 2 else (0.01, 0.50)
        connection = ConnectionPatch(
            xyA=(row.x, row.y),
            coordsA=network_ax.transData,
            xyB=card_anchor,
            coordsB=card_ax.transAxes,
            arrowstyle="-",
            linewidth=0.42,
            color=base.COLORS["gray"],
            alpha=0.28,
            zorder=10,
        )
        fig.add_artist(connection)

    network_ax.text(
        0.015,
        0.985,
        f"{int(summary.sum()):,} unique structures  |  "
        f"{edge_count:,} explicit G0/G1 edge records",
        transform=network_ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.6,
        color=base.COLORS["ink"],
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": base.COLORS["white"],
            "edgecolor": "#d3d9de",
            "linewidth": 0.60,
            "alpha": 0.96,
        },
        zorder=20,
    )
    network_ax.text(
        0.015,
        0.942,
        "G0 648  |  G1 15,801  |  G2 223,823  |  "
        "G3 2,362,766 (exploratory)",
        transform=network_ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color=base.COLORS["gray"],
        bbox={
            "boxstyle": "round,pad=0.32",
            "facecolor": base.COLORS["white"],
            "edgecolor": "#e0e4e8",
            "linewidth": 0.5,
            "alpha": 0.94,
        },
        zorder=20,
    )
    generation_ax = fig.add_subplot(grid[2, 0])
    draw_generation_chart(generation_ax, summary)

    g1_grid = grid[2, 1:].subgridspec(1, 3, wspace=0.07)
    g1_axes = [fig.add_subplot(g1_grid[0, index]) for index in range(3)]
    for ax, row in zip(g1_axes, selected_g1.itertuples(index=False)):
        draw_g1_card(ax, row)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#d92f3d",
            markeredgecolor="white",
            markeredgewidth=0.35,
            markersize=4.6,
            label="Known taxane seed (G0)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=base.COLORS["blue"],
            markeredgecolor="none",
            markersize=4.6,
            label="One-rule intermediate (G1)",
        ),
        Line2D(
            [0],
            [0],
            color=base.COLORS["ink"],
            linewidth=0.8,
            alpha=0.58,
            label="Established G0 edge",
        ),
        Line2D(
            [0],
            [0],
            color=base.COLORS["blue"],
            linewidth=0.7,
            alpha=0.42,
            label="G0 to G1 derivation",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=G2_DISPLAY_COLOR,
            markeredgecolor="none",
            markersize=4.6,
            label="Two-rule intermediate (G2)",
        ),
        Patch(
            facecolor=base.COLORS["g3"],
            edgecolor="none",
            alpha=0.52,
            label="G3 descendant-density background",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=6,
        frameon=False,
        fontsize=5.8,
        bbox_to_anchor=(0.5, 0.014),
        columnspacing=1.35,
        handlelength=1.5,
    )
    fig.suptitle(
        "Multiscale organization of T1-derived taxane chemical space",
        fontsize=15.0,
        fontweight="bold",
        y=0.988,
    )
    fig.text(
        0.5,
        0.959,
        "Node-resolved G0-G2 space, explicit G0/G1 topology, "
        "stippled G3 density, and representative molecular products",
        ha="center",
        va="center",
        fontsize=7.6,
        color=base.COLORS["gray"],
    )
    fig.subplots_adjust(
        left=0.030,
        right=0.985,
        top=0.925,
        bottom=0.060,
    )
    fig.patch.set_facecolor(base.COLORS["white"])

    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    stem = "Figure_2_T1_Taxane_Chemical_Space_V12_Balanced_G2"
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 400}),
    ):
        path = figure_dir / f"{stem}{suffix}"
        fig.savefig(path, **kwargs)
        paths.append(path)
    plt.close(fig)
    return paths


def validate(data: dict[str, object]) -> pd.DataFrame:
    network = data["network"]
    summary = (
        network[network["record_type"].eq("generation_summary")]
        .set_index("generation")["structure_count"]
        .astype(int)
    )
    observed = {
        "G0_unique_structures": int(summary.loc["G0"]),
        "G1_unique_structures": int(summary.loc["G1"]),
        "G2_unique_structures": int(summary.loc["G2"]),
        "G3_unique_structures": int(summary.loc["G3"]),
        "total_unique_structures": int(summary.sum()),
        "explicit_G0_G1_edge_records": int(
            network["record_type"].eq("edge").sum()
        ),
        "G0_molecular_anchor_count": len(data["representative"]),
        "G1_molecular_example_count": len(data["g1_examples"]),
        "G2_recovered_layout_points": len(data["g2_points"]),
        "G2_unique_layout_space_ids": int(
            data["g2_points"]["space_id"].nunique()
        ),
        "G2_unique_parent_links": int(
            data["g2_points"][
                ["parent_space_id", "space_id"]
            ].drop_duplicates().shape[0]
        ),
        "G2_finite_coordinate_rows": int(
            np.isfinite(
                data["g2_points"][["layout_x", "layout_y"]].to_numpy(
                    dtype=float
                )
            ).all(axis=1).sum()
        ),
        "G3_density_stipple_marks": len(data["g3_stipple"]),
        "G3_stipple_non_structure_marks": int(
            (~data["g3_stipple"]["represents_individual_structure"]).sum()
        ),
        "G3_stipple_finite_coordinate_rows": int(
            np.isfinite(
                data["g3_stipple"][["layout_x", "layout_y"]].to_numpy(
                    dtype=float
                )
            ).all(axis=1).sum()
        ),
        "G2_to_G0_marker_diameter_ratio_x1000": int(
            round(np.sqrt(G2_NODE_AREA / G0_G1_NODE_AREA) * 1000)
        ),
        "G0_to_G1_marker_area_ratio_x1000": 1000,
        "layout_uniformity_axes_improved": int(
            (
                data["layout_uniformity_summary"][
                    "marginal_histogram_cv_after"
                ]
                < data["layout_uniformity_summary"][
                    "marginal_histogram_cv_before"
                ]
            ).sum()
        ),
        "display_G2_finite_coordinate_rows": int(
            np.isfinite(
                data["display_g2_points"][
                    ["layout_x", "layout_y"]
                ].to_numpy(dtype=float)
            ).all(axis=1).sum()
        ),
        "G2_relaxed_core_node_count": int(
            data["display_g2_points"]["g2_relaxation_applied"].sum()
        ),
        "G2_nonrelaxed_nodes_unchanged": int(
            np.isclose(
                data["display_g2_points"].loc[
                    ~data["display_g2_points"][
                        "g2_relaxation_applied"
                    ],
                    "relaxation_displacement",
                ].to_numpy(dtype=float),
                0.0,
            ).all()
        ),
        "G2_empty_support_fraction_reduced": int(
            float(
                data["g2_relaxation_summary"].loc[
                    data["g2_relaxation_summary"]["metric"].eq(
                        "empty_support_fraction"
                    ),
                    "after",
                ].iloc[0]
            )
            < float(
                data["g2_relaxation_summary"].loc[
                    data["g2_relaxation_summary"]["metric"].eq(
                        "empty_support_fraction"
                    ),
                    "before",
                ].iloc[0]
            )
        ),
        "G2_support_count_cv_reduced": int(
            float(
                data["g2_relaxation_summary"].loc[
                    data["g2_relaxation_summary"]["metric"].eq(
                        "support_count_cv"
                    ),
                    "after",
                ].iloc[0]
            )
            < float(
                data["g2_relaxation_summary"].loc[
                    data["g2_relaxation_summary"]["metric"].eq(
                        "support_count_cv"
                    ),
                    "before",
                ].iloc[0]
            )
        ),
        "G2_relaxation_within_displacement_cap": int(
            data["display_g2_points"][
                "relaxation_displacement_fraction_of_diagonal"
            ].max()
            <= G2_RELAXATION_MAX_DISPLACEMENT_FRACTION + 1e-12
        ),
    }
    expected = {
        "G0_unique_structures": 648,
        "G1_unique_structures": 15801,
        "G2_unique_structures": 223823,
        "G3_unique_structures": 2362766,
        "total_unique_structures": 2603038,
        "explicit_G0_G1_edge_records": 17661,
        "G0_molecular_anchor_count": 4,
        "G1_molecular_example_count": 3,
        "G2_recovered_layout_points": 223823,
        "G2_unique_layout_space_ids": 223823,
        "G2_unique_parent_links": 223823,
        "G2_finite_coordinate_rows": 223823,
        "G3_density_stipple_marks": G3_STIPPLE_COUNT,
        "G3_stipple_non_structure_marks": G3_STIPPLE_COUNT,
        "G3_stipple_finite_coordinate_rows": G3_STIPPLE_COUNT,
        "G2_to_G0_marker_diameter_ratio_x1000": 500,
        "G0_to_G1_marker_area_ratio_x1000": 1000,
        "layout_uniformity_axes_improved": 2,
        "display_G2_finite_coordinate_rows": 223823,
        "G2_relaxed_core_node_count": 101943,
        "G2_nonrelaxed_nodes_unchanged": 1,
        "G2_empty_support_fraction_reduced": 1,
        "G2_support_count_cv_reduced": 1,
        "G2_relaxation_within_displacement_cap": 1,
    }
    rows = []
    for claim, expected_value in expected.items():
        value = observed[claim]
        rows.append(
            {
                "claim_id": claim,
                "expected_frozen_value": expected_value,
                "observed_value": value,
                "status": "PASS" if value == expected_value else "FAIL",
            }
        )
    audit = pd.DataFrame(rows)
    if not audit["status"].eq("PASS").all():
        raise RuntimeError(
            "Numerical audit failed:\n"
            + audit[audit["status"].ne("PASS")].to_string(index=False)
        )
    if data["g1_examples"]["layout_x"].isna().any():
        raise RuntimeError("One or more selected G1 examples lack layout data")
    bounds = np.asarray(data["density"]["bounds"], dtype=float)
    g2_points = data["g2_points"]
    if not (
        g2_points["layout_x"].between(bounds[0], bounds[1]).all()
        and g2_points["layout_y"].between(bounds[2], bounds[3]).all()
    ):
        raise RuntimeError("Recovered G2 display points exceed frozen layout bounds")
    display_g2 = data["display_g2_points"]
    if not (
        display_g2["layout_x"].between(bounds[0], bounds[1]).all()
        and display_g2["layout_y"].between(bounds[2], bounds[3]).all()
    ):
        raise RuntimeError("Equalized G2 display points exceed layout bounds")
    return audit


def write_release(
    data: dict[str, object],
    audit: pd.DataFrame,
    output: Path,
    script_path: Path,
    figure_paths: list[Path],
) -> None:
    source_dir = output / "source_data"
    workflow_dir = output / "workflow"
    source_dir.mkdir(parents=True, exist_ok=True)
    workflow_dir.mkdir(parents=True, exist_ok=True)

    generation_summary = data["network"][
        data["network"]["record_type"].eq("generation_summary")
    ][["generation", "structure_count"]].copy()
    generation_summary["structure_count"] = (
        pd.to_numeric(generation_summary["structure_count"]).astype(int)
    )
    generation_summary["rendering_mode"] = generation_summary["generation"].map(
        {
            "G0": "node_resolved_shared_density_equalized_display",
            "G1": "node_resolved_shared_density_equalized_display",
            "G2": (
                "node_resolved_shared_density_equalized_and_"
                "bounded_locally_relaxed_display"
            ),
            "G3": (
                "equalized_component_normalized_density_with_"
                "nonstructural_stipple_no_contours_exploratory"
            ),
        }
    )
    generation_summary.to_csv(
        source_dir / "Figure_2_generation_summary.tsv",
        sep="\t",
        index=False,
    )
    data["representative"].to_csv(
        source_dir / "Figure_2_G0_molecular_anchors.tsv",
        sep="\t",
        index=False,
    )
    data["g1_examples"].to_csv(
        source_dir / "Figure_2_G1_molecular_examples.tsv",
        sep="\t",
        index=False,
    )
    data["display_representative"].to_csv(
        source_dir / "Figure_2_G0_molecular_anchors_display.tsv",
        sep="\t",
        index=False,
    )
    data["display_g1_examples"].to_csv(
        source_dir / "Figure_2_G1_molecular_examples_display.tsv",
        sep="\t",
        index=False,
    )
    data["display_g2_points"].to_csv(
        source_dir / "Figure_2_G2_display_layout.tsv",
        sep="\t",
        index=False,
    )
    data["display_g3_stipple"].to_csv(
        source_dir / "Figure_2_G3_density_stipple_display.tsv",
        sep="\t",
        index=False,
    )
    shutil.copy2(
        data["paths"]["network"],
        source_dir / "Figure_2_complete_network_frozen_source.tsv",
    )
    data["display_network"].to_csv(
        source_dir / "Figure_2_complete_network_display_source.tsv",
        sep="\t",
        index=False,
    )
    np.savez_compressed(
        source_dir / "Figure_2_G3_density_layer_frozen.npz",
        bounds=np.asarray(data["density"]["bounds"], dtype=float),
        g3_component_normalized_intensity=np.asarray(
            data["density"]["g3_component_normalized_intensity"],
            dtype=np.float32,
        ),
    )
    np.savez_compressed(
        source_dir / "Figure_2_G3_density_layer_display.npz",
        bounds=np.asarray(data["display_density"]["bounds"], dtype=float),
        g3_component_normalized_intensity=np.asarray(
            data["display_density"]["g3_component_normalized_intensity"],
            dtype=np.float32,
        ),
    )
    data["layout_equalization_knots"].to_csv(
        source_dir / "Figure_2_layout_equalization_knots.tsv",
        sep="\t",
        index=False,
    )
    data["layout_uniformity_summary"].to_csv(
        source_dir / "Figure_2_layout_uniformity_summary.tsv",
        sep="\t",
        index=False,
    )
    data["g2_relaxation_summary"].to_csv(
        source_dir / "Figure_2_G2_relaxation_summary.tsv",
        sep="\t",
        index=False,
    )
    data["g2_relaxation_iterations"].to_csv(
        source_dir / "Figure_2_G2_relaxation_iterations.tsv",
        sep="\t",
        index=False,
    )
    audit.to_csv(output / "NUMERICAL_AUDIT_V12.tsv", sep="\t", index=False)
    shutil.copy2(script_path, workflow_dir / script_path.name)
    shutil.copy2(
        script_path.parent / "build_figure2_figure3_v6.py",
        workflow_dir / "build_figure2_figure3_v6.py",
    )

    caption = """# Figure 2 | Multiscale organization of T1-derived taxane chemical space

Complete integrated representation of the frozen T1-derived taxane reaction-grammar space on a pure-white background. All 648 known taxanes (G0, red), 15,801 unique one-rule products (G1, blue), and 223,823 unique two-rule products (G2, neutral gray) are represented as individual nodes. Ordinary G0 and G1 nodes retain identical marker size, whereas G2 marker diameter is one-half that of G1. G0 is distinguished by saturated red fill, a fine white keyline, and a higher drawing layer; the 1,670 established G0 edge records receive modest additional contrast relative to the 15,991 accepted G0-to-G1 derivation records. A softened empirical-CDF transformation was applied uniformly to G0-G2 nodes, explicit edges, molecular anchors, G3 stipple marks, and the G3 density raster. The 101,943 G2 nodes assigned to the largest established G0 component then underwent a deterministic bounded density-relaxation step to reduce local overplotting and distribute those existing nodes across nearby underoccupied display regions. Node identities, parent links, generation membership, and all explicit G0/G1 edges were unchanged; pre-relaxation coordinates and per-node displacement are retained in the source table. Within the prespecified occupancy support, the empty-cell fraction decreased from 0.1195 to 0.0070 and the local count coefficient of variation decreased from 0.9912 to 0.4465. The 2,362,766 exploratory G3 structures are shown as a component-normalized orange descendant-density background without contour lines. A reproducible 50,000-mark stipple sampled from the frozen G3 density field retains low-density context; these marks encode density and do not represent additional molecular structures. No reaction rule was reapplied and no structure was re-enumerated. Four representative G0 structures are linked to their displayed network coordinates. Lettered markers identify three accepted G1 products exemplifying hydroxylation, O-acetyl transfer, and deacetylation. The inset logarithmic bar chart compares unique structure counts across G0-G3. Layout coordinates organize topology and are not chemical-distance axes.
"""
    (output / "FIGURE_2_CAPTION_V12.md").write_text(
        caption,
        encoding="utf-8",
    )

    manifest = pd.DataFrame(
        [
            {
                "element": "frozen complete network",
                "source_data_file": (
                    "source_data/Figure_2_complete_network_frozen_source.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": "copied without scientific modification",
            },
            {
                "element": "display complete network",
                "source_data_file": (
                    "source_data/Figure_2_complete_network_display_source.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": (
                    "shared monotonic display-coordinate equalization; "
                    "membership and connectivity unchanged"
                ),
            },
            {
                "element": "G2 display node layout",
                "source_data_file": (
                    "source_data/Figure_2_G2_display_layout.tsv"
                ),
                "authoritative_input": str(data["paths"]["frozen_database"]),
                "operation": (
                    "display-only recovery of frozen deterministic coordinates "
                    "followed by bounded relaxation of the largest G0 "
                    "component; identities and parent links retained"
                ),
            },
            {
                "element": "G2 relaxation summary",
                "source_data_file": (
                    "source_data/Figure_2_G2_relaxation_summary.tsv"
                ),
                "authoritative_input": (
                    "source_data/Figure_2_G2_display_layout.tsv"
                ),
                "operation": (
                    "before/after occupancy and displacement audit"
                ),
            },
            {
                "element": "G2 relaxation iteration trace",
                "source_data_file": (
                    "source_data/Figure_2_G2_relaxation_iterations.tsv"
                ),
                "authoritative_input": (
                    "source_data/Figure_2_G2_display_layout.tsv"
                ),
                "operation": (
                    "deterministic iteration-level display-layout QC"
                ),
            },
            {
                "element": "frozen G3 density layer",
                "source_data_file": (
                    "source_data/Figure_2_G3_density_layer_frozen.npz"
                ),
                "authoritative_input": str(data["paths"]["density"]),
                "operation": "display-only extraction of frozen density layer",
            },
            {
                "element": "display G3 density layer",
                "source_data_file": (
                    "source_data/Figure_2_G3_density_layer_display.npz"
                ),
                "authoritative_input": str(data["paths"]["density"]),
                "operation": (
                    "same monotonic display-coordinate equalization as nodes"
                ),
            },
            {
                "element": "display G3 density stipple",
                "source_data_file": (
                    "source_data/Figure_2_G3_density_stipple_display.tsv"
                ),
                "authoritative_input": str(data["paths"]["density"]),
                "operation": (
                    "deterministic display-only sampling from frozen density; "
                    "marks are not individual structures"
                ),
            },
            {
                "element": "layout equalization",
                "source_data_file": (
                    "source_data/Figure_2_layout_equalization_knots.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": "display-only monotonic coordinate transformation",
            },
            {
                "element": "layout uniformity QC",
                "source_data_file": (
                    "source_data/Figure_2_layout_uniformity_summary.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": "display-only marginal occupancy QC",
            },
            {
                "element": "G0 molecular callouts",
                "source_data_file": (
                    "source_data/Figure_2_G0_molecular_anchors.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": "display-only extraction",
            },
            {
                "element": "G0 molecular callout display coordinates",
                "source_data_file": (
                    "source_data/Figure_2_G0_molecular_anchors_display.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": "shared monotonic display-coordinate equalization",
            },
            {
                "element": "G0-G3 generation inset",
                "source_data_file": (
                    "source_data/Figure_2_generation_summary.tsv"
                ),
                "authoritative_input": str(data["paths"]["network"]),
                "operation": "display-only logarithmic rendering",
            },
            {
                "element": "G1 molecular callouts",
                "source_data_file": (
                    "source_data/Figure_2_G1_molecular_examples.tsv"
                ),
                "authoritative_input": str(data["paths"]["projected"]),
                "operation": "display-only extraction",
            },
            {
                "element": "G1 molecular callout display coordinates",
                "source_data_file": (
                    "source_data/Figure_2_G1_molecular_examples_display.tsv"
                ),
                "authoritative_input": str(data["paths"]["projected"]),
                "operation": "shared monotonic display-coordinate equalization",
            },
        ]
    )
    manifest["source_data_sha256"] = manifest["source_data_file"].map(
        lambda value: sha256(output / value)
    )
    manifest["rendering_script"] = f"workflow/{script_path.name}"
    manifest.to_csv(
        output / "FIGURE_2_SOURCE_MANIFEST_V12.tsv",
        sep="\t",
        index=False,
    )

    summary = {
        "release": "Figure 2 V12",
        "scientific_recalculation_performed": False,
        "visual_encoding": {
            "G0": "node_resolved_saturated_red_with_white_keyline",
            "G1": "node_resolved",
            "G2": (
                "neutral_gray_node_resolved_half_G0_G1_marker_diameter_"
                "with_bounded_core_relaxation"
            ),
            "G3": (
                "equalized_component_normalized_density_without_contours"
            ),
        },
        "network_background": "#FFFFFF",
        "G0_G1_ordinary_marker_area_equal": True,
        "G0_edge_visual_emphasis": "modest",
        "G2_display_color": G2_DISPLAY_COLOR,
        "G2_relaxation": {
            "grid_size": G2_RELAXATION_GRID_SIZE,
            "iterations": G2_RELAXATION_ITERATIONS,
            "smoothing_sigma": G2_RELAXATION_SMOOTHING_SIGMA,
            "step_fraction": G2_RELAXATION_STEP_FRACTION,
            "maximum_normalized_displacement": (
                G2_RELAXATION_MAX_DISPLACEMENT_FRACTION
            ),
            "applied_core_node_count": int(
                data["display_g2_points"][
                    "g2_relaxation_applied"
                ].sum()
            ),
        },
        "G3_density_contours_drawn": False,
        "layout_equalization_strength": LAYOUT_EQUALIZATION_STRENGTH,
        "G3_density_stipple_marks": G3_STIPPLE_COUNT,
        "G3_stipple_represents_individual_structures": False,
        "audit_status": "PASS",
        "figures": [
            {
                "path": str(path.relative_to(output)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in figure_paths
        ],
    }
    (output / "BUILD_SUMMARY_V12.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    (output / "README.md").write_text(
        """# Figure 2 V12

This display-only release integrates a compact G0-G3 logarithmic generation
chart and three released G1 molecular examples around the complete network.
G0 and G1 use identical ordinary-node marker size; G2 marker diameter is
one-half that of G1 and is rendered in neutral gray. G0 is made more visible
through saturation, a fine white keyline, and modestly stronger G0 edge
contrast without changing node size. A shared softened density-equalizing
display transform is applied to every plotted layer. Existing G2 nodes in the
largest G0 component then undergo bounded deterministic local relaxation to
reduce empty internal display regions and overplotting. Original coordinates,
parent links, and per-node displacement remain in the source table. The
network background is pure white. G3 remains a density background without
contours, supplemented by reproducible density-derived stippling; these marks
are not individual molecular structures. No reaction rule, structure,
derivation, descriptor, or benchmark was recalculated.
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base.configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    data = load_inputs(args.work)
    audit = validate(data)
    figure_paths = build_figure(data, args.output)
    write_release(
        data,
        audit,
        args.output,
        Path(__file__).resolve(),
        figure_paths,
    )
    print(f"Built Figure 2 V12 at {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
