from __future__ import annotations

from typing import Any

from .utils import clean_text, truthy


REVERSIBLE_VALUES = {"reversible", "both", "bidirectional", "<=>", "<->", "0"}
FORWARD_VALUES = {"forward", "left_to_right", "source_to_target", "substrate_to_product", "=>", "->", "irreversible", "1", "fwd"}
REVERSE_VALUES = {"reverse", "right_to_left", "target_to_source", "product_to_substrate", "<-", "<=", "-1", "rev", "backward", "backwards", "bwd"}


def source_direction_mode(direction: Any = "", is_reversible: Any = "") -> str:
    """Return the direction-evidence mode used by scoring and SMARTS handling.

    `reverse` means the source orientation is opposite to the builder's
    substrate-to-product convention and therefore needs correction before a
    predictive SMARTS is released. It is not the same as biochemical
    reversibility. `reversible` means both directional transforms may be used
    and should be emitted as separate rules.
    """
    direction_text = clean_text(direction).lower()
    reversible_text = clean_text(is_reversible).lower()
    if direction_text in REVERSIBLE_VALUES or truthy(reversible_text):
        return "source_reversible"
    if direction_text in FORWARD_VALUES:
        return "source_forward"
    if direction_text in REVERSE_VALUES:
        return "source_reverse"
    return "unknown"


def direction_mode_from_record(row: dict[str, Any]) -> str:
    return source_direction_mode(row.get("direction", ""), row.get("is_reversible", ""))


def is_source_reverse(direction: Any = "", is_reversible: Any = "") -> bool:
    return source_direction_mode(direction, is_reversible) == "source_reverse"


def is_source_reversible(direction: Any = "", is_reversible: Any = "") -> bool:
    return source_direction_mode(direction, is_reversible) == "source_reversible"


def normalize_direction_label(direction: Any = "", is_reversible: Any = "") -> str:
    """Return forward, reverse, reversible, or unknown.

    This is a human-readable companion to source_direction_mode().
    """
    mode = source_direction_mode(direction, is_reversible)
    if mode == "source_forward":
        return "forward"
    if mode == "source_reverse":
        return "reverse"
    if mode == "source_reversible":
        return "reversible"
    return "unknown"


def direction_qc_from_handling(direction_handling: Any, direction_variant: Any = "") -> tuple[str, str]:
    """Classify whether a left-to-right rule is safe for predictive release.

    The builder's invariant is that released SMARTS are applied left-to-right as
    substrate -> product.  Unknown source orientation is not silently trusted for
    strict releases; it is retained only as exploratory evidence with an explicit
    warning.
    """
    handling = clean_text(direction_handling).lower()
    variant = clean_text(direction_variant).lower()
    if "unknown_direction" in handling or variant == "left_to_right":
        return (
            "direction_qc_unknown_exploratory_only",
            "Source did not provide reliable direction; left-to-right orientation is retained only for exploratory use.",
        )
    if "reversed_from_source" in handling or "source_reverse_corrected" in handling:
        return (
            "direction_qc_reversed_to_substrate_product",
            "Source orientation was corrected so left side is substrate and right side is product.",
        )
    if "split_reversible" in handling or "source_reversible" in handling:
        return (
            "direction_qc_reversible_split",
            "Reversible source evidence is represented as explicit directional rule member(s).",
        )
    if "kept_forward" in handling or "source_columns_as_curated_main_pair" in handling or "exact_pair_left_to_right" in handling:
        return (
            "direction_qc_ok",
            "Source orientation supports left-to-right substrate-to-product direction.",
        )
    return (
        "direction_qc_unknown_exploratory_only",
        "Direction handling was not recognized; rule is restricted to exploratory use.",
    )
