from pathlib import Path

from article_release_builder.build_release import (
    normalize_caption,
    normalize_embedded_heading,
)
from article_release_builder.data_exports import GENERATION_FILES, TABLE_GROUP_TITLES
from article_release_builder.documents import _markdown_blocks
from article_release_builder.figures import SUPPLEMENTARY_FIGURE_GROUPS


def test_generation_file_inventory() -> None:
    assert set(GENERATION_FILES) == {"0", "1", "2", "3"}
    assert GENERATION_FILES["3"].endswith(".csv")


def test_normalize_embedded_heading() -> None:
    assert normalize_embedded_heading("# References\n\n1. Example\n", "## References").startswith(
        "## References\n"
    )


def test_markdown_blocks_join_wrapped_prose() -> None:
    blocks = _markdown_blocks(
        "## Results\n\nThis is source-wrapped\nprose in one paragraph.\n\n"
        "- item with\ncontinuation\n"
    )
    assert blocks == [
        "## Results",
        "This is source-wrapped prose in one paragraph.",
        "- item with continuation",
    ]


def test_caption_heading_is_embeddable() -> None:
    assert normalize_caption("# Figure 1 | Example\n\nCaption body.\n") == (
        "**Figure 1 | Example**\n\nCaption body."
    )


def test_curated_supplementary_inventory() -> None:
    assert set(SUPPLEMENTARY_FIGURE_GROUPS) == set(range(1, 9))
    assert SUPPLEMENTARY_FIGURE_GROUPS[1] == (1, 2)
    assert set(TABLE_GROUP_TITLES) == set(range(1, 12))
