from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from .io import ensure_dir, write_json, write_table


def _scalar(connection: sqlite3.Connection, query: str) -> int:
    return int(connection.execute(query).fetchone()[0])


def validate_generated_space(
    database_path: Path,
    output_dir: Path,
) -> dict[str, Path]:
    """Independently audit identity, endpoints, and generation invariants."""

    output_dir = ensure_dir(output_dir)
    connection = sqlite3.connect(database_path)
    checks: list[dict[str, Any]] = []

    def add(name: str, failed_rows: int, severity: str = "error") -> None:
        checks.append(
            {
                "check": name,
                "severity": severity,
                "failed_rows": int(failed_rows),
                "status": "pass" if failed_rows == 0 else "fail",
            }
        )

    add(
        "full_inchikey_is_unique",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM (
              SELECT full_inchikey FROM nodes
              GROUP BY full_inchikey HAVING COUNT(*) > 1
            )
            """,
        ),
    )
    add(
        "event_source_endpoint_exists",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events e
            LEFT JOIN nodes n ON e.source_space_id = n.space_id
            WHERE n.space_id IS NULL
            """,
        ),
    )
    add(
        "event_target_endpoint_exists",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events e
            LEFT JOIN nodes n ON e.target_space_id = n.space_id
            WHERE n.space_id IS NULL
            """,
        ),
    )
    add(
        "no_identity_self_edges",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events
            WHERE source_space_id = target_space_id
            """,
        ),
    )
    add(
        "source_belongs_to_previous_frontier",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events e
            JOIN nodes n ON e.source_space_id = n.space_id
            WHERE n.generation_first != e.generation - 1
            """,
        ),
    )
    add(
        "new_target_generation_is_event_generation",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events
            WHERE target_is_new = 1
              AND target_generation_first != generation
            """,
        ),
    )
    add(
        "existing_target_not_from_future_generation",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events
            WHERE target_is_new = 0
              AND target_generation_first > generation
            """,
        ),
    )
    add(
        "each_generated_node_has_one_discovery_event",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM (
              SELECT n.space_id,
                     SUM(CASE WHEN e.target_is_new = 1 THEN 1 ELSE 0 END) AS discoveries
              FROM nodes n
              LEFT JOIN derivation_events e ON n.space_id = e.target_space_id
              WHERE n.generation_first > 0
              GROUP BY n.space_id
              HAVING discoveries != 1
            )
            """,
        ),
    )
    add(
        "G0_nodes_are_not_marked_as_new_targets",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM derivation_events e
            JOIN nodes n ON e.target_space_id = n.space_id
            WHERE n.generation_first = 0 AND e.target_is_new = 1
            """,
        ),
    )
    add(
        "processed_application_counts_are_conserved",
        _scalar(
            connection,
            """
            SELECT COUNT(*) FROM application_audit
            WHERE application_status = 'processed'
              AND raw_product_tuple_count
                  != accepted_event_count + rejected_product_count
            """,
        ),
    )
    add(
        "accepted_application_count_matches_event_count",
        abs(
            _scalar(
                connection,
                """
                SELECT COALESCE(SUM(accepted_event_count), 0)
                FROM application_audit
                """,
            )
            - _scalar(connection, "SELECT COUNT(*) FROM derivation_events")
        ),
    )

    inventory = pd.read_sql_query(
        """
        SELECT generation_first AS generation,
               COUNT(*) AS unique_nodes,
               COUNT(DISTINCT connectivity_key) AS unique_connectivity_keys
        FROM nodes
        GROUP BY generation_first
        ORDER BY generation_first
        """,
        connection,
    )
    connection.close()

    checks_frame = pd.DataFrame(checks)
    paths = {
        "checks": output_dir / "generated_space_validation_checks.tsv",
        "generation_inventory": (
            output_dir / "generated_space_generation_inventory.tsv"
        ),
        "summary": output_dir / "generated_space_validation_summary.json",
    }
    write_table(checks_frame, paths["checks"])
    write_table(inventory, paths["generation_inventory"])
    failed_errors = checks_frame[
        (checks_frame["severity"] == "error")
        & (checks_frame["status"] == "fail")
    ]
    summary = {
        "database": str(database_path),
        "validation_status": "pass" if failed_errors.empty else "fail",
        "checks": int(len(checks_frame)),
        "passed_checks": int((checks_frame["status"] == "pass").sum()),
        "failed_error_checks": int(len(failed_errors)),
        "failed_error_names": failed_errors["check"].tolist(),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
