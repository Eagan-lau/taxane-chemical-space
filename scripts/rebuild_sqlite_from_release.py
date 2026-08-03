#!/usr/bin/env python3
"""Rebuild the derivation SQLite database from restored release TSV files."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sqlite3
from typing import Any


TABLE_FILES = {
    "nodes": "chemical_space_nodes.tsv",
    "derivation_events": "derivation_events.tsv",
    "application_audit": "rule_application_audit.tsv",
    "rejection_events": "rejection_events.tsv",
    "generation_parent_progress": "generation_parent_progress.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--space-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=10_000)
    return parser.parse_args()


def convert(value: str, declared_type: str) -> Any:
    normalized_type = declared_type.upper()
    if value == "" and normalized_type in {"INTEGER", "REAL"}:
        return None
    if normalized_type == "INTEGER":
        return int(value)
    if normalized_type == "REAL":
        return float(value)
    return value


def load_table(
    connection: sqlite3.Connection,
    table: str,
    path: Path,
    batch_size: int,
) -> int:
    declared_types = {
        row[1]: row[2]
        for row in connection.execute(f"PRAGMA table_info({table})")
    }
    inserted = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        columns = next(reader)
        missing = [column for column in columns if column not in declared_types]
        if missing:
            raise ValueError(f"{table}: columns absent from schema: {missing}")
        placeholders = ",".join("?" for _ in columns)
        quoted_columns = ",".join(f'"{column}"' for column in columns)
        statement = (
            f'INSERT INTO "{table}" ({quoted_columns}) '
            f"VALUES ({placeholders})"
        )
        batch = []
        for row in reader:
            if len(row) != len(columns):
                raise ValueError(
                    f"{path}: expected {len(columns)} fields, found {len(row)}"
                )
            batch.append(
                tuple(
                    convert(value, declared_types[column])
                    for column, value in zip(columns, row)
                )
            )
            if len(batch) >= batch_size:
                connection.executemany(statement, batch)
                inserted += len(batch)
                batch.clear()
        if batch:
            connection.executemany(statement, batch)
            inserted += len(batch)
    connection.commit()
    return inserted


def main() -> int:
    args = parse_args()
    space_dir = args.space_dir.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    table_schema = space_dir / "database_table_schema.sql"
    index_schema = space_dir / "database_index_schema.sql"
    for path in (table_schema, index_schema):
        if not path.is_file():
            raise FileNotFoundError(path)
    for filename in TABLE_FILES.values():
        if not (space_dir / filename).is_file():
            raise FileNotFoundError(space_dir / filename)

    output.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    with sqlite3.connect(output) as connection:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.executescript(table_schema.read_text(encoding="utf-8"))
        for table, filename in TABLE_FILES.items():
            counts[table] = load_table(
                connection,
                table,
                space_dir / filename,
                args.batch_size,
            )
            print(f"loaded {table}: {counts[table]:,}", flush=True)
        connection.executescript(index_schema.read_text(encoding="utf-8"))
        connection.execute("ANALYZE")
        connection.commit()

    summary = json.loads(
        (space_dir / "chemical_space_build_summary.json").read_text(
            encoding="utf-8"
        )
    )
    expected_nodes = int(summary["final_node_count"])
    expected_events = int(summary["final_derivation_event_count"])
    if counts["nodes"] != expected_nodes:
        raise RuntimeError(
            f"Node count mismatch: {counts['nodes']} != {expected_nodes}"
        )
    if counts["derivation_events"] != expected_events:
        raise RuntimeError(
            "Event count mismatch: "
            f"{counts['derivation_events']} != {expected_events}"
        )
    print(f"rebuilt database: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
