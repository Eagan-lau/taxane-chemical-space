from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from typing import Any

import pandas as pd
import yaml

from .io import ensure_dir, write_json, write_table


TRACKED_DISTRIBUTIONS = (
    "rdkit",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "networkx",
    "pyyaml",
)


def _distribution_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not_installed"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def record_environment(
    output_dir: Path,
    *,
    study_config: Path | None = None,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    rows: list[dict[str, Any]] = [
        {
            "component": "python",
            "version": platform.python_version(),
            "detail": platform.python_implementation(),
        },
        {
            "component": "operating_system",
            "version": platform.platform(),
            "detail": platform.machine(),
        },
    ]
    rows.extend(
        {
            "component": name,
            "version": _distribution_version(name),
            "detail": "Python distribution",
        }
        for name in TRACKED_DISTRIBUTIONS
    )
    environment_path = output_dir / "software_environment.tsv"
    write_table(pd.DataFrame(rows), environment_path)

    input_rows: list[dict[str, Any]] = []
    if study_config is not None:
        configuration = yaml.safe_load(
            study_config.read_text(encoding="utf-8")
        )
        inputs = configuration.get("inputs", {})
        for key, value in inputs.items():
            if key.endswith("_sha256"):
                continue
            path = Path(str(value))
            declared = str(inputs.get(f"{key}_sha256", ""))
            observed = _sha256(path) if path.is_file() else ""
            input_rows.append(
                {
                    "input_name": key,
                    "path": str(path),
                    "exists": path.is_file(),
                    "size_bytes": path.stat().st_size if path.is_file() else 0,
                    "declared_sha256": declared,
                    "observed_sha256": observed,
                    "hash_matches_declared": bool(
                        declared and observed and declared == observed
                    ),
                }
            )
    input_hashes_path = output_dir / "input_file_hashes.tsv"
    write_table(
        pd.DataFrame(
            input_rows,
            columns=[
                "input_name",
                "path",
                "exists",
                "size_bytes",
                "declared_sha256",
                "observed_sha256",
                "hash_matches_declared",
            ],
        ),
        input_hashes_path,
    )

    mismatches = [
        row["input_name"]
        for row in input_rows
        if not row["hash_matches_declared"]
    ]
    summary = {
        "python_executable": sys.executable,
        "tracked_components": len(rows),
        "study_config": str(study_config) if study_config else "",
        "input_files": len(input_rows),
        "input_hash_mismatches": mismatches,
        "status": "pass" if not mismatches else "fail",
        "outputs": {
            "software_environment": str(environment_path),
            "input_file_hashes": str(input_hashes_path),
        },
    }
    summary_path = output_dir / "environment_and_input_audit.json"
    write_json(summary, summary_path)
    if mismatches:
        raise RuntimeError(
            "Input snapshot audit failed: " + ", ".join(mismatches)
        )
    return {
        "software_environment": environment_path,
        "input_file_hashes": input_hashes_path,
        "summary": summary_path,
    }
