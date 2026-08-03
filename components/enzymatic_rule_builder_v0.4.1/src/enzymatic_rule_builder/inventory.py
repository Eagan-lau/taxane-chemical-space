from __future__ import annotations

import platform
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

from .utils import ensure_dir, sha256_file, write_json
from . import __version__


def file_inventory(paths: Iterable[str | Path]) -> pd.DataFrame:
    rows = []
    for path0 in paths:
        path = Path(path0)
        if not path.exists():
            rows.append({"path": str(path), "exists": False})
            continue
        stat = path.stat()
        rows.append({
            "path": str(path.resolve()),
            "exists": True,
            "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": sha256_file(path),
        })
    return pd.DataFrame(rows)


def write_inventory(paths: Iterable[str | Path], output_dir: str | Path) -> pd.DataFrame:
    out = ensure_dir(output_dir)
    df = file_inventory(paths)
    df.to_csv(out / "source_inventory.lock.tsv", sep="\t", index=False)
    env = {
        "package": "enzymatic-rule-builder",
        "package_version": __version__,
        "python_version": sys.version,
        "platform": platform.platform(),
    }
    try:
        import rdkit
        env["rdkit_version"] = rdkit.__version__
    except Exception:
        env["rdkit_version"] = "unavailable"
    write_json(out / "build_environment.json", env)
    return df
