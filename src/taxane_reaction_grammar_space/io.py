from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path, **kwargs: Any) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, dtype=str, **kwargs)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame.to_csv(path, sep=sep, index=False)


def write_json(value: Any, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )

