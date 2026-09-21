"""Shared paths and small helpers for the EJH reproducibility scripts.

Set ``EJH_PROJECT_ROOT`` to a writable analysis directory and
``EJH_UPSTREAM_ROOT`` to the directory containing legally obtained upstream
source data. Institutional row-level data are never distributed here.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
ROOT = Path(os.environ.get("EJH_PROJECT_ROOT", REPOSITORY_ROOT)).resolve()
UPSTREAM = Path(os.environ.get("EJH_UPSTREAM_ROOT", ROOT / "upstream")).resolve()
LI = Path(os.environ.get("EJH_LI_ROOT", UPSTREAM / "li")).resolve()
RESULTS = ROOT / "results"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "figures"


def ensure_dirs() -> None:
    for path in (RESULTS, TABLES, REPORTS, FIGURES):
        path.mkdir(parents=True, exist_ok=True)


def rd(path: str | Path, **kwargs) -> pd.DataFrame:
    path = Path(path)
    if "sep" not in kwargs and path.suffix.lower() in {".tsv", ".txt"}:
        kwargs["sep"] = "\t"
    return pd.read_csv(path, **kwargs)


def say(message: object) -> None:
    print(message, flush=True)


def direction(a: float, b: float) -> bool:
    if pd.isna(a) or pd.isna(b) or a == 0 or b == 0:
        return False
    return bool(np.sign(a) == np.sign(b))


def bh_fdr(p_values) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return result
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def fmt(value, digits: int = 3) -> str:
    if value is None or (isinstance(value, (float, np.floating)) and not np.isfinite(value)):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number != 0 and abs(number) < 10 ** (-digits):
        return f"{number:.{digits}e}"
    return f"{number:.{digits}f}"
