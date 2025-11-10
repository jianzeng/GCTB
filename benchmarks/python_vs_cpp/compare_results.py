#!/usr/bin/env python3
"""
Compare two GCTB output files (e.g., .parRes or .snpRes) and report the maximum
absolute difference across common numeric columns.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep=r"\s+", engine="python")


def numeric_columns(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    numeric = {}
    for col in df.columns:
        try:
            values = pd.to_numeric(df[col], errors="coerce")
        except Exception:
            continue
        if values.notna().any():
            numeric[col] = values.to_numpy(dtype=float)
    return numeric


def compare_tables(py_df: pd.DataFrame, cpp_df: pd.DataFrame) -> Tuple[float, Dict[str, float]]:
    py_numeric = numeric_columns(py_df)
    cpp_numeric = numeric_columns(cpp_df)
    shared_cols = sorted(set(py_numeric) & set(cpp_numeric))
    if not shared_cols:
        raise ValueError("No shared numeric columns to compare.")

    column_diffs = {}
    max_diff = 0.0
    for col in shared_cols:
        py_values = py_numeric[col]
        cpp_values = cpp_numeric[col]
        if py_values.shape != cpp_values.shape:
            raise ValueError(f"Column shape mismatch for '{col}': {py_values.shape} vs {cpp_values.shape}")
        diff = np.abs(py_values - cpp_values)
        column_diffs[col] = float(diff.max())
        max_diff = max(max_diff, column_diffs[col])
    return max_diff, column_diffs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", required=True, type=Path, help="Path to the Python-generated file.")
    parser.add_argument("--cpp", required=True, type=Path, help="Path to the C++-generated file.")
    parser.add_argument("--tolerance", type=float, default=1e-6, help="Absolute tolerance for PASS status.")
    parser.add_argument("--label", default=None, help="Label used in the summary output.")
    parser.add_argument("--json", type=Path, default=None, help="Optional path to write JSON summary.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    py_df = load_table(args.python)
    cpp_df = load_table(args.cpp)
    max_diff, column_diffs = compare_tables(py_df, cpp_df)

    summary = {
        "label": args.label or args.python.stem,
        "python_file": str(args.python),
        "cpp_file": str(args.cpp),
        "max_abs_diff": max_diff,
        "column_max_abs_diff": column_diffs,
        "status": "PASS" if max_diff <= args.tolerance else "FAIL",
        "tolerance": args.tolerance,
    }

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with args.json.open("w") as fh:
            json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


