#!/usr/bin/env python3
"""
Benchmark runner comparing the Python CLI and legacy C++ CLI outputs.

It executes a predefined set of workloads (BayesC and SBayesR by default),
records wall-clock runtime, and compares the resulting output files using
`compare_results.py`.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATASET_DIR = ROOT / "test" / "data"
COMPARE_SCRIPT = Path(__file__).with_name("compare_results.py")


@dataclass
class Workload:
    name: str
    python_cmd: List[str]
    cpp_cmd: List[str]
    comparisons: List[Dict[str, str]]
    prerequisites: List[str] = field(default_factory=list)
    optional: bool = False


def build_workloads(
    dataset_dir: Path,
    python_output_dir: Path,
    cpp_output_dir: Path,
    include_posthoc: bool,
    annotation_file: Optional[Path],
) -> List[Workload]:
    """Create command definitions for each workload."""
    bayes_py_prefix = python_output_dir / "bayesC"
    bayes_cpp_prefix = cpp_output_dir / "bayesC"

    sbayes_py_prefix = python_output_dir / "sbayesR"
    sbayes_cpp_prefix = cpp_output_dir / "sbayesR"

    workloads: List[Workload] = [
        Workload(
            name="bayesC",
            python_cmd=[
                "python3",
                "-m",
                "gctb.cli",
                "bayes",
                "--bfile",
                str(dataset_dir / "uk10k_chr1_1mb"),
                "--pheno",
                str(dataset_dir / "test.phen"),
                "--bayes",
                "C",
                "--chain-length",
                "1000",
                "--burnin",
                "100",
                "--thin",
                "10",
                "--pi",
                "0.01",
                "--hsq",
                "0.5",
                "--out",
                str(bayes_py_prefix),
                "--quiet",
            ],
            cpp_cmd=[
                "./scr/gctb",
                "--bfile",
                str(dataset_dir / "uk10k_chr1_1mb"),
                "--pheno",
                str(dataset_dir / "test.phen"),
                "--bayes",
                "C",
                "--chain-length",
                "1000",
                "--burnin",
                "100",
                "--thin",
                "10",
                "--pi",
                "0.01",
                "--hsq",
                "0.5",
                "--out",
                str(bayes_cpp_prefix),
            ],
            comparisons=[
                {
                    "label": "parRes",
                    "python": f"{bayes_py_prefix}.parRes",
                    "cpp": f"{bayes_cpp_prefix}.parRes",
                },
                {
                    "label": "snpRes",
                    "python": f"{bayes_py_prefix}.snpRes",
                    "cpp": f"{bayes_cpp_prefix}.snpRes",
                },
            ],
        ),
        Workload(
            name="sbayesR",
            python_cmd=[
                "python3",
                "-m",
                "gctb.cli",
                "sbayes",
                "--ldm",
                str(dataset_dir / "test_ldm.ldm.sparse"),
                "--gwas-summary",
                str(dataset_dir / "test_gwas_summary.ma"),
                "--sbayes",
                "R",
                "--chain-length",
                "1500",
                "--burnin",
                "500",
                "--thin",
                "5",
                "--out",
                str(sbayes_py_prefix),
                "--quiet",
            ],
            cpp_cmd=[
                "./scr/gctb",
                "--ldm",
                str(dataset_dir / "test_ldm.ldm.sparse"),
                "--gwas-summary",
                str(dataset_dir / "test_gwas_summary.ma"),
                "--sbayes",
                "R",
                "--chain-length",
                "1500",
                "--burnin",
                "500",
                "--thin",
                "5",
                "--out",
                str(sbayes_cpp_prefix),
            ],
            comparisons=[
                {
                    "label": "parRes",
                    "python": f"{sbayes_py_prefix}.parRes",
                    "cpp": f"{sbayes_cpp_prefix}.parRes",
                },
                {
                    "label": "snpRes",
                    "python": f"{sbayes_py_prefix}.snpRes",
                    "cpp": f"{sbayes_cpp_prefix}.snpRes",
                },
            ],
        ),
    ]

    if include_posthoc:
        if annotation_file is None or not annotation_file.exists():
            raise FileNotFoundError(
                "Post-hoc stratification requested but annotation file not provided."
            )
        strat_py_prefix = python_output_dir / "posthoc"
        strat_cpp_prefix = cpp_output_dir / "posthoc"
        workloads.append(
            Workload(
                name="posthoc_stratify",
                python_cmd=[
                    "python3",
                    "-m",
                    "gctb.cli",
                    "posthoc-stratify",
                    "--ldm",
                    str(dataset_dir / "test_ldm.ldm.sparse"),
                    "--snp-res",
                    f"{sbayes_py_prefix}.snpRes",
                    "--mcmc-prefix",
                    str(sbayes_py_prefix),
                    "--type",
                    "S",
                    "--annotation",
                    str(annotation_file),
                    "--out",
                    str(strat_py_prefix),
                    "--quiet",
                ],
                cpp_cmd=[
                    "./scr/gctb",
                    "--ldm",
                    str(dataset_dir / "test_ldm.ldm.sparse"),
                    "--snp-res",
                    f"{sbayes_cpp_prefix}.snpRes",
                    "--mcmc-sample",
                    str(sbayes_cpp_prefix),
                    "--posthoc-stratify",
                    "S",
                    "--annotation",
                    str(annotation_file),
                    "--out",
                    str(strat_cpp_prefix),
                ],
                comparisons=[
                    {
                        "label": "snpRes",
                        "python": f"{strat_py_prefix}.snpRes",
                        "cpp": f"{strat_cpp_prefix}.snpRes",
                    }
                ],
                prerequisites=["sbayesR"],
                optional=True,
            )
        )

    return workloads


def run_command(cmd: List[str], cwd: Path, log_file: Path) -> float:
    """Execute command and record wall-clock time. Raises on failure."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with log_file.open("w") as log:
        process = subprocess.run(
            cmd,
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    duration = time.perf_counter() - start
    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {process.returncode}: {' '.join(cmd)}\n"
            f"See log: {log_file}"
        )
    return duration


def run_compare(python_path: Path, cpp_path: Path, tolerance: float, out_dir: Path) -> Dict[str, float]:
    """Run the comparison helper and return its JSON summary."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / (python_path.stem + ".json")

    cmd = [
        "python3",
        str(COMPARE_SCRIPT),
        "--python",
        str(python_path),
        "--cpp",
        str(cpp_path),
        "--tolerance",
        str(tolerance),
        "--json",
        str(summary_path),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)
    with summary_path.open() as fh:
        return json.load(fh)


def ensure_cpp_binary():
    """Verify the legacy binary exists; if not, instruct user to build it."""
    binary = ROOT / "scr" / "gctb"
    if binary.exists():
        return
    raise FileNotFoundError(
        "Legacy CLI not found at scr/gctb. Build it with `make -C scr` before running benchmarks."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path("benchmarks/python_vs_cpp/results"), help="Directory to store outputs, logs, and summaries.")
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR, help="Directory containing PLINK and GWAS test data.")
    parser.add_argument("--python-only", action="store_true", help="Only run the Python workflows (skip C++).")
    parser.add_argument("--cpp-only", action="store_true", help="Only run the C++ workflows (skip Python).")
    parser.add_argument("--tolerance", type=float, default=1e-6, help="Absolute tolerance for numeric comparisons.")
    parser.add_argument("--include-posthoc", action="store_true", help="Include post-hoc stratification benchmark.")
    parser.add_argument("--annotation-file", type=Path, default=None, help="Annotation file for post-hoc stratification (required if --include-posthoc).")
    parser.add_argument("--force", action="store_true", help="Delete existing work directory before running.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_cpp_binary()

    work_dir = args.work_dir.resolve()
    python_dir = work_dir / "python"
    cpp_dir = work_dir / "cpp"
    log_dir = work_dir / "logs"
    compare_dir = work_dir / "comparisons"

    if work_dir.exists() and args.force:
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    python_dir.mkdir(parents=True, exist_ok=True)
    cpp_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    compare_dir.mkdir(parents=True, exist_ok=True)

    workloads = build_workloads(
        dataset_dir=args.dataset_dir.resolve(),
        python_output_dir=python_dir,
        cpp_output_dir=cpp_dir,
        include_posthoc=args.include_posthoc,
        annotation_file=args.annotation_file.resolve() if args.annotation_file else None,
    )

    timing_summary = []
    comparison_summary: Dict[str, List[Dict[str, float]]] = {}
    completed = set()

    for workload in workloads:
        missing_prereq = any(prereq not in completed for prereq in workload.prerequisites)
        if missing_prereq:
            print(f"[SKIP] {workload.name} - prerequisites not met.")
            continue

        print(f"[RUN] {workload.name}")
        workload_entry = {"name": workload.name, "python_time": None, "cpp_time": None}

        if not args.cpp_only:
            py_log = log_dir / f"{workload.name}_python.log"
            python_time = run_command(workload.python_cmd, ROOT, py_log)
            workload_entry["python_time"] = python_time
            print(f"  Python time: {python_time:.2f}s")
        else:
            print("  Skipping Python run (cpp-only mode).")

        if not args.python_only:
            cpp_log = log_dir / f"{workload.name}_cpp.log"
            cpp_time = run_command(workload.cpp_cmd, ROOT, cpp_log)
            workload_entry["cpp_time"] = cpp_time
            print(f"  C++ time: {cpp_time:.2f}s")
        else:
            print("  Skipping C++ run (python-only mode).")

        timing_summary.append(workload_entry)

        if not args.python_only and not args.cpp_only:
            comp_results = []
            for comp in workload.comparisons:
                py_path = Path(comp["python"])
                cpp_path = Path(comp["cpp"])
                if not py_path.exists():
                    raise FileNotFoundError(f"Missing Python output: {py_path}")
                if not cpp_path.exists():
                    raise FileNotFoundError(f"Missing C++ output: {cpp_path}")
                summary = run_compare(py_path, cpp_path, args.tolerance, compare_dir / workload.name)
                summary["file"] = comp["label"]
                comp_results.append(summary)
            comparison_summary[workload.name] = comp_results

        completed.add(workload.name)

    # Write summary
    summary = {
        "workloads": timing_summary,
        "comparisons": comparison_summary,
        "tolerance": args.tolerance,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with (work_dir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print("\n=== Summary ===")
    for entry in timing_summary:
        name = entry["name"]
        py_time = entry["python_time"]
        cpp_time = entry["cpp_time"]
        print(f"{name:>12}: python={py_time!s:>8}  cpp={cpp_time!s:>8}")
        if name in comparison_summary:
            for comp in comparison_summary[name]:
                status = "PASS" if comp["max_abs_diff"] <= args.tolerance else "FAIL"
                print(f"    {comp['file']}: max_abs_diff={comp['max_abs_diff']:.3e} [{status}]")

    print(f"\nDetailed artifacts stored in: {work_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())


