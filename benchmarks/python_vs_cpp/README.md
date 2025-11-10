# Python vs C++ Benchmark Suite

This suite runs the core GCTB workflows with both the Python extension and
the original C++ executable, then compares their outputs and wall-clock
runtime.  It is intended to provide confidence that the ported workflows
behave identically and perform within the expected range.

## Prerequisites

1. Build the Python extension (editable install):
   ```bash
   pip install -e .
   ```
2. Build the legacy C++ CLI (if not already built):
   ```bash
   make -C scr
   ```
3. Ensure test data are available under `test/data/` (they are tracked in the repo).

## Quick Start

```bash
python benchmarks/python_vs_cpp/run_benchmarks.py \
    --work-dir /tmp/gctb_benchmarks
```

This will:

1. Run the BayesC and SBayesR workflows via the Python CLI.
2. Run the same workflows via the legacy `scr/gctb` executable.
3. Compare the resulting `.parRes` and `.snpRes` files with a numerical tolerance.
4. Record timing information and save a JSON summary to the work directory.

All logs, outputs, and comparisons are stored under `--work-dir` (default
`benchmarks/python_vs_cpp/results` inside the repository).

## Options

```
usage: run_benchmarks.py [-h] [--work-dir PATH] [--dataset-dir PATH]
                        [--python-only] [--cpp-only] [--tolerance TOL]
                        [--include-posthoc] [--annotation-file PATH]
```

- `--python-only` / `--cpp-only`: run a single side for debugging.
- `--include-posthoc`: include the post-hoc stratification workflow (requires
  `--annotation-file` pointing to a valid annotation file).
- `--tolerance`: absolute numeric tolerance used by `compare_results.py`.
- `--dataset-dir`: override the default `test/data` location.

## Output

Each run produces:

- `python/` and `cpp/` directories with raw output files.
- `logs/` containing stdout/stderr and per-command timing JSON.
- `summary.json` consolidating runtimes and comparison results.
- `comparisons/*.json` for individual file comparisons.

## Extending the Suite

- To add an additional workload (e.g., BayesB or multi-chain runs), edit
  `run_benchmarks.py` and extend the `WORKLOADS` list.
- To compare additional artefacts (e.g., window-level summaries), update both
  the workload definitions and `compare_results.py`.

## Known Limitations

- Post-hoc stratification is optional because it requires an annotation file
  that is not currently bundled with the repo.
- The suite focuses on correctness and wall-clock time for small reference
  datasets.  For full-scale benchmarks, adjust the chain lengths and data
  paths accordingly.


