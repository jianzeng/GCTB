"""
High-level workflows built on top of the GCTB Python bindings.

These helpers orchestrate multi-chain MCMC runs and other composite
tasks that previously lived in the C++ `GCTB` controller.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Sequence, Tuple

import numpy as np

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
except ImportError:  # pragma: no cover - optional dependency
    sp = None
    spla = None

from . import _core as gctb


def _default_output_freq(chain_length: int) -> int:
    return max(chain_length // 10, 100)


def run_multi_chain_sbayes(
    data: gctb.Data,
    sbayes_type: str,
    *,
    heritability: float,
    pi: float,
    num_chains: int = 4,
    chain_length: int = 3000,
    burnin: int = 1000,
    thin: int = 10,
    output_prefix: str = "sbayes",
    random_start: bool = True,
    verbose: bool = True,
) -> List[List[gctb.McmcSamples]]:
    """
    Run multiple SBayes chains and compute Gelman–Rubin diagnostics.

    Parameters
    ----------
    data : gctb.Data
        Prepared data object (LD matrix and summary statistics loaded).
    sbayes_type : str
        One of ``"C"``, ``"R"``, ``"S"`` (case insensitive).
    heritability : float
        Prior/initial heritability used to initialise model variances.
    pi : float
        Prior probability of non-zero SNP effects.
    num_chains : int, default 4
        Number of independent chains to run.
    chain_length : int, default 3000
        MCMC chain length for each chain.
    burnin : int, default 1000
        Burn-in iterations.
    thin : int, default 10
        Thinning interval.
    output_prefix : str, default ``"sbayes"``
        Prefix for per-chain output artefacts and diagnostics.
    random_start : bool, default True
        Whether to randomise starting values for chains after the first.
    verbose : bool, default True
        Emit model initialisation logs from the first chain.

    Returns
    -------
    list[list[gctb.McmcSamples]]
        A list of chains, each containing the usual collection of
        ``McmcSamples`` objects returned by ``MCMC.run``.
    """
    if num_chains < 2:
        raise ValueError("num_chains must be at least 2 for multi-chain workflow")

    chains: List[List[gctb.McmcSamples]] = []
    models: List[gctb.Model] = []

    sbayes_type = sbayes_type.upper()

    for chain_idx in range(num_chains):
        model = gctb.build_model_summary(
            data,
            sbayes_type,
            heritability,
            pi,
            random_start=bool(random_start and chain_idx > 0),
            verbose=bool(verbose and chain_idx == 0),
        )
        models.append(model)

        runner = gctb.MCMC()
        title = f"{output_prefix}_chain{chain_idx + 1}"
        samples = runner.run(
            model=model,
            chain_length=chain_length,
            burnin=burnin,
            thin=thin,
            output_freq=_default_output_freq(chain_length),
            title=title,
        )
        chains.append(samples)

    # Gelman–Rubin diagnostics (written to <output_prefix>.gelman)
    diag = gctb.MCMC()
    diag_title = str(Path(output_prefix))
    diag.gelman_rubin(models[0], chains, diag_title)

    return chains


def compute_window_pip(
    data: gctb.Data,
    *,
    snp_results: str,
    window_width: int = 1_000_000,
    step_size: int | None = None,
) -> List[Dict[str, Any]]:
    """
    Compute window-level probabilities of association using SNP PIPs.

    Parameters
    ----------
    data : gctb.Data
        Data object with LD metadata already initialised.
    snp_results : str
        Path to a ``*.snpRes`` file containing SNP posterior inclusion probabilities.
    window_width : int, default 1_000_000
        Physical width (bp) for constructing windows.
    step_size : int | None, default None
        Optional stride (bp) for overlapping windows.  When omitted the windows
        are non-overlapping.

    Returns
    -------
    list of dict
        Per-window summaries including window-level PIP and per-SNP details.
    """
    data.input_new_snp_results(snp_results)
    if step_size is None:
        data.get_nonoverlap_window_info(window_width)
    else:
        data.get_overlap_windows(window_width, step_size)

    snp_info = data.get_snp_info_vec()
    window_starts = data.window_starts
    window_sizes = data.window_sizes

    windows: List[Dict[str, Any]] = []

    for idx, (start, size) in enumerate(zip(window_starts, window_sizes)):
        indices = range(start, start + size)
        snps = [
            {
                "id": snp_info[i].ID,
                "pip": float(snp_info[i].pip),
                "chrom": int(snp_info[i].chrom),
                "position": int(snp_info[i].physPos),
            }
            for i in indices
        ]
        pip_product = float(np.prod([max(1.0 - snp["pip"], 1e-12) for snp in snps])) if snps else 1.0
        windows.append(
            {
                "index": idx + 1,
                "size": size,
                "start": snps[0]["id"] if snps else None,
                "end": snps[-1]["id"] if snps else None,
                "window_pip": 1.0 - pip_product,
                "snps": snps,
            }
        )

    return windows


def compute_window_credible_sets(
    data: gctb.Data,
    *,
    snp_results: str,
    window_width: int = 1_000_000,
    pip_threshold: float = 0.9,
) -> Dict[str, Any]:
    """
    Build non-overlapping genomic windows and derive credible sets using SNP PIPs.

    This mirrors the behaviour of the legacy ``GCTB::calcCredibleSets`` routine
    but operates entirely from Python, relying on the SNP PIP values loaded
    from ``snp_results``.  The per-window probability of association is
    approximated as ``1 - prod(1 - pip_i)`` under independence.

    Parameters
    ----------
    data : gctb.Data
        Data object with LD metadata already initialised.
    snp_results : str
        Path to a ``*.snpRes`` file produced by ``data.outputSnpResults`` or
        the Python equivalents.
    window_width : int, default 1_000_000
        Physical window width (bp) used to create non-overlapping windows.
    pip_threshold : float, default 0.9
        Cumulative PIP threshold used to define credible sets.

    Returns
    -------
    dict
        Dictionary containing per-window summaries and genome-wide credible sets.
    """
    windows = compute_window_pip(
        data,
        snp_results=snp_results,
        window_width=window_width,
        step_size=None,
    )

    genomewide: List[Dict[str, Any]] = []
    enriched_windows: List[Dict[str, Any]] = []

    for window in windows:
        snps = window["snps"]
        credible: List[Dict[str, Any]] = []
        cumulative = 0.0

        for snp in sorted(snps, key=lambda s: s["pip"], reverse=True):
            credible.append(snp)
            cumulative += snp["pip"]
            if cumulative >= pip_threshold:
                break

        enriched = dict(window)
        enriched["credible_set"] = credible
        genomewide.extend(credible)
        enriched_windows.append(enriched)

    return {
        "windows": enriched_windows,
        "genomewide_credible_set": genomewide,
        "window_width": window_width,
        "pip_threshold": pip_threshold,
    }


def pip_to_pvalues(pips: Sequence[float], prop_null: float) -> List[float]:
    """
    Convert posterior inclusion probabilities to p-value style metrics.

    Mirrors the legacy ``GCTB::pip2p`` helper which assumes a uniform prior
    on causal variants and estimates the expected number of null SNPs.
    """
    if not 0.0 < prop_null <= 1.0:
        raise ValueError("prop_null must be within (0, 1]")

    pip_array = np.asarray(pips, dtype=float)
    if pip_array.ndim != 1:
        raise ValueError("pips must be a one-dimensional sequence")

    sorted_idx = np.argsort(pip_array)[::-1]
    sorted_pip = pip_array[sorted_idx]

    cumsum = np.cumsum(sorted_pip)
    num_null = len(pip_array) * prop_null
    p_values = (np.arange(1, len(sorted_pip) + 1) - cumsum) / num_null
    p_values = np.clip(p_values, 0.0, None)

    restored = np.empty_like(p_values)
    restored[sorted_idx] = p_values
    return restored.tolist()


def mcmc_samples_sparse_matrix(samples: gctb.McmcSamples) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Tuple[int, int]]:
    """
    Retrieve sparse COO data for an MCMC samples object.

    Returns
    -------
    rows : ndarray[int]
    cols : ndarray[int]
    data : ndarray[float]
    shape : tuple[int, int]
    """
    if samples.storage_mode != "sparse":
        raise ValueError("Expected sparse storage mode for MCMC samples")

    rows, cols, data, shape = samples.sparse_data()
    return np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32), np.asarray(data, dtype=np.float32), tuple(shape)


def mcmc_samples_to_csr(samples: gctb.McmcSamples):
    """
    Convert sparse MCMC samples to a SciPy CSR matrix if SciPy is available.
    """
    if sp is None:
        raise RuntimeError("scipy is required to materialise sparse MCMC samples; install scipy to use this helper.")

    rows, cols, data, shape = mcmc_samples_sparse_matrix(samples)
    return sp.csr_matrix((data, (rows, cols)), shape=shape)


# ---------------------------------------------------------------------------
# Credible set helpers
# ---------------------------------------------------------------------------

_THRESHOLD_GRID = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], dtype=float)


def _sum_columns(matrix: Any, indices: Sequence[int]) -> np.ndarray:
    if not indices:
        return np.zeros(matrix.shape[0], dtype=float)
    if isinstance(matrix, np.ndarray):
        return matrix[:, indices].sum(axis=1).astype(float, copy=False)
    sub = matrix[:, indices]
    return np.asarray(sub.sum(axis=1)).ravel().astype(float, copy=False)


def _column_vector(matrix: Any, index: int) -> np.ndarray:
    if isinstance(matrix, np.ndarray):
        return matrix[:, index].astype(float, copy=True)
    col = matrix.getcol(index)
    return np.asarray(col.toarray()).ravel().astype(float, copy=False)


def _normalised_mcmc_matrix(samples: gctb.McmcSamples) -> tuple[Any, Any, np.ndarray]:
    if samples.storage_mode == "dense":
        mat = samples.dense_matrix().astype(float, copy=True)
        betas_sq = mat ** 2
        total_var = betas_sq.sum(axis=1)
        normalised = np.zeros_like(betas_sq, dtype=float)
        mask = total_var > 0
        if mask.any():
            normalised[mask] = betas_sq[mask] / total_var[mask, None]
        return normalised, betas_sq, total_var

    if samples.storage_mode == "sparse":
        csr = mcmc_samples_to_csr(samples)
        csr_sq = csr.copy()
        csr_sq.data = csr_sq.data.astype(float) ** 2
        total_var = np.asarray(csr_sq.sum(axis=1)).ravel()
        inv_total = np.zeros_like(total_var)
        mask = total_var > 0
        if mask.any():
            inv_total[mask] = 1.0 / total_var[mask]
        diag = sp.diags(inv_total) if mask.any() else sp.csr_matrix(np.zeros((csr_sq.shape[0], csr_sq.shape[0])))
        normalised = diag.dot(csr_sq)
        return normalised, csr_sq, total_var

    raise ValueError(f"Unsupported storage mode {samples.storage_mode!r}")


def compute_credible_sets(
    data: Any,
    mcmc_samples: Any,
    *,
    pip_threshold: float = 0.9,
    pep_threshold: float = 0.9,
    window_width: int | None = None,
    snp_results_path: str | None = None,
    unconverged_snplist: str | None = None,
) -> Dict[str, Any]:
    """
    Python reimplementation of ``GCTB::calcCredibleSets``.
    """
    if not (0 < pip_threshold <= 1):
        raise ValueError("pip_threshold must be in (0, 1]")
    if not (0 < pep_threshold <= 1):
        raise ValueError("pep_threshold must be in (0, 1]")

    if window_width is not None and hasattr(data, "get_nonoverlap_window_info"):
        data.get_nonoverlap_window_info(int(window_width))

    if snp_results_path and hasattr(data, "input_new_snp_results"):
        data.input_new_snp_results(snp_results_path)

    if unconverged_snplist and hasattr(data, "read_unconverged_snplist"):
        data.read_unconverged_snplist(unconverged_snplist)

    window_starts = list(getattr(data, "window_starts", []))
    window_sizes = list(getattr(data, "window_sizes", []))
    num_windows = len(window_starts)
    if num_windows == 0:
        raise ValueError("No windows defined. Provide window_width or pre-computed windows.")

    snp_info_vec = list(data.get_snp_info_vec())
    num_snps = getattr(mcmc_samples, "ncol", len(mcmc_samples.posterior_mean))

    posterior_mean = np.asarray(mcmc_samples.posterior_mean, dtype=float)
    posterior_sqr_mean = np.asarray(mcmc_samples.posterior_sqr_mean, dtype=float)
    pip = np.asarray(mcmc_samples.pip, dtype=float) if getattr(mcmc_samples, "pip", None) is not None else np.zeros(num_snps, dtype=float)

    normalised, squared_matrix, total_var = _normalised_mcmc_matrix(mcmc_samples)
    n_iterations = normalised.shape[0]
    if n_iterations == 0:
        raise ValueError("MCMC samples contain no iterations.")

    if isinstance(normalised, np.ndarray):
        var_explained = normalised.sum(axis=0) / n_iterations
    else:
        var_explained = np.asarray(normalised.sum(axis=0)).ravel() / n_iterations

    records_by_index: Dict[int, Dict[str, Any]] = {}
    snp_records: List[Dict[str, Any]] = []
    improper_components = np.zeros(num_windows, dtype=float)

    for snp in snp_info_vec:
        idx = int(getattr(snp, "index", 0)) - 1
        if idx < 0 or idx >= num_snps:
            continue
        effect = float(posterior_mean[idx])
        snp.effect = effect
        snp.pip = float(pip[idx]) if pip.size else float(getattr(snp, "pip", 0.0))
        snp.varExplained = float(var_explained[idx])

        record = {
            "snp": snp,
            "index": idx,
            "id": str(snp.ID),
            "chrom": int(snp.chrom),
            "position": int(snp.physPos),
            "pip": float(pip[idx]) if pip.size else float(getattr(snp, "pip", 0.0)),
            "effect": effect,
            "var_explained": float(var_explained[idx]),
            "af": float(getattr(snp, "af", 0.0)),
            "a1": str(getattr(snp, "a1", "")),
            "a2": str(getattr(snp, "a2", "")),
            "unconverged": bool(getattr(snp, "unconverged", False)),
        }
        records_by_index[idx] = record
        snp_records.append(record)

    if not snp_records:
        raise ValueError("No SNP records available for credible set calculation.")

    # Window initial metrics
    window_entries: List[Dict[str, Any]] = []
    for w_idx, (start, size) in enumerate(zip(window_starts, window_sizes), start=1):
        window_snps: List[Dict[str, Any]] = []
        column_indices: List[int] = []
        effect_improper = 0.0
        for offset in range(size):
            snp = snp_info_vec[start + offset]
            idx = int(getattr(snp, "index", 0)) - 1
            if idx < 0 or idx >= num_snps:
                continue
            record = records_by_index.get(idx)
            if record is None:
                continue
            window_snps.append(record)
            column_indices.append(idx)
            af = record["af"]
            effect = record["effect"]
            effect_improper += 2.0 * af * (1.0 - af) * effect * effect

        if not column_indices:
            continue

        initial_vec = _sum_columns(normalised, column_indices)
        window_entries.append(
            {
                "index": w_idx,
                "size": len(column_indices),
                "column_indices": column_indices,
                "snps": window_snps,
                "initial_vec": initial_vec.copy(),
                "residual_vec": initial_vec.copy(),
                "effect_improper": effect_improper,
            }
        )

    if not window_entries:
        raise ValueError("Unable to build windows with current configuration.")

    num_windows_effective = len(window_entries)
    average_window_share = 1.0 / num_windows_effective

    for entry in window_entries:
        initial_vec = entry["initial_vec"]
        entry["initial_prop_gen_var"] = float(initial_vec.mean())
        entry["initial_gen_var_enrich"] = float(initial_vec.mean() * num_windows_effective)
        entry["initial_gen_var_enrich_pp"] = float(np.mean(initial_vec > average_window_share))

    vg = sum(entry["effect_improper"] for entry in window_entries)
    if vg > 0:
        for entry in window_entries:
            entry["pve_improper"] = entry["effect_improper"] / vg
    else:
        for entry in window_entries:
            entry["pve_improper"] = 0.0

    credible_sets: List[Dict[str, Any]] = []
    num_cs = 0
    num_single = 0
    sum_cs_size = 0
    cumsum_pip = 0.0
    cumsum_prop_var = 0.0

    for entry in window_entries:
        window_vec = entry["residual_vec"]
        sum_cs_snp_window = 0
        cumulative_pip = 0.0
        cumulative_prop_var = 0.0
        top_snps: List[Dict[str, Any]] = []

        for record in sorted(entry["snps"], key=lambda r: r["pip"], reverse=True):
            if record["unconverged"]:
                continue

            varj = _column_vector(normalised, record["index"])
            window_vec -= varj
            np.maximum(window_vec, 0.0, out=window_vec)

            prop_gen_var = float(window_vec.mean())
            gen_var_enrich = float(prop_gen_var * num_windows_effective)
            gen_var_enrich_pp = float(np.mean(window_vec > average_window_share))

            if record["pip"] > pip_threshold:
                num_cs += 1
                num_single += 1
                sum_cs_size += 1
                cumsum_pip += record["pip"]
                cumsum_prop_var += record["var_explained"]
                credible_sets.append(
                    {
                        "index": num_cs,
                        "size": 1,
                        "sum_pip": record["pip"],
                        "prop_var": record["var_explained"],
                        "wind_size": entry["size"] - sum_cs_snp_window,
                        "wind_prop_gen_var": prop_gen_var,
                        "wind_gen_var_enrich": gen_var_enrich,
                        "wind_gen_var_enrich_pp": gen_var_enrich_pp,
                        "snps": [
                            {
                                "id": record["id"],
                                "pip": record["pip"],
                                "var_explained": record["var_explained"],
                                "effect": record["effect"],
                                "chrom": record["chrom"],
                                "position": record["position"],
                            }
                        ],
                    }
                )
                sum_cs_snp_window += 1
                continue

            if gen_var_enrich_pp > pep_threshold:
                cumulative_pip += record["pip"]
                cumulative_prop_var += record["var_explained"]
                top_snps.append(record)
                if cumulative_pip > pip_threshold:
                    num_cs += 1
                    sum_cs_size += len(top_snps)
                    cumsum_pip += cumulative_pip
                    cumsum_prop_var += cumulative_prop_var
                    credible_sets.append(
                        {
                            "index": num_cs,
                            "size": len(top_snps),
                            "sum_pip": cumulative_pip,
                            "prop_var": cumulative_prop_var,
                            "wind_size": entry["size"] - sum_cs_snp_window,
                            "wind_prop_gen_var": prop_gen_var,
                            "wind_gen_var_enrich": gen_var_enrich,
                            "wind_gen_var_enrich_pp": gen_var_enrich_pp,
                            "snps": [
                                {
                                    "id": snp_rec["id"],
                                    "pip": snp_rec["pip"],
                                    "var_explained": snp_rec["var_explained"],
                                    "effect": snp_rec["effect"],
                                    "chrom": snp_rec["chrom"],
                                    "position": snp_rec["position"],
                                }
                                for snp_rec in top_snps
                            ],
                        }
                    )
                    top_snps = []
                    cumulative_pip = 0.0
                    cumulative_prop_var = 0.0
                    break
            else:
                break

    pip_values = np.array([rec["pip"] for rec in snp_records], dtype=float)
    nnz = float(len(pip_values)) * float(pip_values.mean()) if pip_values.size else 0.0
    estimated_power = (cumsum_pip / nnz) if nnz else 0.0
    average_cs_size = (sum_cs_size / num_cs) if num_cs else 0.0

    sorted_snps = sorted(snp_records, key=lambda r: r["pip"], reverse=True)
    target = pip_threshold * nnz
    cumulative = 0.0
    genomewide_top: List[Dict[str, Any]] = []
    for rec in sorted_snps:
        genomewide_top.append(
            {
                "id": rec["id"],
                "pip": rec["pip"],
                "var_explained": rec["var_explained"],
            }
        )
        cumulative += rec["pip"]
        if cumulative > target:
            break

    threshold_curve: List[Dict[str, Any]] = []
    for threshold in _THRESHOLD_GRID:
        limit = threshold * nnz
        cum = 0.0
        prop_var = 0.0
        count = 0
        for rec in sorted_snps:
            cum += rec["pip"]
            prop_var += rec["var_explained"]
            count += 1
            if cum > limit:
                break
        threshold_curve.append(
            {
                "threshold": threshold,
                "cs_size": count,
                "prop_hsq": prop_var,
            }
        )

    return {
        "credible_sets": credible_sets,
        "summary": {
            "pip_threshold": pip_threshold,
            "pep_threshold": pep_threshold,
            "num_cs": num_cs,
            "num_single": num_single,
            "num_multi": num_cs - num_single,
            "sum_cs_size": sum_cs_size,
            "nnz": nnz,
            "estimated_power": estimated_power,
            "estimated_identified": cumsum_pip,
            "estimated_prop_var": cumsum_prop_var,
            "average_cs_size": average_cs_size,
        },
        "global_top_snps": genomewide_top,
        "threshold_curve": threshold_curve,
        "window_metrics": [
            {
                "index": entry["index"],
                "size": entry["size"],
                "prop_gen_var": entry["initial_prop_gen_var"],
                "gen_var_enrich": entry["initial_gen_var_enrich"],
                "gen_var_enrich_pp": entry["initial_gen_var_enrich_pp"],
                "pve_improper": entry["pve_improper"],
            }
            for entry in window_entries
        ],
        "pip_threshold": pip_threshold,
        "pep_threshold": pep_threshold,
        "window_width": window_width,
        "num_windows": num_windows_effective,
    }


def compute_window_pip_trace(
    data: gctb.Data,
    snp_effects: gctb.McmcSamples,
    *,
    window_width: int,
    step_size: int,
    snp_results: str,
) -> Dict[str, Any]:
    """
    Python analogue of ``GCTB::getWindowPIP``.
    """
    data.input_new_snp_results(snp_results)
    data.get_overlap_windows(window_width, step_size)

    num_windows = len(data.window_starts)
    if snp_effects.storage_mode != "sparse":
        raise ValueError("Expected sparse storage mode for window PIP trace")

    rows, cols, values, shape = mcmc_samples_sparse_matrix(snp_effects)
    num_iters = shape[0]
    window_delta = np.zeros((num_iters, num_windows), dtype=float)
    snp_pip_counts = np.zeros(shape[1], dtype=float)

    window_assignments = [getattr(data.get_snp_info_vec()[i], "window", -1) for i in range(shape[1])]
    iter_offsets = rows
    snp_indices = cols

    for iter_idx, snp_idx, value in zip(iter_offsets, snp_indices, values):
        win_idx = window_assignments[snp_idx]
        if 0 <= win_idx < num_windows:
            window_delta[iter_idx, win_idx] = 1.0
        snp_pip_counts[snp_idx] += 1.0

    window_pip = window_delta.mean(axis=0)
    snp_pip = snp_pip_counts / snp_effects.nrow if snp_effects.nrow else np.zeros_like(snp_pip_counts)

    credible_windows: Dict[int, List[Dict[str, Any]]] = {}
    snp_info = data.get_snp_info_vec()
    for win_idx, pip in enumerate(window_pip):
        if pip <= 0.9:
            continue
        start = data.window_starts[win_idx]
        size = data.window_sizes[win_idx]
        indices = range(start, start + size)
        snp_pip_window = sorted(
            [(i, snp_pip[i]) for i in indices],
            key=lambda t: t[1],
            reverse=True,
        )
        cumulative = 0.0
        credible: List[Dict[str, Any]] = []
        for snp_idx, pip_val in snp_pip_window:
            snp = snp_info[snp_idx]
            credible.append({"id": snp.ID, "pip": float(pip_val)})
            cumulative += pip_val
            if cumulative > 0.9:
                break
        credible_windows[win_idx + 1] = credible

    return {
        "window_delta": window_delta,
        "window_pip": window_pip,
        "snp_pip": snp_pip,
        "credible_windows": credible_windows,
    }


def solve_snp_effects_cg(
    data: gctb.Data,
    *,
    lambda_: float,
    output_path: str,
) -> Dict[str, Any]:
    """
    Python implementation of ``GCTB::solveSnpEffectsByConjugateGradientMethod``.
    """
    if sp is None or spla is None:
        raise RuntimeError("scipy is required for conjugate gradient solver; install scipy to use this helper.")

    rows, cols, values, shape = data.get_zpz_sparse_matrix()
    mat = sp.csr_matrix((np.asarray(values, dtype=np.float64), (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))), shape=shape, dtype=np.float64)
    if lambda_ != 0.0:
        mat = mat + lambda_ * sp.identity(mat.shape[0], format="csr", dtype=np.float64)

    zpy = np.asarray(data.zpy, dtype=np.float64)
    if zpy.shape[0] != mat.shape[0]:
        raise ValueError(f"ZPy length {zpy.shape[0]} does not match matrix dimension {mat.shape[0]}")

    sol, info = spla.cg(mat, zpy, atol=0)
    if info != 0:
        raise RuntimeError(f"Conjugate gradient did not converge (info={info})")

    snps = data.get_incd_snp_info_vec()
    if len(snps) != len(sol):
        raise ValueError("Number of SNPs does not match solution length")

    results = []
    for idx, snp in enumerate(snps):
        af = float(getattr(snp, "af", 0.0))
        sqrt2pq = np.sqrt(max(2.0 * af * (1.0 - af), 1e-12))
        effect = float(sol[idx])
        flipped = bool(getattr(snp, "flipped", False))
        allele1 = getattr(snp, "a1", "A")
        allele2 = getattr(snp, "a2", "G")
        if flipped:
            allele1, allele2 = allele2, allele1
            freq = 1.0 - af
            effect = -effect
        else:
            freq = af
        results.append(
            {
                "index": idx + 1,
                "id": getattr(snp, "ID", f"snp_{idx+1}"),
                "chrom": int(getattr(snp, "chrom", 0)),
                "position": int(getattr(snp, "physPos", 0)),
                "a1": allele1,
                "a2": allele2,
                "freq": freq,
                "effect": effect,
                "scaled_effect": effect / sqrt2pq if sqrt2pq > 0 else effect,
            }
        )

    header = "Index\tName\tChrom\tPosition\tA1\tA2\tA1Frq\tA1Sol\n"
    with open(output_path, "w") as fout:
        fout.write(header)
        for rec in results:
            fout.write(
                f"{rec['index']}\t{rec['id']}\t{rec['chrom']}\t{rec['position']}\t"
                f"{rec['a1']}\t{rec['a2']}\t{rec['freq']:.6f}\t{rec['effect']:.6f}\n"
            )

    return {
        "lambda": lambda_,
        "num_snps": len(sol),
        "matrix_nnz": int(mat.nnz),
        "output": output_path,
    }


def run_posthoc_stratify(
    ldm_prefix: str,
    snp_results: str,
    mcmc_prefix: str,
    bayes_type: str,
    *,
    annotation_file: str | None = None,
    annotation_transpose: bool = False,
    continuous_annotation: str | None = None,
    flank: int = 0,
    eqtl_file: str | None = None,
    genetic_map_file: str | None = None,
    gen_map_n: float = 60_000.0,
    gwas_summary: str | None = None,
    pvalue_threshold: float = 1.0,
    impute_n: bool = True,
    multi_ldm: bool = False,
    chain_length: int = 1500,
    burnin: int = 500,
    thin: int = 5,
    output_prefix: str = "posthoc_stratify",
) -> List[gctb.McmcSamples]:
    """
    Reimplementation of the post-hoc stratified workflow from ``GCTB::stratify``.
    """
    bayes_type = bayes_type.upper()
    data = gctb.Data()

    # Load LD metadata
    if multi_ldm:
        data.read_multi_ld_matrix_info_file(ldm_prefix)
    else:
        info_path = f"{ldm_prefix}.info"
        try:
            data.read_ld_matrix_info_file(info_path)
        except Exception:
            data.read_ld_matrix_info_file(ldm_prefix)

    data.input_snp_info_and_results(snp_results, bayes_type)

    if annotation_file:
        data.read_annotation_file(annotation_file, annotation_transpose, True)
    elif continuous_annotation:
        data.read_annotation_file_format2(continuous_annotation, flank * 1000, eqtl_file or "")

    if gwas_summary:
        data.read_gwas_summary_file(
            gwas_file=gwas_summary,
            af_diff=1.0,
            maf_min=0.0,
            maf_max=0.0,
            pvalue_threshold=pvalue_threshold,
            impute_n=impute_n,
            remove_outlier_n=True,
        )

    data.include_matched_snp()

    if genetic_map_file:
        if multi_ldm:
            data.read_multi_ld_matrix_bin_file_and_shrink(ldm_prefix, gen_map_n)
        else:
            bin_path = f"{ldm_prefix}.bin"
            try:
                data.read_ld_matrix_bin_file_and_shrink(bin_path)
            except Exception:
                data.read_ld_matrix_bin_file_and_shrink(ldm_prefix)
    else:
        if multi_ldm:
            data.read_multi_ld_matrix_bin_file(ldm_prefix)
        else:
            bin_path = f"{ldm_prefix}.bin"
            try:
                data.read_ld_matrix_bin_file(bin_path)
            except Exception:
                data.read_ld_matrix_bin_file(ldm_prefix)

    data.build_sparse_mme(sample_overlap=False, noscale=True)
    data.make_annowise_sparse_ldm()

    # Load MCMC samples
    snp_effects = gctb.McmcSamples("SnpEffects")
    snp_effects.read_data_bin(mcmc_prefix)

    hsq = gctb.McmcSamples("hsq")
    hsq.read_data_txt(f"{mcmc_prefix}.Par", "hsq")
    hsq_hat_vec = hsq.mean()
    hsq_hat = float(hsq_hat_vec[0]) if hsq_hat_vec.size > 0 else 0.0

    delta_s = None
    if bayes_type == "SMIX":
        delta_s = gctb.McmcSamples("DeltaS")
        delta_s.read_data_bin(mcmc_prefix)

    model = gctb.build_posthoc_model(
        data,
        bayes_type,
        snp_effects,
        hsq,
        thin=thin,
        hsq_hat=hsq_hat,
        delta_s=delta_s,
    )

    runner = gctb.MCMC()
    samples = runner.run(
        model=model,
        chain_length=chain_length,
        burnin=burnin,
        thin=thin,
        output_freq=_default_output_freq(chain_length),
        title=output_prefix,
    )

    return {
        "samples": samples,
        "data": data,
        "snp_effects": snp_effects,
        "hsq": hsq,
        "delta_s": delta_s,
    }


def tune_eigen_cutoff(
    data: gctb.Data,
    eigen_prefix: str,
    cutoffs: Sequence[float],
    *,
    heritability: float,
    prop_var_random: float,
    pi: float,
    chain_length: int = 150,
    burnin: int = 100,
    thin: int = 1,
    noscale: bool = False,
    make_pseudo_summary: bool = False,
    runner_factory: Callable[[Any], Any] | None = None,
    model_factory: Callable[[gctb.Data, float], Any] | None = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Python reimplementation of ``GCTB::tuneEigenCutoff``.
    """
    cutoffs = [float(c) for c in cutoffs]
    if not cutoffs:
        raise ValueError("At least one eigen cutoff must be supplied")

    original_ngwas = np.asarray(data.n_gwas_block, dtype=float).copy()
    pseudo_ngwas = np.asarray(data.pseudo_gwas_ntrn_block, dtype=float)

    correlations: List[float] = []
    records: List[Dict[str, Any]] = []
    warning: str | None = None
    base_corr: float | None = None

    data.n_gwas_block = pseudo_ngwas

    try:
        for cutoff in cutoffs:
            if verbose:
                print(f"[tune_eigen_cutoff] Evaluating cutoff {cutoff}")
            data.read_eigen_matrix_binary_file_and_make_wq(
                eigen_prefix,
                cutoff,
                noscale=noscale,
                make_pseudo_summary=make_pseudo_summary,
            )
            data.init_variances(heritability, prop_var_random)
            if model_factory:
                model = model_factory(data, cutoff)
            else:
                model = gctb.build_model_summary(
                    data,
                    "R",
                    heritability=heritability,
                    pi=pi,
                )
            runner = runner_factory(model) if runner_factory else gctb.MCMC()
            samples = runner.run(
                model=model,
                chain_length=chain_length,
                burnin=burnin,
                thin=thin,
                output_freq=_default_output_freq(chain_length),
                title=f"{eigen_prefix}_cutoff_{cutoff:.3f}",
            )
            snp_effects = next((res for res in samples if res.label == "SnpEffects"), None)
            if snp_effects is None:
                raise RuntimeError("SnpEffects samples not produced during eigen cutoff tuning")

            beta_mean = np.asarray(snp_effects.posterior_mean, dtype=float)
            b_val = np.asarray(data.b_val, dtype=float)
            var_pheno = float(data.var_phenotypic)
            denom = np.sqrt(np.dot(beta_mean, beta_mean) * var_pheno)
            corr = float(beta_mean.dot(b_val) / denom) if denom else 0.0
            if base_corr is None:
                base_corr = corr if corr != 0 else 1.0
            rel = float(corr / base_corr) if base_corr else 0.0
            correlations.append(corr)
            records.append(
                {
                    "cutoff": cutoff,
                    "correlation": corr,
                    "relative": rel,
                }
            )
    finally:
        data.n_gwas_block = original_ngwas

    best_idx = int(np.argmax(correlations)) if correlations else -1
    best_cutoff = cutoffs[best_idx] if best_idx >= 0 else None

    if best_cutoff is not None and np.isclose(best_cutoff, min(cutoffs)):
        warning = (
            "Best eigen cutoff equals the minimum candidate; consider extending the range lower."
        )

    return {
        "records": records,
        "correlations": correlations,
        "best_cutoff": best_cutoff,
        "warning": warning,
    }

