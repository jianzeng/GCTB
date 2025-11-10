"""
High-level workflows built on top of the GCTB Python bindings.

These helpers orchestrate multi-chain MCMC runs and other composite
tasks that previously lived in the C++ `GCTB` controller.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

try:
    import scipy.sparse as sp
except ImportError:  # pragma: no cover - optional dependency
    sp = None

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

