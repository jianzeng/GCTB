# GCTB Python Extension Reference

## Overview

The Python package under `python/gctb` exposes the high-performance C++ implementation of GCTB through a pybind11 module named `_core`. The bindings live in `python/bindings/bindings.cpp` and wrap classes from the C++ core (`scr/*.cpp`, mirrored under `src/`). On top of the extension module the package provides pure-Python orchestration layers (`workflows.py`, `diagnostics.py`, `cli.py`) that offer higher-level workflows, diagnostics, and a command-line interface.

This document enumerates every function, method, and property made available through the Python bindings and their associated pure-Python helpers, and explains how they interoperate with the C++ backend.

```26:1146:python/bindings/bindings.cpp
PYBIND11_MODULE(_core, m) {
    m.doc() = "GCTB: Genome-wide Complex Trait Bayesian Analysis - C++ Core";
    ...
```

## C++ Extension Module `_core`

### `SnpInfo`

Wraps the C++ `SnpInfo` struct defined in `scr/data.hpp`. Instances represent a SNP record and are primarily returned by `Data.get_snp_info_vec()` and related methods.

- Constructor arguments mirror the C++ struct fields: `idx`, `id`, `allele1`, `allele2`, `chr`, `gpos`, `ppos`.
- Read-only attributes: `ID`, `a1`, `a2`, `chrom`, `physPos`.
- Mutable attributes mirror mutable members in C++: `genPos`, `index`, `af`, `twopq`, `included`, `unconverged`, `effect`, `pip`, `varExplained`, `gwas_b`, `gwas_se`, `gwas_n`, `gwas_pvalue`.

### `IndInfo`

Python view of `IndInfo` from `scr/data.hpp`, representing individuals in PLINK files.

- Constructor arguments: `idx`, `fid`, `pid`, `dad`, `mom`, `sex`.
- Read-only attributes: `famID`, `indID`, `catID`, `fatherID`, `motherID`, `sex`.
- Mutable attributes: `index`, `kept`, `phenotype`.

### `Data`

Central data container (`Data` in `scr/data.cpp`). All methods call the identically named C++ member, with exception wrappers to translate errors into Python `RuntimeError`.

#### Construction & metadata

- `Data()` – constructs an empty container.
- Read-only counts updated by C++: `num_snps`, `num_inds`, `num_incd_snps`, `num_kept_inds`, `num_annos`, `num_ld_blocks`.
- Read-only variances: `var_phenotypic`, `var_genotypic`, `var_residual`.
- Read-write label: `title`.

#### PLINK, phenotype, and auxiliary readers

Each method maps directly to a `Data::read*` C++ method.

- `read_fam_file(path)`
- `read_bim_file(path)`
- `read_bed_file(noscale, path)`
- `read_phenotype_file(path, mphen)`
- `read_covariate_file(path)`
- `read_random_covariate_file(path)`
- `read_residual_diag_file(path)`
- `read_gwas_summary_file(path, af_diff, maf_min, maf_max, pvalue_threshold, impute_n, remove_outlier_n)`
- `read_genetic_map_file(path)`
- `read_ld_block_info_file(path)`
- `read_ldscore_file(path)`
- `read_plink_af_file(path)`
- `read_plink_ld_txt_file(path)`
- `read_plink_ld_bin_file(path)`

#### Summary statistics & LD metadata readers

- `read_ld_matrix_info_file(path)` / `read_ld_matrix_bin_file(path)` – for single-chromosome LD matrices.
- `read_multi_ld_matrix_info_file(path)` / `read_multi_ld_matrix_bin_file(path)` – directory-based multi-chromosome LD matrices.
- `read_ld_matrix_txt_file(path)` – text-form LD matrix.
- `read_ld_matrix_bin_file_and_shrink(path)` / `read_multi_ld_matrix_bin_file_and_shrink(path, gen_map_n)` – load LD and apply shrinkage.
- `read_eigen_matrix_binary_file_and_make_wq(dirname, eigen_cutoff, noscale=False, make_pseudo_summary=False)` – load eigen-decomposed LD matrices and construct W/Q structures; uses cached pseudo GWAS data on the C++ side.
- `read_window_file(path)` – load pre-defined LD windows.
- `read_unconverged_snplist(path)` – flag SNPs with convergence issues.
- `input_new_snp_results(path)` – hydrate SNP info with results saved to disk.
- `input_snp_info_and_results(path, bayes_type)` – combined metadata/results import for post-hoc stratification workflows.

#### Annotation ingestion

- `read_annotation_file(path, transpose=False, allow_multi_anno=True)` – categorical annotations.
- `read_annotation_file_format2(path, flank, eqtl_file="")` – continuous annotations.
- `set_annotation_info()` – finalise C++ annotation metadata.
- `make_annowise_sparse_ldm()` – build per-annotation sparse LD matrices from current state.

#### Individual and SNP selection

These wrap filtering utilities inside `Data`.

- `include_snp(path)` / `exclude_snp(path)`
- `include_chr(chrom)` / `include_block(block)`
- `exclude_region(path)` / `exclude_mhc()` / `exclude_ambiguous_snp()`
- `include_skeleton_snp(path)` – restrict to skeleton SNP list.
- `include_matched_snp()` – populate included SNP set after metadata load.
- `build_kept_individuals()` – copy all individuals as “kept”.
- `keep_matched_ind(path="", keep_ind_max=UINT_MAX)` – filter individuals using phenotype/covariate checks.
- `filter_snp_by_ld_rsq(rsq_threshold)` – prune by LD r².
- `bin_snp_by_ld_rsq(rsq_threshold, title)` – bucket SNPs by LD r².
- `bin_snp_by_window_id()` – assign SNPs to current windows for downstream use.

#### Variance initialisation

- `init_variances(heritability, prop_var_random)` – initialises variance components and caches them on the C++ object.

#### LD matrix construction & manipulation

- `build_sparse_mme(sample_overlap=False, noscale=False)` – instantiate sparse mixed-model equations from summary statistics; uses previously loaded LD matrix.
- `part_ld_matrix(part_param, out_filename, ldmat_type)` – partition LD matrix per config file.
- `make_ld_matrix(bed_file, ldmat_type, chisq_threshold, ld_threshold, window_width, snp_range, filename, write_ldm_txt=False)` – construct LD matrix from genotypes.
- `make_shrunk_ld_matrix(...)` – LD shrinkage build variant.
- `make_block_ld_matrix(bed_file, ldmat_type, block, filename, write_ldm_txt=False, ld_block_region_wind=0)`
- `resize_ld_matrix(ldmat_type, chisq_threshold, window_width, ld_threshold, effpop_ne, cutoff, gen_map_n)`
- `output_ld_matrix(ldmat_type, filename, write_ldm_txt=False)`
- `direct_prune_ld_matrix(ldm_file, out_ldmat_type, chisq_threshold, title, write_ldm_txt=False)`
- `jackknife_ld_matrix(ldm_file, out_ldmat_type, title, write_ldm_txt=False)`
- `merge_ldm_info(out_ldmat_type, dirname)` – merge per-chromosome `.info` files.

#### Output writers

- `output_snp_results(posterior_mean, posterior_sqr_mean, pip, noscale, filename)`
- `output_snp_results(posterior_mean, posterior_sqr_mean, last_sample, pip, noscale, filename)` – overload including final sample.
- `output_fixed_effects(effects, filename)`
- `output_random_effects(effects, filename)`
- `output_window_results(posterior_mean, filename)`

Each mirrors a `Data::output*` function emitting text files identical to the legacy CLI.

#### Matrix accessors & derived data

- `get_snp_info_vec()` / `get_incd_snp_info_vec()` / `get_ind_info_vec()` – return vectors of `SnpInfo`/`IndInfo` (reference semantics).
- `get_zpz_sparse_matrix()` – export sparse ZᵀZ (COO arrays + shape). Calls `Data::getZPZspmat()` first.
- `wind_size` (property) – returns window sizes per SNP (vector clone of `windSize`).
- `zpy` (property) – right-hand side vector `ZPy`.
- `n_gwas_block` (property with setter) – GWAS sample sizes per block (used during eigen tuning).
- `pseudo_gwas_ntrn_block`, `pseudo_gwas_effect_trn`, `b_val`, `var_phenotypic` – expose internal vectors used for pseudo-summary workflows.

#### Window utilities

- `get_overlap_windows(window_width, step_size)` – build overlapping windows.
- `get_nonoverlap_window_info(window_width)` – build disjoint windows.
- `num_windows` (property) – count of current windows.
- `window_starts` / `window_sizes` (properties) – parallel arrays describing each window.

### `Model` and factory helpers

`Model` objects are opaque wrappers around Bayesian models. Only `num_snps` is exposed for inspection.

- `build_model(data, bayes_type, heritability=0.5, pi=0.01)` – instantiates individual-level BayesC/B/R/S models. The lambda selects the correct C++ subclass (`BayesC`, `BayesB`, `BayesR`, `BayesS`) after calling `Data::initVariances`.
- `build_model_summary(data, sbayes_type, heritability=0.5, pi=0.01, random_start=False, verbose=True)` – constructs summary-statistics counterparts (`ApproxBayesC/R/S`) using `data.lowRankModel`.
- `build_posthoc_model(data, bayes_type, snp_effects, hsq, thin=1, hsq_hat=0.0, delta_s=None)` – wraps C++ post-hoc stratification models (`PostHocStratifyS`, `PostHocStratifySMix`).

### `McmcSamples`

Container for MCMC samples (`scr/mcmc.cpp`).

- Constructor: `McmcSamples(label)`.
- Read-only metadata: `label`, `chain_length`, `burnin`, `thin`, `nrow`, `ncol`, `posterior_mean`, `posterior_sqr_mean`, `pip`, `last_sample`.
- `storage_mode` – returns `"dense"` or `"sparse"` according to the underlying enum.
- `dense_matrix()` – materialise samples as a NumPy array (requires dense storage).
- `sparse_data()` – return COO components and shape, mirroring `datMatSp`.
- `read_data_bin(title)` – load sparse binary `.mcmcsamples.bin`.
- `read_data_txt(title)` – load dense text `.mcmcsamples.txt` (overload for combined files: `read_data_txt(filename, label)`).
- `mean()` – compute posterior mean across iterations.
- `to_dict()` – convenience dictionary with key statistics.

### `MCMC` and helpers

Wraps the sampler in `scr/mcmc.cpp`.

- `MCMC()` – constructor.
- `run(model, chain_length=3000, burnin=1000, thin=10, output_freq=100, title="gctb")` – executes sampling; releases the GIL inside the binding.
- `gelman_rubin(model, chains, title)` – compute Gelman–Rubin diagnostics across chains.
- Module-level `run_mcmc(model, ...)` – convenience function that instantiates `MCMC` internally.

### Utilities

- `Timer` mirrors `Gadget::Timer` from `scr/gadgets.hpp` with methods `set_time()`, `get_time()`, `get_date()`, `get_elapse()`, `format(seconds)`.
- Module attributes: `__version__ = "1.2.0"`, `__author__` list.

## Pure Python Package Layer (`python/gctb`)

### `__init__.py`

Imports the `_core` module and re-exports `SnpInfo`, `IndInfo`, `Data`, `Model`, `MCMC`, `McmcSamples`, `build_model`, `build_model_summary`, `run_mcmc`, `Timer`, along with workflow helpers (`compute_credible_sets`, `run_posthoc_stratify`, `solve_snp_effects_cg`, `tune_eigen_cutoff`). If SciPy is present it also exposes the `diagnostics` module. The top-level `__version__` mirrors the extension metadata.

### `workflows.py`

Higher-level orchestrators that compose `_core` functionality:

- `_default_output_freq(chain_length)` – internal helper matching C++ defaults.
- `run_multi_chain_sbayes(data, sbayes_type, *, heritability, pi, num_chains=4, chain_length=3000, burnin=1000, thin=10, output_prefix="sbayes", random_start=True, verbose=True)` – repeatedly calls `build_model_summary` and `MCMC.run`, then `MCMC.gelman_rubin` for diagnostics.
- `compute_window_pip(data, *, snp_results, window_width=1_000_000, step_size=None)` – loads SNP results via `Data.input_new_snp_results` and window metadata, then aggregates PIPs.
- `compute_window_credible_sets(data, *, snp_results, window_width=1_000_000, pip_threshold=0.9)` – builds windows and derives credible sets from per-SNP PIPs.
- `pip_to_pvalues(pips, prop_null)` – converts PIPs to pseudo p-values following the legacy heuristic.
- `mcmc_samples_sparse_matrix(samples)` – wrapper around `McmcSamples.sparse_data()` returning NumPy arrays.
- `mcmc_samples_to_csr(samples)` – materialises SciPy CSR matrices when SciPy is available.
- `aggregate_parameter_blocks(samples)` – summarises non-SNP parameter blocks from a list of `McmcSamples`.
- Internal helpers for credible-set maths: `_sum_columns`, `_column_vector`, `_normalised_mcmc_matrix`.
- `compute_credible_sets(...)` – Python port of `GCTB::calcCredibleSets`, working with `Data` instance and `McmcSamples`.
- `compute_window_pip_trace(data, snp_effects, *, window_width, step_size, snp_results)` – builds window-level PIP traces from sparse SNP samples.
- `solve_snp_effects_cg(data, *, lambda_, output_path)` – uses `Data.get_zpz_sparse_matrix()` and SciPy CG solver to recover SNP effects; writes a tab-separated report.
- `run_posthoc_stratify(...)` – orchestrates post-hoc annotation analysis by loading LD matrices, annotations, and existing MCMC samples, then invoking `build_posthoc_model` and `MCMC.run`.
- `tune_eigen_cutoff(...)` – automates eigen cutoff tuning by repeatedly calling `Data.read_eigen_matrix_binary_file_and_make_wq`, `Data.init_variances`, and `build_model_summary`.

### `diagnostics.py`

Optional diagnostics (SciPy-dependent):

- `geweke_test(samples, first=0.1, last=0.5)` – computes Geweke z-scores.
- `effective_sample_size(samples, max_lag=None)` – ESS estimate via autocorrelation sums.
- `check_convergence(results, chain_length, burnin, verbose=True)` – iterates over `McmcSamples` objects returned by the sampler, reporting ESS and Geweke stats (skips SNP effect blocks).
- `recommend_chain_length(diagnostics)` – simple rule-based advisory on extending chain length.

### `cli.py`

Implements the `gctb` command-line interface backed by `_core`.

Helpers:

- `load_plink_data(...)` – repeatedly calls `Data.read_*` methods to load individual-level data and optional annotations.
- `_parse_float_list(value)` – utility for comma-separated float arguments.
- `load_summary_data(...)` – summary-statistics analogue of `load_plink_data`, handling LD data and filtering.
- `save_parameter_results(results, output_prefix)` – writes non-SNP parameter summaries by reading fields from `McmcSamples`.
- `save_snp_results(results, data, output_prefix)` – uses `Data.output_snp_results` or pandas-style writing fallback logic.
- `write_credible_set_outputs(output_prefix, results)` – dumps credible set summaries produced by `workflows.compute_credible_sets`.

Entry points:

- `main()` – top-level Click group.
- `ldmatrix` – subgroup for LD utilities.
  - `make_ldmatrix(...)` – delegates to `load_plink_data` and `Data.make_ld_matrix`/related functions to construct LD matrices.
- `bayes(...)` – individual-level Bayes analysis; builds model via `_core.build_model` and runs sampling (optionally with progress bar), then saves outputs and diagnostics.
- `sbayes(...)` – summary statistics analysis. Uses `load_summary_data`, `build_model_summary`, optionally runs multiple chains (`workflows.run_multi_chain_sbayes`), and integrates credible set outputs.
- `posthoc_stratify(...)` – wraps `workflows.run_posthoc_stratify` for post-hoc annotation analysis.
- `tune_eigen(...)` – command-line interface to `workflows.tune_eigen_cutoff`.
- `benchmarks/python_vs_cpp/run_benchmarks.py` – orchestrates end-to-end comparisons between the Python CLI and the legacy C++ executable.

All CLI commands surface their parameters to Click and ultimately call the underlying C++ wrappers provided by `Data`, `Model`, and `MCMC`.

## Interactions and Typical Workflow

1. **Data preparation**: `load_plink_data` or `load_summary_data` create a `_core.Data` object by invoking the `read_*` functions above.
2. **Model build**: `build_model` / `build_model_summary` compute an appropriate C++ model based on the data configuration.
3. **Sampling**: `MCMC.run` (or `run_mcmc`) iterates in C++, returning lists of `McmcSamples` objects.
4. **Post-processing**: Python helpers in `workflows` and CLI functions consume `Data` and `McmcSamples` to write outputs (`Data.output_*`), compute credible sets, or perform diagnostics (`diagnostics.py`).
5. **Utilities**: LD matrix creation, pruning, and tuning call back to the corresponding `Data` methods to reuse the C++ implementations.

## Suggestions for Further Refinement

- Consider auto-generating API documentation from the binding definitions (e.g. with `sphinx-autodoc` plus custom templates) to keep the Python and C++ function lists synchronised.
- Many binding wrappers catch `std::string` exceptions only to convert them to `RuntimeError`; adding richer exception types (e.g. custom Python exceptions) could improve downstream error handling.
- High-level workflows could benefit from type annotations and docstrings that reference this document to help users discover the correct call order.
- If additional bindings are added, grouping related `Data` methods into smaller façade objects in Python might reduce the cognitive load for users who only need specific workflows.


