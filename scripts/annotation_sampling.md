# Annotation sampling

The checkpoint before these changes is `2917ea3`.

## Usage and defaults

Single-site annotation Gibbs is the default for RC, RD and RD_indep, including
GWFM. Add `--joint-annot` to explicitly enable full blocking for RC:
RD uses BayesR Gibbs regardless of this flag (see below).

```
gctb --gwfm RC ... --joint-annot
gctb --sbayes RD ... --joint-annot
gctb --sbayes RD_indep ... --joint-annot
```

An option file can instead specify `jointAnno Yes` (`No` is the default).
The flag requires annotation data and annotation-effect estimation enabled.
The old experimental `GCTB_RC_ANNO_BLOCKED` and `GCTB_RD_ANNO_TGS` environment
switches no longer control sampling: the command-line/option-file choice is
authoritative. Original annotation ordering and output labels are preserved.

## RC

With `--joint-annot`, RC draws ALL Gaussian annotation coefficients jointly,
conditional on the probit latent responses. The intercept is included with zero prior
precision. The other coefficients have their existing component-specific
Gaussian prior. Cholesky solves are performed in double precision; the
existing float coefficient/residual storage is retained.

Omit `--joint-annot` to use the original scalar Gibbs sweep for comparison.
This does not change the prior or likelihood.

## RD and RD_indep

RD and its RD_indep alias now use an annotation BayesR prior:

```
alpha_jc ~ w_0 delta_0 + w_1 N(0,0.0001) + w_2 N(0,0.001) + w_3 N(0,0.01) + w_4 N(0,0.1)
w_c ~ Dirichlet(1,1,1,1,1)
```

These are fixed absolute variances on the probit coefficient scale, not
multipliers of a learned base variance. Intercepts retain their flat prior.
Weights start at (0.9, 0.025, 0.025, 0.025, 0.025); every sweep samples all five
assignments using integrated Gaussian likelihoods and then updates weights
from assignment counts, excluding intercepts. Empty probability levels draw
non-intercept coefficients from the prior rather than forcing them to zero.

The default independent model learns separate weights for each probability
level. `AnnoMix_p*` records the five weights, `AnnoPi_p*` is their total
nonzero mass, and `AnnoPIP_p*` records conditional nonzero probabilities.
The `RD_Joint` model (`--sbayes RD_Joint`) passes
`sampleAnnoEffectsIndepOpt=false` and shares an assignment and weight vector across
levels; effects are conditionally independent Gaussian draws with the
selected variance. Its weight output is `AnnoMix`.

With `--write-mcmc-txt`, sampled annotation inclusion indicators are saved
in `<out>.mcmcsamples/AnnoDelta.mcmcsamples.txt` for RD_Joint, or
`AnnoDelta_p2.mcmcsamples.txt`, etc., for independent RD/RD_indep. Each
entry is 0 for the exact-zero assignment or 1 for a nonzero assignment;
these are sampled binary indicators, not conditional inclusion probabilities.
The intercept is the first annotation column and is always 1.

Delta traces follow the same burn-in, thinning and row order as annotation
coefficient traces. Each delta trace is a tab-separated table with a header
`Iteration`, `Chain`, followed by the annotation names. Every data row
contains the iteration and chain followed by its binary inclusion indicators.
Indices are one-based; multi-chain rows are ordered by iteration, then chain.
There are no separate row or column metadata files.

All other parameter text traces written with `--write-mcmc-txt` use this
same format, including annotation effects, mixture weights and
`CoreParameters.mcmcsamples.txt`. Every retained chain draw occupies one
row. This changes the legacy headerless parameter-set text format; the
simulation analysis scripts accept both formats through `trace_values`.

Binary SNP-effect samples retain their sparse binary encoding. Whenever
they are written (`--write-mcmc-bin`, or automatically under `--gwfm`),
`SnpEffects.mcmcsamples.bin.columns.txt` maps one-based `Column` to SNP
`Name`, and `SnpEffects.mcmcsamples.bin.rows.txt` maps one-based `Row` to
`Iteration` and `Chain`. The mapping includes all retained rows, including
rows without nonzero effects. Sparse indices inside the binary remain
zero-based. The metadata is written even when text output is disabled.
This preserves joint inclusion draws needed for annotation credible sets.
The indicators are recorded directly from mixture assignments without extra
random draws. Run `python3 scripts/test_annotation_delta_output.py
test/annotation_delta_output` to check single-chain and multi-chain alignment
for both independent and joint RD using the existing null-annotation dataset.

RD no longer updates SigmaSqAnno. RC retains its existing variance prior.
The existing annotation TGS and swap kernels assume a single slab and are
not used with BayesR. `--joint-annot` therefore leaves RD on ordinary Gibbs
and prints an explanatory message. The SNP sampler is unaffected.
The old single-slab kernels are retained for regression tests only.

### Legacy single-slab annotation TGS

The following describes the retained test kernels, not the BayesR sampler.

Within the full block of non-intercept annotations:

1. Calculate each exact conditional inclusion log-odds `l_j` and conditional
   Gaussian slab. Default RD applies this separately to each probability
   level. The retained internal shared-indicator kernel instead sums the
   component-specific log Bayes factors and updates the shared indicator and
   its whole effect vector together.
2. Define `g_j` by replacing `l_j` with `l_j/2`, retaining the exact conditional
   Gaussian slab. Select a coordinate with weight `w_j = g_j(x_j|x_-j) /
   p_j(x_j|x_-j)`. Gaussian densities cancel in this weight.
3. Draw the selected coordinate from `g_j`, then recompute the weights at
   the proposed state. Accept with `min(1, Z(x)/Z(x'))`, where `Z=sum(w_j)`.
   The posterior and proposal density ratios cancel to this expression.
4. Perform one proposal per coordinate in the block. Only accepted states
   update the residual. The intercept is left for its ordinary Gibbs update.

This is a posterior-invariant MH kernel, not uncorrected importance-tempered
sampling. No importance weights are needed for posterior summaries. A final
ordinary Gibbs sweep alone would NOT justify uncorrected tempered draws.

## Cost and limitations

Full cross-products for the first mixture component are cached per chain,
because this component always uses all SNP rows. Later components' selected
SNP rows change, so their cross-products are recomputed. Cached storage is
O(number of annotations squared) per chain. RC factorization costs O(p^3);
RD stores component-specific p-by-p cross-products and performs p proposals,
each scanning p conditional weights. There is no additional double-precision
SNP-by-annotation design matrix. TGS proposals operate on cross-products/residual projections;
they do not rescan all SNPs for each proposed state.

There are no between-block correlations to manage, so correlation clustering
is unnecessary. Full blocking may be expensive for many annotations and many
chains. No runtime speedup is assumed: compare effective samples
per second and between-chain diagnostics on representative data.

Both kernels are independent of float/q8/q16 LD storage. The existing handling
of mixture components with no eligible SNPs is unchanged.

## Tests

`scripts/test_annotation_sampler.cpp` checks conditional Gaussian means and
covariances, scalar and multi-coefficient blocks, and the corrected TGS kernel
alone against exactly enumerated two-annotation posteriors for independent
and shared multi-component selection. This includes PIPs and first/second
moments, multiple prior inclusion probabilities, and residual consistency.
Compile it with the project's Eigen/Boost include paths and `scr/stat.cpp`,
without `-DNDEBUG` (the tests use assertions).

`scripts/test_annotation_integration.cpp` exercises default and full-block
updates with 1, 4 and 35 annotations, checking a single full-sized cache,
empty mixture components, duplicated columns, residual-related summaries and
nonzero counts. `scripts/test_joint_anno_options.cpp` checks parsing/defaults.
`scripts/test_rd_default.sh` compares RD with its RD_indep alias for one and
two chains, with single-site Gibbs and full annotation TGS, and checks the
per-level output labels.

## BayesR validation

After building with `make -C scr -j4`, run `bash scripts/test_rd_bayesr.sh`.
This checks integrated assignment probabilities, prior mixture moments,
independent/shared updates, empty levels, collinear annotations, intercept-only
inputs, and the retained RC/legacy RD kernels.
