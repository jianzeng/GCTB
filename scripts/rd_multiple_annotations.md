# RD simulation with multiple annotation effects

Run the model-matched simulation and both fits with:

```sh
python3 scripts/simulate_rd_annotations.py test/rd_multiple_annotations --include-rc
Rscript scripts/plot_rd_annotation_recovery.R test/rd_multiple_annotations
```

The executed simulation uses seed 20261011, 3,000 independent SNPs and 40
overlapping binary annotations plus an intercept. There are 160 non-intercept
annotation effects across four sequential probit probability levels. At each
level, annotation weights are drawn from Dirichlet(1,1,1,1,1); assignments and
coefficients are drawn from RD's exact-zero/normal mixture with variances
0, 0.0001, 0.001, 0.01 and 0.1. This produces positive, negative and zero
effects, with different weight vectors at different levels.

SNP membership follows those annotation-dependent probabilities. SNP effects
are drawn from normal components with variance multipliers 0, 0.001, 0.01,
0.1 and 1. A fixed base variance targets expected genetic variance 0.5;
drawn SNP effects are not rescaled. Residual variance is 0.5. Gaussian GWAS
summary statistics use total phenotype variance / N noise, matching the sparse
GWAS likelihood, with effective N=5,000,000. Exact identity LD is supplied.
This is a summary-data simulation; individual genotypes and phenotypes are
not generated. The binary LD encoding, annotation dimensions, mixture weights,
zero-effect assignments and summary noise variance were checked.

Both fits use `--hsq-percentage-model false` to match the standard SNP
variance mixture used for generation. The prior comparison used the default
percentage model, so its results are not directly comparable to this study.
Each model has four independently seeded chains of 10,000 iterations,
2,000 burn-in and thinning by 5.

| Level | RD effect RMSE | RC effect RMSE | RD effect correlation with truth |
|---|---:|---:|---:|
| p2 | 0.0655 | 0.0913 | 0.5797 |
| p3 | 0.1324 | 0.2230 | 0.8242 |
| p4 | 0.1442 | 0.1902 | 0.7693 |
| p5 | 0.1366 | 0.6359 | 0.4378 |

RD heritability is 0.486414 versus the realized truth 0.488775. Its nonzero
SNP count is 1204.66 versus the realized truth 1205. Non-intercept annotation
coefficient split R-hat values are at most 1.0097; intercept R-hat values are
at most 1.0213. Neither model has nonfinite SNP PIPs in this simulation.

RD's coefficient RMSE is smaller at every level in this realization. However,
RC's p5 coefficients mix poorly (maximum split R-hat 1.1925), so that level's
RC error estimate is provisional. This is one dataset generated under RD,
not a general performance benchmark.

Mixture-weight recovery is incomplete. For example, the p3 generating weight
on variance 0.1 is 0.528 and its posterior mean is 0.542. At p4 and p5, weights
remain near a uniform mixture despite nonuniform generating weights. Tiny
variance components are difficult to distinguish, and each weight vector
has only 40 annotation assignments to learn from.

Saved outputs in `test/rd_multiple_annotations/` include `truth.json`,
`annotation_truth.tsv`, `snp_truth.tsv`, the annotation and GWAS files,
`annotation_recovery.tsv`, `recovery.md`, and the recovery plot in PDF/PNG.
The recovery report includes generating weights, realized assignment fractions,
posterior weights and chain diagnostics. Intercepts are excluded from coefficient
RMSE and correlation. Its approximate coverage uses mean ± 1.96 SD, not
posterior quantiles.

## Added RC-generated scenario

```sh
python3 scripts/simulate_rd_annotations.py test/rc_multiple_annotations \
  --annotation-model RC --annotation-reference test/rd_multiple_annotations/sim.annot \
  --include-rc
Rscript scripts/plot_rd_annotation_recovery.R test/rc_multiple_annotations
```

This reuses exactly the same 3,000-by-40 annotation matrix. Each level now
has a single Gaussian annotation distribution: `alpha_jc ~ N(0,tau_c^2)`,
with true variances 0.01, 0.03, 0.06 and 0.1 at p2–p5. All 160 coefficients
are Gaussian draws; there is no exact-zero spike or generating mixture
weight vector. These are fixed generating variance parameters, not draws
from RC's inverse-chi-square hyperprior. RC learns its variances in the fit
using its existing hyperprior. The other simulation and fit settings match
the RD-generated scenario. `--rc-variances` can specify other positive
variance values.

| Level | True Gaussian variance | RC effect RMSE | RD effect RMSE |
|---|---:|---:|---:|
| p2 | 0.01 | 0.0838 | 0.0823 |
| p3 | 0.03 | 0.3157 | 0.1363 |
| p4 | 0.06 | 0.2463 | 0.2027 |
| p5 | 0.10 | 0.6441 | 0.2602 |

True heritability is 0.48370; posterior means are 0.481234 for RC and
0.481235 for RD. Neither model has nonfinite SNP PIPs. RD's non-intercept
coefficient split R-hat values are at most 1.0087, and its intercept values
are at most 1.0272. RC has poor p5 mixing: maximum coefficient split R-hat
1.1822 and intercept split R-hat 1.5590. Its p4 intercept also has split
R-hat 1.1050. Consequently RC's higher-level errors remain provisional.

RD's smaller RMSE in this realization does not establish general superiority:
this is one dataset with fixed, relatively small generating variances;
RC's production variance hyperprior also influences its fit, and some
RC chains have not mixed well. This comparison does not establish equivalence
between the learned Gaussian variance model and the fixed-variance mixture.

New truths, data, recovery tables and plots are saved separately in
`test/rc_multiple_annotations/`. The original RD scenario is preserved.
The updated simulator was checked to reproduce the original RD annotation,
GWAS, truth and binary LD files exactly with the original seed.
