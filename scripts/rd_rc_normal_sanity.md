# Normal-prior RD sanity check

The isolated test disables RD's annotation BayesR prior, fixes all annotation
inclusion probabilities at 1, and restores RC's learned annotation variance
update (including its prior and non-intercept coefficient count). It uses the
retained RD independent annotation sampler, rather than redirecting RD to RC.
The intercept keeps its flat prior. Production sources and defaults are unchanged.

The simulation uses seed 20261010, 1,000 individuals, 200 independent population
SNPs with finite-sample LD, and 27 realized causal effects. The phenotype has
target heritability 0.5. One binary annotation increases the probability of a
causal effect from 0.08 to 0.35. GWAS summary statistics and the complete sample
LD matrix are computed from the same simulated genotypes.

## Identifiable two-component SNP mixture

With SNP variance multipliers 0 and 1, each model ran four independent chains
of 30,000 iterations, discarding 6,000 and saving every fifth iteration.

| Quantity | RC | Normal-prior RD | Difference / combined MCSE |
|---|---:|---:|---:|
| Heritability | 0.350950 | 0.351185 | 0.665 |
| SNP variance parameter | 0.0112481 | 0.0112714 | 0.533 |
| Number of nonzero SNPs | 33.3045 | 33.3093 | 0.043 |
| Annotation intercept | -1.38952 | -1.38612 | 0.706 |
| Enriched annotation coefficient | 1.13303 | 1.12441 | 1.267 |

Posterior SNP effect means correlate at 0.99999523; SNP PIPs correlate at
0.99999753, with maximum absolute PIP difference 0.0023185. Ordinary split
R-hat values for the compared parameters are at most 1.0010. MCSEs use 20
contiguous batches per chain. These fits agree within Monte Carlo uncertainty;
they are not identical sample paths because RD consumes additional random draws.
This is an equivalence check, not a claim of accurate recovery of true heritability.

## Default five-component SNP mixture

Four chains per model of 10,000 iterations (2,000 burn-in) are inconclusive.
SNP effect means correlate at 0.99997551, but RC's nonzero SNP count has split
R-hat 1.8617 and RD's has 1.1832. The pooled nonzero counts differ substantially,
and 250 of the 400 pooled per-model SNP PIPs are nonfinite. Heritability means
are 0.353466 and 0.351359. This run does not establish equivalence for the
default five-component configuration. The observed numerical symptoms and
mixing need a separate investigation; the simpler mixture result does not
resolve them.

## Production RD with unknown annotation inclusion mass

The production binary was then run on exactly the same simulated data, with
the same four seeds and iteration settings. It learns separate annotation
mixture weights at each conditional probability level. Annotation pi is the
total nonzero mass, `1 - w_zero`, with a Dirichlet(1,1,1,1,1) prior on weights.

For the converged two-component SNP comparison (30,000 iterations per chain):

| Quantity | RC | Normal RD, pi=1 | Production RD, learned pi |
|---|---:|---:|---:|
| Heritability | 0.350950 | 0.351185 | 0.349175 |
| Nonzero SNP count | 33.3045 | 33.3093 | 32.4189 |
| SNP variance parameter | 0.0112481 | 0.0112714 | 0.0115345 |
| Annotation intercept | -1.38952 | -1.38612 | -1.09816 |
| Enriched annotation coefficient | 1.13303 | 1.12441 | 0.317349 |

Production RD's annotation pi averages 0.818108; the enriched annotation's
conditional nonzero PIP averages 0.911491. These are different quantities.
With just one non-intercept annotation, the mixture weights have very little
data to learn from. The pooled weights are approximately
(0.1819, 0.1813, 0.1824, 0.1873, 0.2671), in ascending variance order.

Compared with normal RD, production RD's SNP effect means correlate at
0.99871873 and SNP PIPs at 0.95876601. Mean absolute PIP difference is 0.05353;
maximum difference is 0.19554. The reported production RD parameters have
ordinary split R-hat values at most 1.0006. Thus annotation effects change
substantially, while the overall heritability and SNP effect estimates remain
close. This does not isolate the effect of learning pi: production RD also
replaces the learned single-normal annotation prior with fixed variances up
to 0.1, which imposes stronger shrinkage.

With the default five-component SNP mixture (10,000 iterations per chain),
production RD remains inconclusive: its nonzero SNP count has split R-hat
2.7719 and the first annotation intercept has split R-hat 4.9456. All 200 pooled
production RD SNP PIPs are nonfinite. Its heritability mean is 0.350209, but
component allocations and annotation coefficients cannot be reliably compared
from these chains.

Reproduce the learned-pi comparisons after creating the baseline fits above:

```sh
python3 scripts/test_rd_unknown_pi.py test/rd_rc_normal_two test/rd_unknown_pi_two --two-components
python3 scripts/test_rd_unknown_pi.py test/rd_rc_normal_five test/rd_unknown_pi_five --iterations 10000
```

Detailed outputs from the executed comparisons are in
`test/rd_unknown_pi/two_component/comparison.md` and
`test/rd_unknown_pi/five_component/comparison.md`.

### Longer multi-component RD run

The full five-SNP-component production RD model was also run for 100,000
iterations per chain, with four independent seeds, 20,000 burn-in and thin 5.
The enriched annotation coefficients at levels p2 through p5 averaged
0.285947, 0.028482, 0.001738 and 0.000877. Their annotation PIPs averaged
0.899717, 0.809118, 0.800731 and 0.801472, respectively. The higher levels stay
close to the annotation prior's expected nonzero mass of 0.8; these PIPs
are not evidence of strong effects at those levels.

Heritability averaged 0.349284, but the nonzero SNP count still had split
R-hat 1.1860 and the first annotation intercept had split R-hat 1.2428.
All 200 pooled SNP PIPs were nonfinite. Increasing chain length therefore
does not resolve the multi-component comparison on this small dataset.
Detailed per-level diagnostics and the reproduction command are saved in
`test/rd_unknown_pi/five_component_long/multi_component.md`.

RC was subsequently rerun with the same 100,000-iteration settings. The
enriched annotation coefficient means compare as follows:

| Probability level | RC | Production RD | RC coefficient split R-hat | RD coefficient split R-hat |
|---|---:|---:|---:|---:|
| p2 | 0.785661 | 0.285947 | 1.1487 | 1.0109 |
| p3 | 0.359244 | 0.028482 | 1.0928 | 1.0338 |
| p4 | 0.052256 | 0.001738 | 1.0032 | 1.0000 |
| p5 | 0.037758 | 0.000877 | 1.0021 | 1.0001 |

RD shrinks these coefficient means substantially. However, this is a
descriptive comparison of unconverged chains: the maximum intercept split
R-hat across probability levels is 3.2550 for RC and 2.2275 for RD. RC's
annotation signal can move between probability levels across chains. The
converged two-component comparison remains the stronger evidence of the
shrinkage difference (1.13303 for RC versus 0.317349 for RD).

The matched long-chain annotation comparison is saved in
`test/rd_unknown_pi/five_component_long/rc_rd_annotation_comparison.md`.
The separate `comparison.md` now uses long RC and production RD fits but
the earlier, shorter normal-RD baseline; it is not a converged comparison.

## Reproduce

After building `scr/gctb`, with `EIGEN3_INCLUDE_DIR` and `BOOST_LIB` configured:

```sh
python3 scripts/test_rd_rc_normal.py test/rd_rc_normal_two --two-components --iterations 30000
python3 scripts/test_rd_rc_normal.py test/rd_rc_normal_five --iterations 10000
```

The script writes an isolated modified `model_normal.cpp`, compiles and links
it with the existing project objects, generates PLINK genotypes and summary
statistics, and runs full GCTB fits. It saves sampled parameters, logs, posterior
summaries, and `comparison.md`. The other linked objects must come from the
current build. The RD constructor's production AnnoMix and AnnoPi output labels
are stale in this test build and are excluded from the comparison; the actual
inclusion vector passed to the annotation sampler is fixed at 1. The script
works around a legacy command-line parsing omission for RD's `--pi` by selecting
RC while parsing the SNP mixture, then selecting RD before running the fit.
