# Calibration with wholly null annotations

Added 20 overlapping binary null annotations to each existing simulation, retaining the original 40 causal annotations, original coefficients, and GWAS summary statistics. Each new annotation has exactly zero effect at all four probability levels. These shared nulls are a stress test: default RD models inclusion separately at each level rather than imposing shared annotation-level sparsity.

Both generating scenarios were fitted using RD and RC, with four chains per model, 10,000 iterations, 2,000 burn-in iterations and thinning of 5. Each dataset has 3,000 independent SNPs and a summary-data sample size of 5,000,000. The RC generating scenario uses fixed normal variances 0.01, 0.03, 0.06 and 0.1, rather than drawing the variance from the RC hyperprior.

## Results

RD inclusion PIP here means the probability that at least one of the four annotation coefficients is exactly nonzero, calculated from synchronized posterior draws.

| Generating effects | Mean RD PIP on null annotations | False positives at PIP ≥0.95 | Observed false discovery proportion | PIP-predicted FDR |
|---|---:|---:|---:|---:|
| RD plus shared nulls | 0.9983 | 20/20 | 33.33% | 0.12% |
| RC plus shared nulls | 0.9988 | 20/20 | 33.33% | 0.08% |

All 60 annotations exceed PIP 0.99 in both scenarios. At individual levels, mean RD PIPs on wholly null annotations range from 0.7750 to 0.8355 in the RD generating scenario and from 0.7403 to 0.8464 in the RC generating scenario. Thus the issue is present within levels and is amplified by taking the union across four levels.

The fitted RD weights put substantial probability on tiny nonzero variance components. These components can shrink effects strongly while counting as nonzero for inclusion. In these data, effect shrinkage therefore does not imply reliable causal-annotation PIPs. Inclusion calibration must be assessed independently of correlation or RMSE.

RC has a continuous normal prior and no zero-effect spike, so its implied inclusion probability is 1 for every annotation. It does not provide a directly comparable variable-selection PIP. A separate unadjusted sign screen (any level with posterior sign probability ≥0.975 or ≤0.025) selects 17 annotations (16 causal, one null) under RD generation and 23 (22 causal, one null) under RC generation. This is neither an inclusion probability nor an FDR-controlled procedure.

These are descriptive results from one dataset per scenario. RD coefficient maximum R-hat is below 1.017. RC p5 mixing is poor (maximum coefficient R-hat approximately 1.24 and 1.50 respectively), so its sign-screen results are provisional. More replicated simulations and a prior representing shared annotation-level sparsity would be needed to establish calibration more generally.

## Reproduction

```sh
python3 scripts/simulate_rd_annotations.py test/rd_null_annotations --append-null-annotations 20 --include-rc
python3 scripts/simulate_rd_annotations.py test/rc_null_annotations --annotation-model RC --annotation-reference test/rd_multiple_annotations/sim.annot --append-null-annotations 20 --include-rc
python3 scripts/check_annotation_pip_calibration.py test/rd_null_annotations test/rc_null_annotations
Rscript scripts/plot_annotation_pip_calibration.R test/rd_null_annotations test/rc_null_annotations
Rscript scripts/plot_rd_annotation_recovery.R test/rd_null_annotations
Rscript scripts/plot_rd_annotation_recovery.R test/rc_null_annotations
```

Detailed calibration tables and plots are in each output folder as `pip_calibration.md`, `annotation_pip_calibration.tsv`, and `pip_calibration.pdf/png`; effect recovery and chain diagnostics are in `recovery.md`.

## Shared-indicator RD comparison

The existing joint sampler is now exposed by `--sbayes RD_Joint`, which passes `sampleAnnoEffectsIndepOpt=false` in both single-chain and multi-chain paths. It shares one mixture assignment across all four annotation coefficients. Thus a zero assignment sets every coefficient to zero, and a nonzero assignment supplies a common variance for four conditionally independent Gaussian effects. `--joint-annot` is a separate sampling option and is not needed to select this model. RD and RD_indep retain independent assignments.

Four joint chains were run for each of the same two datasets, using the original seeds and iteration settings. All eight fits completed; all saved coefficient draws switch inclusion together across levels; SNP PIPs are finite. Maximum coefficient split R-hat is 1.0143 under RD generation and 1.0104 under RC generation. A short two-chain smoke fit also verified the multi-chain route and shared output labels (`AnnoMix`, `AnnoPIP`).

| Generating effects | Mean null PIP, independent / joint | Brier score, independent / joint | Joint calls at PIP ≥0.95 | Joint false positives |
|---|---:|---:|---:|---:|
| RD plus shared nulls | 0.9983 / 0.8335 | 0.33223 / 0.24388 | 9/40 causal | 0/20 null |
| RC plus shared nulls | 0.9988 / 0.7792 | 0.33252 / 0.21712 | 12/40 causal | 0/20 null |

Joint RD improves screening and probability prediction relative to independent RD's any-level PIP, but remains overconfident for many annotations: at threshold 0.5, it still calls all 60, with observed false discovery proportion 33.3% versus predicted FDR 13.3% and 15.4%. At threshold 0.9, RD generation yields 14 causal and one null call, while RC generation yields 15 causal and no null calls. No false positives at 0.95 in these two datasets does not establish general calibration.

These simulations were not generated from the full joint RD prior: their active annotations can have different component assignments or scales across levels. That distinction matters because the joint model shares the variance assignment as well as inclusion. Joint PIPs use posterior means of conditional inclusion probabilities; independent any-level PIPs use synchronized coefficient draws.

Reproduce fits and plots with:

```sh
python3 scripts/test_rd_joint_calibration.py test/rd_null_annotations test/rc_null_annotations
Rscript scripts/plot_rd_joint_calibration.R test/rd_null_annotations test/rc_null_annotations
```

Each folder contains `rd_joint_calibration.md`, `rd_joint_calibration.tsv`, and `rd_joint_calibration.pdf/png`. The comparison suggests testing a separate sparsity prior for shared inclusion next, while retaining the mixture for effect estimation; replicated calibration experiments are still necessary.
