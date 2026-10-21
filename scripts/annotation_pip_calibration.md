# Annotation inclusion calibration checks

RC has a continuous Gaussian annotation prior, so its posterior probability of an exactly nonzero coefficient is 1. It does not output a nontrivial annotation inclusion PIP. Its PostProbAboveZero is a sign probability, not an inclusion probability. RD explicitly reports per-level inclusion PIPs.

Brier score is mean squared error of predicted binary probabilities; it combines calibration and discrimination, so it is not a standalone calibration measure. The RC values below use its mathematically implied PIP=1, not a sign-probability output.

These are two individual datasets with only 40 annotations each, and no wholly null annotations. They cannot establish a general winner for annotation-level inclusion calibration. Mean PIP agreement can hide differences between probability levels.

## Generated under RD

Causal coefficient means an exactly nonzero generating effect at that probability level. All 40 annotations are causal at some level; there are no wholly null annotations.

Across 160 level-specific effects: observed causal fraction 0.83125; mean RD PIP 0.83493; RD Brier score 0.13357; RC implicit nonzero-probability-1 Brier score 0.16875.

| Level | Observed causal fraction | Mean RD annotation PIP |
|---|---:|---:|
| p2 | 0.975 | 0.822 |
| p3 | 0.975 | 0.904 |
| p4 | 0.775 | 0.802 |
| p5 | 0.600 | 0.810 |

RD posterior probability of a nonzero effect at any level, computed from synchronized coefficient draws: mean 0.99932, range 0.99813–1.00000.

## Generated under RC

Causal coefficient means an exactly nonzero generating effect at that probability level. All 40 annotations are causal at some level; there are no wholly null annotations.

Across 160 level-specific effects: observed causal fraction 1.00000; mean RD PIP 0.86189; RD Brier score 0.02255; RC implicit nonzero-probability-1 Brier score 0.00000.

| Level | Observed causal fraction | Mean RD annotation PIP |
|---|---:|---:|
| p2 | 1.000 | 0.806 |
| p3 | 1.000 | 0.904 |
| p4 | 1.000 | 0.887 |
| p5 | 1.000 | 0.851 |

RD posterior probability of a nonzero effect at any level, computed from synchronized coefficient draws: mean 0.99966, range 0.99875–1.00000.
