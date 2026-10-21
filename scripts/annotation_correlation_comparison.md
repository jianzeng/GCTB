# True versus estimated annotation effect correlations

Pearson correlations across 40 non-intercept annotation coefficients per probability level. Estimates are posterior means pooled over four independently seeded chains. Both scenarios use the same annotation matrix and simulation sizes.

| Generating annotation model | Fitted model | p2 | p3 | p4 | p5 |
|---|---|---:|---:|---:|---:|
| RD | RC | 0.6346 | 0.8063 | 0.5916 | 0.4473 |
| RD | RD | 0.5797 | 0.8242 | 0.7693 | 0.4378 |
| RC | RC | 0.7302 | 0.6444 | 0.6566 | 0.4982 |
| RC | RD | 0.6708 | 0.7746 | 0.7151 | 0.5324 |

Higher-level RC fits have mixing problems, especially p5 (maximum coefficient split R-hat 1.1925 under RD generation and 1.1822 under RC generation). Those RC correlations are descriptive and provisional. RD coefficient split R-hat values are below 1.01 in both scenarios; some intercept values are modestly higher.

Correlation measures ordering and linear association, not coefficient calibration. It should be read alongside RMSE and chain diagnostics in rd_multiple_annotations.md. Strong shrinkage can preserve correlation while changing coefficient magnitude. These are two individual simulated datasets, not a replicated simulation study.

## Why RC has higher p2 Pearson correlation

Under RD generation, true p2 effects have SD 0.06859. RC posterior means have
SD 0.11628, correlation 0.63459 and RMSE 0.09134. RD posterior means have
SD 0.07546, correlation 0.57967 and RMSE 0.06550. RC preserves more of the
linear association but over-disperses the effect estimates. RD is more accurate
under squared-error loss despite lower correlation.

The fitted-versus-true regression slopes are 1.07576 for RC and 0.63768 for RD
under RD generation; under RC generation they are 0.81507 and 0.33556.
RD's mixture can shrink different coefficients by different amounts,
including towards exact zero, altering their linear relationship to truth.
Uniform positive rescaling alone cannot explain the difference, because it
does not change Pearson correlation.

The p2 Pearson differences are stable across the four individual chains:
under RD generation RC ranges 0.62969–0.63991 and RD 0.57697–0.58299;
under RC generation RC ranges 0.72468–0.73427 and RD 0.66225–0.67640.
These ranges describe Monte Carlo variation on the same datasets, not
uncertainty across repeated simulated datasets.

Pearson correlation should also not be interpreted as a ranking metric.
Under RC generation, p2 Spearman rank correlation is higher for RD (0.78462)
than RC (0.73471), despite RD having lower Pearson correlation. Under RD
generation the Spearman values are 0.60038 for RD and 0.62008 for RC.
