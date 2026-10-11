# Interval Predictor Model (ipm.jl, figures_ipm.jl)

## What the IPM is (UQ.jl)

`IntervalPredictorModel(df, out, basis, inputs)` fits two models on the same basis, a lower and an upper one,
`lb(x) = b(x)'β_lb` and `ub(x) = b(x)'β_ub`. It solves a linear program (JuMP + Clarabel):
minimise the mean width `ub - lb` over the training points, subject to **every** training point lying
inside `[lb, ub]`.

- Basis: `MonomialBasis(length(inputs), 1)`, i.e. linear, 16 terms with 15 inputs (2n = 32 coefficients in total).
- `evaluate!(ipm, df)` overwrites the output column with `Interval` values: evaluate on a copy and read
  `lb`/`ub` from it (`getproperty.(pred[!, out], :lb)`).
- `propagate_intervals!` is only needed for interval-valued inputs. With precise inputs it just calls `evaluate!`.
- Scaling the inputs to [-1, 1] (`scaler` in settings.jl) matters for the IPM, unlike least squares:
  κ ~ 1e-13 next to ρ ~ 2700 makes the LP badly conditioned, and the scaled IPM gives narrower intervals.

## Files

- Fit: `scripts/surrogates/ipm.jl` → `data/surrogates/ipm_<scaled|not_scaled>_<TAG>.jld2`
  (model, test, lb, ub; the scaled file also has lo, hi). Dataset and split from `settings.jl`.
- Plots: `scripts/auxiliary/figures_ipm.jl` → `data/figures/ipm_<scaled|not_scaled>_<TAG>_<plot>.pdf`.

## The three plots

1. **Parity** (`_parity.pdf`): x = OGS crossing year, y = IPM interval (midpoint with the interval as error bar).
   Blue = test point inside the interval, red = outside. The title shows the test coverage.
2. **Band** (`_band.pdf`): the test points sorted by OGS value, plotting `lb - y` and `ub - y` (interval minus
   OGS value). The OGS value is inside when the band crosses 0; red crosses mark the misses.
3. **Reliability** (`_reliability.pdf`): the theoretical guarantee (below), with the measured test miss rate
   as a red dashed line. The ε axis adapts to N: it runs up to where 1 - β ≈ 1.

## Meaning of the reliability plot

Take a new OGS run with random inputs from the same distributions. Let ε be the probability that its
crossing year falls **outside** the interval. The true ε is unknown, because the interval comes from a random
training sample. The curve shows, for each ε:

    confidence 1 - β = P(the IPM misses at most a fraction ε of new runs)

In scenario theory only a few training points actually decide the bounds, at most 2n (32 here). So

    β = reliability(ipm, ε) = cdf(Binomial(N, ε), 2n - 1)

depends only on N (training points) and n (terms), not on the test data.

- Pick a confidence (e.g. 0.95) and read off ε: the guaranteed worst-case miss rate.
- The red line (measured miss rate) should lie at or left of where the curve approaches 1.
- The guarantee is conservative: the real miss rate is usually lower.

| ε (max miss rate) | v1, N = 775 | v2, N = 98 |
|---|---|---|
| 0.05 | 0.89 | 0.00 |
| 0.10 | 1.00 | 0.00 |
| 0.20 | 1.00 | 0.002 |
| 0.30 | 1.00 | 0.32 |
| 0.40 | 1.00 | 0.95 |
| 0.50 | 1.00 | 1.00 |

## Results (seed 1, linear basis, test 20%)

| | train N | scaled: coverage / mean width | not scaled: coverage / mean width | ε at 95% confidence |
|---|---|---|---|---|
| v1_969 | 775 | 98.5% / 0.22 y | 99.0% / 0.41 y | 5.4% |
| v2_123 | 98 | 80% / 0.18 y | 84% / 0.24 y | 40% |

- v1: "≥ 90% of new runs inside" holds with confidence ≈ 1; the measured miss rate (1.5%) agrees.
- v2: the best statement is "with 95% confidence, ≥ 60% of new runs inside". The measured miss rate (20%) is
  consistent with it but the guarantee is weak. The predictions themselves are good (points on the diagonal);
  the interval is just too tight for so few points.

## Precision vs certainty (more terms)

About 0.2 y of the v1 width comes from the linear form (LBFM deg 2 reaches 0.007 y), so a richer basis gives
narrower intervals. But each term adds two coefficients that can decide the bounds, so the guarantee weakens:

| IPM basis | ε at 95% confidence, v1 (N = 775) | v2 (N = 98) |
|---|---|---|
| deg 1 (16 terms) | 5.4% | 40% |
| deg 2 (136 terms) | 38% | not possible |

A top-k degree-2 basis (see `lbfm.md`) is the candidate compromise on v1; on v2 keep degree 1.
