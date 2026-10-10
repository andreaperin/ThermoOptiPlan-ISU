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
- Scaling the inputs to [-1, 1] (training min/max) matters: κ ~ 1e-13 next to ρ ~ 2700 makes the LP
  badly conditioned. The scaled IPM gives narrower intervals.

## Settings and files

- `scripts/surrogates/settings.jl`: `DATASET`, `TEST_PERCENT`, `SEED`, `TAG = "<DATASET>_test<TEST_PERCENT>"`,
  `train_test_split()`. It is shared by ipm.jl, lbfm.jl and figures_ipm.jl.
- Fit: `scripts/surrogates/ipm.jl` → `data/surrogates/ipm_<scaled|not_scaled>_<TAG>.jld2`
  (model, test, lb, ub; the scaled file also has lo, hi).
- Plots: `scripts/auxiliary/figures_ipm.jl` → `data/figures/ipm_<scaled|not_scaled>_<TAG>_<plot>.pdf`.

## The three plots

1. **Parity** (`_parity.pdf`): x = OGS crossing year, y = IPM interval (midpoint with the interval as error bar).
   Blue = test point inside the interval, red = outside. Points on the diagonal with short bars = accurate
   and tight. The title shows the test coverage.
2. **Band** (`_band.pdf`): the test points sorted by OGS value, plotting `lb - y` and `ub - y` (interval minus
   OGS value). The OGS value is inside when the band crosses 0. Red crosses mark the misses. It shows how
   wide the interval is and by how much the misses fall outside.
3. **Reliability** (`_reliability.pdf`): the theoretical guarantee (see below), with the measured test miss
   rate as a red dashed line.

## Meaning of the reliability plot

Take a new OGS run with random inputs from the same distributions. Let ε be the probability that its
crossing year falls **outside** the interval. The true ε is unknown, because the interval comes from a random
training sample. The curve shows, for each ε:

    confidence 1 - β = P(the IPM misses at most a fraction ε of new runs)

In scenario theory only a few training points actually decide the bounds, at most 2n = 32 here. So

    β = reliability(ipm, ε) = cdf(Binomial(N, ε), 2n - 1)

depends only on N (training points) and n (terms), not on the test data.

- How to read it: pick a confidence (e.g. 0.95) and read off the ε where the curve reaches it. That ε is
  the guaranteed worst-case miss rate.
- The curve rises sooner with more training points or fewer terms.
- The red line (measured test miss rate) should lie at or to the left of where the curve approaches 1.
- The guarantee is conservative, so the real miss rate is usually lower.
- The ε axis adapts to N: it runs up to where 1 - β ≈ 1.

| ε (max miss rate) | v1, N = 775 | v2, N = 98 |
|---|---|---|
| 0.05 | 0.89 | 0.00 |
| 0.10 | 1.00 | 0.00 |
| 0.20 | 1.00 | 0.002 |
| 0.30 | 1.00 | 0.32 |
| 0.40 | 1.00 | 0.95 |
| 0.50 | 1.00 | 1.00 |

## Results (seed 1, linear basis)

| | train N | scaled: coverage / mean width | not scaled: coverage / mean width | β(ε = 0.1) |
|---|---|---|---|---|
| v1_969, test 20% | 775 | 98.5% / 0.22 y | 99.0% / 0.41 y | 3.5e-10 |
| v1_969, test 90% | 97 | 72.7% | 79.8% | ≈ 1 |
| v2_123, test 20% | 98 | 80% / 0.18 y | 84% / 0.24 y | ≈ 1 |

- v1 with 775 points: "≥ 90% of new runs inside" holds with confidence ≈ 1. The measured miss rate
  (1.5%) agrees.
- v2 with 98 points: the best statement is "with 95% confidence, ≥ 60% of new runs inside" (ε = 0.4).
  The measured miss rate of 20% is consistent with it, but the guarantee is weak. The predictions themselves
  are good (points on the diagonal); the interval is just too tight. More runs, or fewer terms, tighten it.
- v2_123 = the 123 complete runs of the 128 on the cluster. 5 OGS runs stopped early (rows 32, 73, 75, 78, 86
  of `samples` in `dataset_3layers_v2_128.jld2`). Their inputs are not extreme, so dropping them does not bias
  the fit.
