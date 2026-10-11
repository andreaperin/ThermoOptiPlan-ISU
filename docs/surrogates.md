# Surrogates of the 3-layer OGS model: overview

Output: `crossing_year`, the year the weighted extraction temperature has dropped by 1 °C (about 54-63 y).
Inputs: 15 random material parameters, λ / cp / ρ / φ / κ for the top, middle and bottom layer
(`models/README.md`, distributions in `src/inputs.jl`).

## Datasets

| Dataset | Runs | Distributions | Train / test (20%) |
|---|---|---|---|
| `v1_969` | 969 complete of 1024 (NAS run, June) | old, `inputs_v1()` | 775 / 194 |
| `v2_123` | 123 complete of 128 (cluster, October) | current per-layer, `inputs_v2()` | 98 / 25 |

The 5 failed v2 runs (rows 32, 73, 75, 78, 86 of `samples` in `dataset_3layers_v2_128.jld2`) stopped early
in OGS; their inputs are not extreme, so dropping them does not bias the fit.
v1 and v2 have different input distributions: every model is fitted per dataset, never mixed.

## Settings shared by all surrogates

`scripts/surrogates/settings.jl`:

- `DATASET` (`"v1_969"` or `"v2_123"`), `TEST_PERCENT = 20` (always), `SEED = 1`, `TAG = "<DATASET>_test20"`
- `train_test_split()`: same split for every surrogate
- `scaler(train, inputs)`: inputs to [-1, 1] with the training min/max
- `input_distributions()`: `inputs_v1()` or `inputs_v2()`, from `DATASET`

Fits in `scripts/surrogates/`, results in `data/surrogates/<model>_<variant>_<TAG>.jld2`,
figures from `scripts/auxiliary/figures_*.jl` in `data/figures/`.

## Results (test set, 20%)

| Surrogate | v1_969 RMSE [y] | v2_123 RMSE [y] | Notes |
|---|---|---|---|
| LBFM monomial deg 1 (16 terms) | 0.050 | 0.062 | R² 0.9989 / 0.9982 |
| LBFM monomial deg 2 (136 terms) | **0.007** | not possible (136 > 98) | R² 0.99998 |
| Polyharmonic spline k = 1 / 2 / 3 | 0.024 / 0.014 / 0.009 | 0.057 / 0.052 / **0.047** | interpolates the training points |
| IPM deg 1, scaled | coverage 98.5%, width 0.22 | coverage 80%, width 0.18 | ε at 95% confidence: 5.4% / 40% |
| PCE deg 4 (earlier runs 128/256/512, v1) | not on this test set yet | - | 3876 coefficients for at most 512 runs |

Exploratory (scratch analysis, not in the scripts yet): degree 2 only on the top-k inputs, linear on the rest:
top-6 (37 terms) gives 0.018 y on v1 and **0.026 y on v2**, better than the spline on v2.

## Key message

The response is almost quadratic. On v1, a parametric degree-2 model is already very good (0.007 y).
On v2, the parametric model is limited by the number of terms (degree 2 cannot be fitted with 98 points),
so non-parametric models are the natural next step: the spline now, the GP next time.
Knowing which inputs matter (sensitivity) also lets a small parametric model work on v2 (top-k degree 2).

## Parametric vs non-parametric

| | parametric (LBFM, IPM, PCE) | non-parametric (spline, GP) |
|---|---|---|
| number of parameters | fixed by the chosen degree / basis | grows with N (spline: N + 16) |
| more data gives | more precise coefficients | a more flexible model |
| main choice | form / degree | kernel / shape |
| typical failure | too simple a form, or too many terms for the data | follows noise, cost ~N³ |
| model after fitting | coefficients only | coefficients + all training points |

## Reporting

- test RMSE [y]: the physical error (0.05 y ≈ 18 days, 0.007 y ≈ 2.5 days)
- test R² (Q²): unitless share of variance explained on unseen runs; always high here because the crossing year
  varies by ~9 y, so it needs the RMSE next to it
- fitted R² (on training data) only says whether a model form is adequate (used for SRC², see `sensitivity.md`)
- parity / error plots: where the errors are

## Next steps

Sobol indices on the surrogates, top-k degree-2 basis for LBFM and IPM, spline on the top-k inputs,
the old PCE on the v1 test set, a comparison table with fit times, and the GP.
