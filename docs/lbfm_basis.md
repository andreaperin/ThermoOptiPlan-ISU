# Choosing the basis for the LinearBasisFunctionModel (lbfm.jl)

`LinearBasisFunctionModel(df, out, basis, inputs)` (UQ.jl) fits the coefficients by least squares, `β = b(X)' \ y`.

## Available bases

| Basis | Constructor | Terms with 15 inputs |
|---|---|---|
| polynomial | `MonomialBasis(15, p)` | p = 1 → 16, p = 2 → 136, p = 3 → 816 |
| Gaussian radial | `GaussianRadialBasis(C, ϵ)` | one per center (columns of `C`, 15 × m) |
| polyharmonic radial | `PolyharmonicRadialBasis(C, k)` | one per center |

## Rule: number of terms vs training points

Least squares needs clearly more points than terms (at least 2-3 points per term).

| | v1_969 (775 train) | v2_123 (98 train) |
|---|---|---|
| degree 1 (16 terms) | yes | yes (~6 points per term) |
| degree 2 (136 terms) | yes (~5.7 per term) | no: 136 > 98, `\` interpolates and generalises badly |
| degree 3 (816 terms) | no | no |

## Plan

1. `MonomialBasis(length(inputs), 1)`: same basis as the IPM, the fair comparison
   (IPM = interval around a linear model, LBFM = best single linear model).
   The response is close to linear (IPM bands ~0.2 y wide over an 8 y range).
2. Degree 2 on v1 only: if the test RMSE drops clearly vs degree 1, there is curvature/interactions;
   otherwise degree 1 is enough. Always decide on test error, never training error.
3. Skip radial bases for now: they have no constant/linear term (a near-linear response needs
   many centers); you must choose centers (subset of training points, 15 × m) and, for the
   Gaussian, ϵ ≈ 1 / typical distance between points. All training points as centers = interpolation,
   which is PolyharmonicSpline / GP territory (own GP code).

## Practical points for lbfm.jl

- Scale inputs to [-1, 1] as for the scaled IPM: for degree 2 it is essential
  (κ² ~ 1e-26 next to ρ² ~ 7e6). Copy `lo`/`hi`/`scale` from ipm.jl or move them into settings.jl.
- `evaluate!` overwrites the output column with numbers: evaluate on a copy (`pred = scale(test)`).
- Same split: `train, test = train_test_split()` from settings.jl.
- Print: test RMSE and max |error|, RMSE of the IPM midpoint for comparison,
  and the fraction of LBFM predictions inside the IPM interval.
- Save as `lbfm_deg<p>_$(TAG).jld2` so degree 1 and 2 don't overwrite each other.
