# Polyharmonic spline (spline.jl, figures_spline.jl)

## What it is (UQ.jl)

`PolyharmonicSpline(df, k, output)` = linear model + one radial function per training point:

    s(x) = Σᵢ wᵢ φ(‖x − xᵢ‖) + v₀ + v₁x₁ + … + v₁₅x₁₅,   φ(r) = r^k (k odd) or r^k log r (k even)

It **interpolates**: it passes exactly through every training point, so the training error is 0 and only the
test error says something. It uses every column except `output` as an input, so it gets `scale(train)`
(15 scaled inputs + `crossing_year`); scaling matters because the distances r mix all inputs.

## Why there is no term limit

| | unknowns | equations |
|---|---|---|
| radial weights wᵢ | N | N: s(xᵢ) = yᵢ |
| linear coefficients v | 16 | 16 side conditions (Σ wᵢ = 0, Σ wᵢ xᵢ = 0) |

The system is always square (N + 16), whatever N: it is a non-parametric model, its flexibility grows with
the data (like a GP). K does not add terms, it only changes the shape φ. Cost ~N³ (791 × 791 on v1, instant).

## Why K = 1, 2, 3 only

The system is guaranteed solvable only if the polynomial part is rich enough for φ:

| K | φ(r) | minimum polynomial part | UQ.jl (linear) |
|---|---|---|---|
| 1 | r | constant | ok |
| 2 | r² log r | linear | ok |
| 3 | r³ | linear | ok |
| 4 | r⁴ log r | quadratic | no guarantee |
| 5 | r⁵ | quadratic | no guarantee |

## Files

- Fit: `scripts/surrogates/spline.jl`, `K` and `VARIANT = "k<K>"` at the top, prints the max training error
  (interpolation check) → `data/surrogates/spline_k<K>_<TAG>.jld2` (spline, lo, hi, test, ŷ, rmse, maxerr, r2).
- Plots: `scripts/auxiliary/figures_spline.jl` (copy of figures_lbfm.jl): parity and error for every spline of `TAG`.

## Results (test 20%, RMSE [y])

| | v1_969 (775 train) | v2_123 (98 train) |
|---|---|---|
| LBFM deg 1 (reference) | 0.050 | 0.062 |
| LBFM deg 2 (reference) | 0.007 | not possible |
| spline k = 1 | 0.024 | 0.057 |
| spline k = 2 | 0.014 | 0.052 |
| spline k = 3 | **0.009** | **0.047** |

- k = 3 is the best spline on both datasets.
- v1: close to LBFM deg 2 without being told the model form.
- v2: 25% better than the linear model, where degree 2 cannot be fitted.
- Limits: it cannot smooth (follows any noise in `crossing_year`), all 15 inputs weigh equally in the
  distances (unimportant inputs blur them: a spline on the top-k inputs is a next step), and it gives no
  uncertainty. A GP has a noise term, one length scale per input and a predictive variance: why GP next.
