# Linear Basis Function Model (lbfm.jl, figures_lbfm.jl)

`LinearBasisFunctionModel(df, out, basis, inputs)` (UQ.jl) fits the coefficients by least squares, `β = b(X)' \ y`.
`ResponseSurface` in UQ.jl is the same model with a monomial basis.

## Files

- Fit: `scripts/surrogates/lbfm.jl`, with `DEGREE` and `BASIS = "monomial_deg<DEGREE>"` at the top
  → `data/surrogates/lbfm_<BASIS>_<TAG>.jld2` (lbfm, lo, hi, test, ŷ, rmse, maxerr, r2).
- Plots: `scripts/auxiliary/figures_lbfm.jl` plots every LBFM of the current `TAG`:
  parity (OGS vs LBFM) and error (LBFM - OGS, test points sorted by OGS value, ±RMSE lines).
- Inputs are scaled to [-1, 1] (`scaler`). For least squares this does not change the predictions,
  only the conditioning, which matters a lot for degree 2 (κ² ~ 1e-26 next to ρ² ~ 7e6).
- A guard stops the fit when the basis has more than N/2 terms (at least 2 training points per term).

## Choosing the basis

| Basis | Constructor | Terms with 15 inputs |
|---|---|---|
| polynomial | `MonomialBasis(15, p)` | p = 1 → 16, p = 2 → 136, p = 3 → 816 |
| Gaussian radial | `GaussianRadialBasis(C, ϵ)` | one per center (columns of `C`, 15 × m) |
| polyharmonic radial | `PolyharmonicRadialBasis(C, k)` | one per center |

Least squares needs clearly more points than terms (at least 2-3 per term):

| | v1_969 (775 train) | v2_123 (98 train) |
|---|---|---|
| degree 1 (16 terms) | yes | yes (~6 points per term) |
| degree 2 (136 terms) | yes (~5.7 per term) | no: 136 > 98 |
| degree 3 (816 terms) | no | no |

Radial bases are not used: the response is smooth and global (polynomials fit it well), 15 dimensions with
98 points is far too sparse for radial functions, and the UQ.jl radial bases have no constant/linear term
inside the LBFM (predictions fall back towards 0 away from the centers). `PolyharmonicSpline` (with a linear
term) is the sensible radial option, as a separate model (`spline.md`).

## Results (test 20%)

| | terms | test RMSE | max error | R² |
|---|---|---|---|---|
| monomial deg 1, v2_123 | 16 | 0.062 y | 0.196 y | 0.9982 |
| monomial deg 1, v1_969 | 16 | 0.050 y | 0.175 y | 0.9989 |
| monomial deg 2, v1_969 | 136 | **0.007 y** | 0.026 y | 0.99998 |

Degree 2 is 7× more accurate than degree 1 on v1: the response has real curvature/interactions, and the small
remaining error shows the crossing year is barely noisy.

## Smaller quadratic bases (exploratory, not in lbfm.jl yet)

The degree-2 v1 model is 99.8% linear; the curvature lives in a few interactions (cp2·ρ2, cp2·φ3, cp3·ρ3),
pure squares carry only 0.02%. Degree 2 on the k most important inputs (ranking from `sensitivity.md`),
linear on the rest:

| Basis | terms | v1 RMSE | v2 RMSE |
|---|---|---|---|
| linear | 16 | 0.050 | 0.062 |
| linear + squares | 31 | 0.049 | 0.074 (useless, dropped) |
| top-4 deg 2 | 26 | 0.024 | 0.035 |
| top-6 deg 2 | 37 | 0.018 | **0.026** |
| full deg 2 | 136 | 0.007 | not possible |

On v2 top-6 degree 2 beats the spline (0.047 y): a small, well-chosen parametric model works on small data.
The top inputs were ranked on v1 and reused on v2 (prior knowledge; v2 has only 25 test points).
Implementing it needs a custom basis (a struct `<: AbstractBasis`), since `MonomialBasis` cannot do
"degree 2 on a subset"; the same basis would also work in the IPM.
A physics-based input (layer heat capacity C = (1-φ)ρcp + φρw·cw, 9 inputs) did worse and was dropped.
