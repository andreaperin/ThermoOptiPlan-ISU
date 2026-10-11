# Sensitivity of the crossing year (scripts/sensitivity)

## SRC² from the data (src.jl, figures_src.jl)

A linear regression on **all** the runs of `DATASET` (no test set needed), with standardized inputs:

    y ≈ β₀ + Σ βᵢ zᵢ,   zᵢ = (xᵢ − mean) / std,   SRCᵢ = βᵢ / std(y)

- SRCᵢ (standardized regression coefficient): how many standard deviations of y move when xᵢ moves by one of
  its own standard deviations. Unitless, so comparable across inputs; the sign is the direction of the effect.
- SRCᵢ²: share of the variance of y caused by xᵢ alone.
- Standardizing also keeps the least squares well conditioned (κ ~ 1e-13 next to ρ ~ 2700).

`src.jl` saves `data/sensitivity/src_<DATASET>.jld2` (inputs, src, src2, r2); `figures_src.jl` reads both
datasets and saves `data/figures/sensitivity_src.pdf` (bars per input, v1 and v2 side by side, R² in the legend).

### R² and Σ SRC²

R² (fitted, on all runs) is the share of the variance of y explained by the whole linear model.

    R² = Σ SRCᵢ² + Σᵢ≠ⱼ SRCᵢ · SRCⱼ · corr(xᵢ, xⱼ)

Σ SRC² counts each input's own contribution; R² also includes the cross terms from sample correlations
between inputs. Independent inputs have them only by chance, so they shrink with N.

| | R² | Σ SRC² |
|---|---|---|
| v1 (969 runs) | 0.99879 | 0.99888 |
| v2 (123 runs) | 0.99891 | 0.98366 |

SRC² ≈ first-order Sobol index when both hold: R² ≈ 1 (linear model enough: 1 − R² is the nonlinear part)
and Σ SRC² ≈ R² (inputs uncorrelated in the sample). Both hold on v1; on v2 roughly (±0.01-0.02), ranking unchanged.

### Results

| Input | Layer | SRC² v1 | SRC² v2 | sign |
|---|---|---|---|---|
| cp (heat capacity) | middle | **0.782** | **0.774** | + |
| cp | bottom | 0.100 | 0.087 | + |
| ρ (density) | middle | 0.081 | 0.077 | + |
| φ (porosity) | middle | 0.017 | 0.024 | + |
| ρ | bottom | 0.010 | 0.010 | + |
| λ (all layers), φ bottom | | ≈ 0.002 each | ≈ 0.002-0.003 each | λ top −, rest + |
| cp, φ, ρ top | top | < 0.001 | < 0.001 | + |
| κ (all layers) | | ≈ 0 | ≈ 0 | |

- The heat stored in the rock (cp, ρ) dominates: more heat capacity, slower cooling, later crossing.
- The middle layer matters most: it carries the largest pumping share (43%); the top layer (24.5%) barely matters.
- κ ≈ 0 also in v2, where it varies per layer: the flow split between layers is fixed by the pumping shares in
  the post-processing (option A), so κ cannot change how much water each layer carries. With flow shares
  from OGS (option B, future work) κ could matter.

## The old PCE Sobol indices: do not use them

The 512-run PCE (degree 4: 3876 coefficients for 512 runs) gives flat indices: first order ≈ 0.01 for every
input (sum ≈ 0.18), total effect ≈ 0.2 for every input. That is the signature of an overfitted,
underdetermined PCE, in contradiction with SRC² (one input explains 78%). Show them, at most, as evidence
that the PCE was underdetermined.

## Next: Sobol indices on the surrogates

`sobolindices(models, inputs, output, MonteCarlo(n))` (UQ.jl, Saltelli estimators, first order and total
effect with bootstrap errors) on the best surrogate per dataset (LBFM deg 2 for v1, spline k = 3 or top-k
degree 2 for v2), sampling from `input_distributions()`. The surrogates are fitted on scaled inputs, so the
model must be wrapped: `Model(df -> (s = scale(df); evaluate!(surrogate, s); s[!, :crossing_year]), :crossing_year)`.
Expected to agree with SRC² (first order), with total effects showing the small interactions.
