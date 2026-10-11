# Next steps (code), after the 6th internal presentation (2026-10-11)

Open points of the surrogate / sensitivity work. All are listed as "Next Steps" in the slides
(`~/Documents/academic/slides/ThermoOptiPlan/6_presentation_INTERNAL/slides.typ`, 24 slides, done).
Background: `surrogates.md` (overview and results), `lbfm.md`, `spline.md`, `ipm.md`, `sensitivity.md`.

Conventions to keep (see `scripts/surrogates/settings.jl`):

- `DATASET` (`"v1_969"` or `"v2_123"`), `TEST_PERCENT = 20`, `SEED = 1`, `TAG = "<DATASET>_test20"`
- fit scripts in `scripts/surrogates/`, figures in `scripts/auxiliary/`, sensitivity in `scripts/sensitivity/`
- inputs scaled to [-1, 1] with `scale, lo, hi = scaler(train, inputs)`; save `lo, hi` with the model
- files `data/surrogates/<model>_<variant>_<TAG>.jld2` with `test, ŷ, rmse, maxerr, r2`
- input distributions: `input_distributions()` (v1 or v2 from `DATASET`, defined in `src/inputs.jl`)

Suggested order: 1 → 2 → 3 → 4 → 5 → 6. The GP (7) is the big one, for the next presentation.

---

## 1. Top-k degree-2 basis for the LBFM (and IPM)

**Why.** Full degree 2 needs 136 terms, impossible on v2 (98 runs). Degree 2 only on the k most important
inputs and linear on the rest keeps the curvature where the variance is. Exploratory result (throwaway
scripts, not in the repo, to be reproduced):

| basis | terms | v1 RMSE | v2 RMSE |
|---|---|---|---|
| monomial deg 1 | 16 | 0.050 y | 0.062 y |
| top-4 deg 2 | 26 | 0.024 y | 0.035 y |
| top-6 deg 2 | 37 | 0.018 y | **0.026 y** (spline K = 3: 0.047 y) |
| monomial deg 2 | 136 | 0.007 y | not possible |

Terms: 1 + 15 linear + k(k+1)/2 quadratic (squares and products of the top k) → k = 6: 1 + 15 + 21 = 37.

**Top inputs** (SRC² on v1, `data/sensitivity/src_v1_969.jld2`; same order on v2):

1. `specific_heat_capacity_sandstone2` (c_p middle, 0.78)
2. `specific_heat_capacity_sandstone3` (c_p bottom, 0.10)
3. `density_sandstone2` (ρ middle, 0.081)
4. `sandstone_porosity_parameter_3` (φ middle, 0.017)
5. `density_sandstone3` (ρ bottom, 0.010)
6. `thermal_conductivity_sandstone3` (λ bottom, 0.002; λ middle is equal, any of the λ is fine)

Use the v1 ranking also for v2 (prior knowledge; v2 has only 25 test runs).

**How.** `MonomialBasis` cannot do "degree 2 on a subset", so write a small basis type, e.g. in `src/`:

```julia
struct TopKQuadraticBasis <: UncertaintyQuantification.AbstractBasis
    n::Int              # number of inputs (15)
    top::Vector{Int}    # indices of the top-k inputs in `inputs`
end
# evaluate: [1, x..., x[i]*x[j] for i ≤ j in top]
```

Check how UQ.jl evaluates a basis (`MonomialBasis` source: which method `LinearBasisFunctionModel` and
`IntervalPredictorModel` call, and `length(basis)`), and implement the same methods.
Then in `lbfm.jl`: a `BASIS` switch (`"monomial_deg1"`, `"monomial_deg2"`, `"top6_deg2"`), same file names.
The guard `length(basis) >= nrow(train)/2` stays (37 < 49 on v2: ok).

**IPM.** The same basis works in `ipm.jl`. Trade-off precision vs certainty: more terms → narrower band but
weaker guarantee (ε at 95% confidence: deg 1 → v1 5.4%, v2 40%; full deg 2 on v1 → 38%).
Check top-6 on v1 (37 terms) as a compromise; keep deg 1 on v2.

**Done when:** `lbfm_top6_deg2_<TAG>.jld2` for v1 and v2, numbers close to the table above, figures with
`figures_lbfm.jl` (works unchanged: it loops over `lbfm_*_<TAG>.jld2`).

---

## 2. Sobol indices on the surrogates (`scripts/sensitivity/sobol.jl`)

**Why.** SRC² is only first order and assumes a linear model. Sobol indices on a good surrogate give first
order and total effects (interactions), and cross-check SRC². The old PCE Sobol indices are meaningless
(flat: S1 ≈ 0.01, ST ≈ 0.2 each, overfitted PCE), do not reuse them.

**Which surrogate.** v1: LBFM deg 2 (`lbfm_monomial_deg2_v1_969_test20.jld2`).
v2: spline K = 3 (`spline_k3_v2_123_test20.jld2`), or top-6 deg 2 once step 1 is done.

**How.** The surrogates were fitted on scaled inputs, so wrap them in a `Model` that scales first:

```julia
include(joinpath(@__DIR__, "..", "surrogates", "settings.jl"))
res = load(joinpath(SURROGATES_DIR, "lbfm_monomial_deg2_$(TAG).jld2"))
lbfm, lo, hi = res["lbfm"], res["lo"], res["hi"]
inputs = input_distributions()
# lo, hi: Dicts column => training min / max (from scaler); same formula as scaler's scale(d)
scale_df(df) = (s = copy(df); for c in keys(lo)
    s[!, c] = 2 .* (df[!, c] .- lo[c]) ./ (hi[c] .- lo[c]) .- 1 end; s)
surrogate = Model(df -> (s = scale_df(df); evaluate!(lbfm, s); s[!, :crossing_year]), :crossing_year)
si = sobolindices([surrogate], inputs, :crossing_year, MonteCarlo(10^5))
```

`lo`/`hi` are `Dict`s keyed by column name (`scaler` in `settings.jl`); the `scale` closure itself is not
saved, so rebuild it from them as above. Note the Monte Carlo samples cover the whole distribution
support, slightly beyond the training min/max (scaled values a bit outside [-1, 1]: mild extrapolation).
For the spline the output column is `spline.output`
(`lbfm.out` for the LBFM). Sample many points (10⁴-10⁵): the surrogate is cheap.

**Expected.** First order ≈ SRC² (c_p middle ≈ 0.78, …), total effects slightly above (small interactions,
the deg-2 model is 99.8% linear). Save to `data/sensitivity/sobol_<TAG>.jld2`, bar chart with first order and
total effect next to SRC² in `scripts/auxiliary/figures_sobol.jl`.

---

## 3. Spline on the top-k inputs

**Why.** The spline distance treats all 15 inputs equally; the unimportant ones (κ ≈ 0, λ) only blur the
distances. A spline on the top-k inputs (k = 4-6) should be more accurate on v2. It is also the bridge to the
GP, which learns one length scale per input.

**How.** In `spline.jl`, fit `PolyharmonicSpline` on `scale(train)[:, [top...; output]]` (it uses every
column except the output as input) and evaluate the test set on the same columns. Variant name e.g.
`"k3_top6"`. Compare with K = 3 on all inputs (v1 0.009 y, v2 0.047 y).

---

## 4. Old PCE on the v1 test set

**Why.** For the comparison table, a PCE RMSE on the same test set as the other surrogates (the old PCE
results are MSE on 12 runs only).

**How.** `data/surrogates/2026_05_29_21_00_3layers_sobolsampling_512_deg_4.jld2` has key `res`
(`Vector{Any}`, the PCE is `res[1]`). Load it with `using UncertaintyQuantification` first (JLD2 prints
type-reconstruction warnings otherwise; if it still cannot rebuild the type, refit is not an option, so
check the UQ.jl version used in May 2026). The PCE was fitted on **unscaled** inputs with the v1
distributions: evaluate it on the v1_969 test set (`DATASET = "v1_969"`) without `scale`.
The PCE runs (128/256/512, dataset `v1_pce`) and v1_969 (the NAS 1024 run) are different runs, so the
v1_969 test set is unseen by the PCE.

**Expected.** Much worse than everything else (underdetermined: 3876 coefficients for 512 runs).

---

## 5. Fit times

Wrap every fit in `@elapsed` and save `fit_time` in the jld2 (ipm, lbfm, spline, and the new ones).
Needed for the comparison table. The OGS run time (> 6 h per run, 29 in parallel on the cluster) is the
reference cost.

---

## 6. Comparison table / figure

`scripts/auxiliary/compare_surrogates.jl`: read all `data/surrogates/*_<TAG>.jld2` for both datasets and
print/save one table: model, variant, terms, train N, test RMSE, max error, R², fit time; for the IPM
coverage, mean width, ε at 95% confidence. Optionally a bar chart of RMSE per model and dataset.
Optional: learning curve (test RMSE vs N on subsets of v1, 50-775 runs) for LBFM deg 1/2, top-6 and spline,
to show where each model becomes possible / better.

---

## 7. Gaussian Process (next presentation)

UQ.jl has no GP. Options: GaussianProcesses.jl, AbstractGPs.jl (+ KernelFunctions.jl, ParameterHandling.jl
for the hyperparameters). Same structure as the other surrogates (`scripts/surrogates/gp.jl`,
`figures_gp.jl`).
- kernel: squared exponential or Matérn 5/2 with **one length scale per input** (ARD), constant or linear mean
- noise term (small: OGS is deterministic, but it regularizes)
- hyperparameters by maximum marginal likelihood
- report RMSE, max error, R², and the predictive std (coverage of ±2σ on the test set, like the IPM)
- the learned length scales are a sensitivity measure: compare with SRC² (κ should get a long length scale)

---

## Model and inputs (later, see slides "Next Steps: Model and inputs")

- **Dependent inputs**: so far all 15 inputs are independent. Add correlations (deck 5 values:
  ρ(κ, φ) = 0.8, ρ(ρ, φ) = −0.8, ρ(κ, ρ) = −0.7) with a Gaussian copula (`JointDistribution` in UQ.jl);
  sensitivity then with Kucherenko indices instead of Sobol. Needs new OGS runs.
- **Flow split from OGS** (option B) instead of the fixed pumping shares (0.245 / 0.430 / 0.325): then κ could
  matter (now SRC² ≈ 0 because the split is fixed in the post-processing).
- **IPM on real / measured data**, where scatter and imprecision exist; link to the ieBN risk assessment.
