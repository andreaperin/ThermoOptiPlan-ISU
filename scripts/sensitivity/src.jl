# Sensitivity from a linear regression on the data: squared standardized regression coefficients (SRC²).
# Valid when the linear fit explains almost all the variance (R² ≈ 1): then SRC² ≈ first-order Sobol index.

using Statistics

include(joinpath(@__DIR__, "..", "surrogates", "settings.jl")) # DATASET

# all the runs of the dataset: no test set is needed for a sensitivity analysis
df = load(joinpath(DATASETS_DIR, "dataset_3layers_$(DATASET).jld2"), "df")
output = :crossing_year
inputs = Symbol.(names(df, Not(output)))

# standardize the inputs (mean 0, std 1): the coefficients are then directly comparable,
# and the least squares problem is well conditioned (κ ~ 1e-13 next to ρ ~ 2700)
X = Matrix(df[:, inputs])
Z = (X .- mean(X; dims = 1)) ./ std(X; dims = 1)
y = df[!, output]

# linear least squares with intercept
A = [ones(nrow(df)) Z]
β = A \ y
r2 = 1 - sum((y .- A * β) .^ 2) / sum((y .- mean(y)) .^ 2)

# standardized regression coefficients: sign = direction of the effect, square = share of variance
src = β[2:end] ./ std(y)
src2 = src .^ 2

println("SRC² on $DATASET ($(nrow(df)) runs): R² of the linear fit = $(round(r2, digits = 5)), Σ SRC² = $(round(sum(src2), digits = 5))")
for i in sortperm(src2; rev = true)
    println(rpad(inputs[i], 36), "  SRC = ", lpad(round(src[i], digits = 4), 8), "   SRC² = ", round(src2[i], digits = 5))
end
mkpath(SENSITIVITY_DIR)
jldsave(joinpath(SENSITIVITY_DIR, "src_$(DATASET).jld2"); inputs, src, src2, r2)
