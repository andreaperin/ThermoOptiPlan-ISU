using UncertaintyQuantification

include(joinpath(@__DIR__, "settings.jl")) # DATASET, TEST_PERCENT, TAG, train_test_split
DEGREE = 2
BASIS = "monomial_deg$(DEGREE)" # name of the basis in the file names (other bases later)

train, test = train_test_split()

output = :crossing_year
inputs = Symbol.(names(train, Not(output)))

scale, lo, hi = scaler(train, inputs)

basis = MonomialBasis(length(inputs), DEGREE)
# add check for degree vs training size
if length(basis) >= nrow(train) / 2
    error("degree $DEGREE has $(length(basis)) terms for $(nrow(train)) training points: too few points")
end
println("LBFM degree $DEGREE on $DATASET: $(length(basis)) terms, $(nrow(train)) training points")

# fit: least squares on the scaled training set
@time lbfm = LinearBasisFunctionModel(scale(train), output, basis, inputs)

# predict the test set (evaluate! overwrites the output column, so on a copy)
pred = scale(test)
evaluate!(lbfm, pred)
ŷ = pred[!, output]
y = test[!, output]

# errors on the test set [years]
err = ŷ .- y
rmse = sqrt(mean(err .^ 2))
maxerr = maximum(abs.(err))
r2 = 1 - sum(err .^ 2) / sum((y .- mean(y)) .^ 2)

println("test RMSE:      ", rmse, " years")
println("test max error: ", maxerr, " years")
println("test R²:        ", r2)

mkpath(SURROGATES_DIR)
jldsave(joinpath(SURROGATES_DIR, "lbfm_$(BASIS)_$(TAG).jld2"); lbfm, lo, hi, test, ŷ, rmse, maxerr, r2)
