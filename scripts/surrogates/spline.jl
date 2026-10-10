using UncertaintyQuantification

include(joinpath(@__DIR__, "settings.jl")) # DATASET, TEST_PERCENT, TAG, train_test_split, scaler

K = 1 # polyharmonic degree, φ(r) = r^k (k odd) or r^k log(r) (k even)
VARIANT = "k$(K)" # name of the variant in the file names

train, test = train_test_split()

output = :crossing_year
inputs = Symbol.(names(train, Not(output)))

# scale into [-1,1] all the inputs (min/max of the training set): the distances r mix all inputs
scale, lo, hi = scaler(train, inputs)

# fit: linear term + one radial function per training point (uses all columns except output)
println("polyharmonic spline k = $K on $DATASET: $(nrow(train)) training points")
@time spline = PolyharmonicSpline(scale(train), K, output)

# check: the spline interpolates, so the training error should be ~0
fit = scale(train)
evaluate!(spline, fit)
println("max training error: ", maximum(abs.(fit[!, output] .- train[!, output])), " years")

# predict the test set (evaluate! overwrites the output column, so on a copy)
pred = scale(test)
evaluate!(spline, pred)
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
jldsave(joinpath(SURROGATES_DIR, "spline_$(VARIANT)_$(TAG).jld2"); spline, lo, hi, test, ŷ, rmse, maxerr, r2)
