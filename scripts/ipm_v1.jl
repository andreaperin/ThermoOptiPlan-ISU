using UncertaintyQuantification
using JLD2, DataFrames, Random

include(joinpath(@__DIR__, "..", "src", "paths.jl"))

df = load(joinpath(DATASETS_DIR, "dataset_3layers_v1_969.jld2"), "df")

output = :crossing_year
inputs = Symbol.(names(df, Not(output)))

# train/test split: TEST_PERCENT of the runs are used for testing (e.g. 20 or 90)
TEST_PERCENT = 20
suffix = "v1_969_test$(TEST_PERCENT)" # file names: ipm_<scaled|not_scaled>_<suffix>.jld2

Random.seed!(1)
idx = shuffle(1:size(df, 1))
ntest = round(Int, TEST_PERCENT / 100 * size(df, 1))
test, train = df[idx[1:ntest], :], df[idx[(ntest + 1):end], :]

# IPM not SCALED

@time ipm_not_scaled = IntervalPredictorModel(train, output, MonomialBasis(length(inputs), 1), inputs)
pred_ns = copy(test)
@time evaluate!(ipm_not_scaled, pred_ns)
lb, ub = getproperty.(pred_ns[!, output], :lb), getproperty.(pred_ns[!, output], :ub)
y = test[!, output]
inside = lb .<= y .<= ub

println("test coverage: ", mean(inside))
println("mean width:    ", mean(ub .- lb), " years")
println("β for ≥90% coverage: ", reliability(ipm_not_scaled, 0.1))

mkpath(SURROGATES_DIR)
jldsave(joinpath(SURROGATES_DIR, "ipm_not_scaled_$(suffix).jld2"); ipm_not_scaled, test, lb, ub)

# IPM scaled

# scale into [-1,1] all the inputs
lo = Dict(c => minimum(train[!, c]) for c in inputs)
hi = Dict(c => maximum(train[!, c]) for c in inputs)
function scale(d)
    s = copy(d)
    for c in inputs
        s[!, c] = 2 .* (d[!, c] .- lo[c]) ./ (hi[c] .- lo[c]) .- 1
    end
    return s
end

@time ipm_scaled = IntervalPredictorModel(scale(train), output, MonomialBasis(length(inputs), 1), inputs)

pred = scale(test)
@time evaluate!(ipm_scaled, pred)
lb, ub = getproperty.(pred[!, output], :lb), getproperty.(pred[!, output], :ub)
y = test[!, output]
inside = lb .<= y .<= ub

println("test coverage: ", mean(inside))
println("mean width:    ", mean(ub .- lb), " years")
println("β for ≥90% coverage: ", reliability(ipm_scaled, 0.1))

mkpath(SURROGATES_DIR)
jldsave(joinpath(SURROGATES_DIR, "ipm_scaled_$(suffix).jld2"); ipm_scaled, lo, hi, test, lb, ub)
