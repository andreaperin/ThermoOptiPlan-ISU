# Inputs + crossing_year of all OGS runs in the PCE result files, for training a GP.

include(joinpath(@__DIR__, "..", "..", "tools", "dataset.jl"))

files = filter(endswith(".jld2"), readdir(SURROGATES_DIR; join=true))
df = dataframe_from_pce(files)

mkpath(DATASETS_DIR)
jldsave(joinpath(DATASETS_DIR, "dataset_3layers.jld2"); df)
println(size(df, 1), " runs saved to ", joinpath(DATASETS_DIR, "dataset_3layers.jld2"))
