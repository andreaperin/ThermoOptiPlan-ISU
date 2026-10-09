# Settings shared by ipm.jl, lbfm.jl and the figures in scripts/auxiliary:
# change them here, rerun the fit, then the figures.

using JLD2, DataFrames, Random

include(joinpath(@__DIR__, "..", "..", "src", "paths.jl"))

DATASET = "v2_123" # uses data/datasets/dataset_3layers_<DATASET>.jld2
TEST_PERCENT = 20  # % of the runs kept for testing (e.g. 20 or 90)
SEED = 1           # same seed -> same train/test split for every surrogate

TAG = "$(DATASET)_test$(TEST_PERCENT)" # in the names of the saved surrogates and figures

"""
    train_test_split()

Load the dataset `DATASET` and return `(train, test)`, with `TEST_PERCENT` % of the runs in `test`.
"""
function train_test_split()
    df = load(joinpath(DATASETS_DIR, "dataset_3layers_$(DATASET).jld2"), "df")
    Random.seed!(SEED)
    idx = shuffle(1:size(df, 1))
    ntest = round(Int, TEST_PERCENT / 100 * size(df, 1))
    return df[idx[(ntest + 1):end], :], df[idx[1:ntest], :]
end
