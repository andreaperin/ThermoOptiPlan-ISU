# All paths of the project. Generated data lives under DATA_DIR, which defaults to
# <repo>/data and can be moved (e.g. to /work on the cluster) with THERMOOPTIPLAN_DATA.

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const MODELS_DIR = joinpath(ROOT_DIR, "models")

const DATA_DIR = get(ENV, "THERMOOPTIPLAN_DATA", joinpath(ROOT_DIR, "data"))
const RUNS_DIR = joinpath(DATA_DIR, "runs")             # raw OGS output
const SURROGATES_DIR = joinpath(DATA_DIR, "surrogates") # fitted PCE/GP objects
const DATASETS_DIR = joinpath(DATA_DIR, "datasets")     # clean input/output tables
