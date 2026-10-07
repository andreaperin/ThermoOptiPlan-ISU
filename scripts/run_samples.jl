# Run the 3-layer OGS model on N Sobol samples (on SLURM) and save inputs + outputs for GP training.

using Dates
using JLD2
using DataFrames
using UncertaintyQuantification

include(joinpath(@__DIR__, "..", "src", "model_3layers.jl"))

const N_SAMPLES = 128
const cleanup = true   # false keeps the raw .vtu files (~1.7 GB per run)

const WORK_DIR = joinpath(RUNS_DIR, "multilayer_irz")

options = Dict(
    "job-name" => "run_samples",
    "account" => "andrea.perin",
    "ntasks" => "1",
    "cpus-per-task" => "1",
    "mem-per-cpu" => "2G"
)

slurm = SlurmInterface(
    options;
    throttle = 29,
    extras = ["module load OpenGeoSys", "export OMP_NUM_THREADS=1"],
)

ext = ogs_model(WORK_DIR; cleanup = cleanup, scheduler = slurm)
models = [ext; postprocessing_models]

samples = sample(inputs, SobolSampling(N_SAMPLES))
@time evaluate!(models, samples)

# GP table: 15 inputs + crossing_year, crashed runs removed
input_names = [string(rv.name) for rv in inputs]
gp_data = samples[.!isnan.(samples.crossing_year), [input_names; "crossing_year"]]

mkpath(DATASETS_DIR)
name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_3layers_sobol_" * string(N_SAMPLES) * ".jld2"
jldsave(joinpath(DATASETS_DIR, name); samples, gp_data)
println(size(gp_data, 1), " / ", N_SAMPLES, " complete runs saved to ", joinpath(DATASETS_DIR, name))
