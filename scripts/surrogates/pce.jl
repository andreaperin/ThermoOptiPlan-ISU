using Dates
using UncertaintyQuantification
using JLD2

include(joinpath(@__DIR__, "..", "..", "src", "model_3layers.jl"))

const TRAIN_SAMPLES = 1024
const PCE_DEGREE = 4

const WORK_DIR = joinpath(RUNS_DIR, "multilayer_irz")
const cleanup = true

const path_to_pce = SURROGATES_DIR

options = Dict(
               "job-name" => "pce_sobol",
               "account" => "andrea.perin",
               "ntasks" => "1",
               "cpus-per-task" => "1",
               "mem-per-cpu" => "2G"
              )

slurm = SlurmInterface(
                       options;
                       throttle=29,
                       extras=["module load OpenGeoSys", "export OMP_NUM_THREADS=1"],
                      )

ext = ogs_model(WORK_DIR; cleanup=cleanup, scheduler=slurm)

models = [ext; postprocessing_models]

println("Preparing PCE with TRAIN_SAMPLES=$TRAIN_SAMPLES, degree=$PCE_DEGREE")

bases = choose_basis.(inputs)
Ψ = PolynomialChaosBasis(bases, PCE_DEGREE)

est = LeastSquares(SobolSampling(TRAIN_SAMPLES))

println("Running polynomial chaos construction (this may run external model per sample)...")

mkpath(path_to_pce)

@show("start pce analysis with simulation: $(est)")
@time pce, samples, mse = polynomialchaos(inputs, models, Ψ, :crossing_year, est)
res = [pce, samples, mse]

name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_3layers" * "_sobolsampling" * "_" * string(TRAIN_SAMPLES) * "_deg" * "_" * string(PCE_DEGREE) * ".jld2"
@save joinpath(path_to_pce, name) res
