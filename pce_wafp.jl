using UncertaintyQuantification
using JLD2

include("src/model_3layers.jl")

const TRAIN_SAMPLES = 128
const PCE_DEGREE = 4

const WORK_DIR = "/work/andrea.perin/ThermoOptiPlan/output/Model_ML_IRZ_3layers"
const cleanup = false

wafp_inputs = [kappa_sandstone2, kappa_sandstone3]

ext = ogs_model(WORK_DIR; cleanup=cleanup)

models = [ext; postprocessing_models]

println("Preparing PCE with TRAIN_SAMPLES=$TRAIN_SAMPLES, degree=$PCE_DEGREE")

bases = choose_basis.(wafp_inputs)
Ψ = PolynomialChaosBasis(bases, PCE_DEGREE)

wafp = WeightedApproximateFetekePoints(MonteCarlo(5000), fadd=10, fmult=2)
# polynomialchaos(wafp_inputs, models, Ψ, :crossing_year, wafp)


wafp_inputs = UncertaintyQuantification.wrap(wafp_inputs)
outputs = UncertaintyQuantification.wrap(:simsda)

samples = sample(wafp_inputs, wafp.sim)
random_inputs = filter(i -> isa(i, RandomUQInput), wafp_inputs)
random_names = names(random_inputs)
UncertaintyQuantification.to_standard_normal_space!(random_inputs, samples)
x = UncertaintyQuantification.map_to_bases(Ψ, Matrix(samples[:, random_names]))

random_inputs = filter(i -> isa(i, RandomUQInput), wafp_inputs)


Np = length(Ψ.α)
n = wafp.sim.n
rest = min((Np - 1) * wafp.fmult, wafp.fadd, n - Np)
rest = max(rest, 0)

A = Matrix{Float64}(undef, n, Np)
for i in 1:n
    A[i, :] .= evaluate(Ψ, x[i, :])
end
w = UncertaintyQuantification.norm.(eachrow(A)) .^ (-2.0)
B = A' .* reshape(w .^ (1 / 2), 1, :)
_, _, p = UncertaintyQuantification.qr(B, UncertaintyQuantification.ColumnNorm())
pout = zeros(Int, Np + rest)
pout[1:Np] .= p[1:Np]
Ginv = inv(B[:, p] * B[:, p]')
for i in 1:rest
    val, j = findmax(j -> B[:, p[j]]' * Ginv * B[:, p[j]], Np+i:n)
    pout[Np+i] = p[j]
    Ginv .-= ((Ginv * B[:, p[j]]) * (B[:, p[j]]' * Ginv)) ./ (1 + val)
end

samples = samples[pout, :]
UncertaintyQuantification.to_physical_space!(random_inputs, samples)
w = w[pout]

# est = LeastSquares(SobolSampling(TRAIN_SAMPLES))

# println("Running polynomial chaos construction (this may run external model per sample)...")


# path_to_pce = joinpath("/work/andrea.perin/ThermoOptiPlan/results/pce")
# mkpath(path_to_pce)

# @show("start pce analysis with simulation: $(est)")
# @time pce, samples, mse = polynomialchaos(inputs, models, Ψ, :crossing_year, est)
# res = [pce, samples, mse]

# name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_3layers" * "_sobolsampling" * "_" * string(TRAIN_SAMPLES) * "_deg" * "_" * string(PCE_DEGREE) * ".jld2"
# @save joinpath(path_to_pce, name) res
