using Dates
using UncertaintyQuantification
using JLD2

include("extractor.jl")

const TRAIN_SAMPLES = 1024
const PCE_DEGREE = 4

const OGS_CMD = "ogs"
const WORK_DIR = "/work/andrea.perin/ThermoOptiPlan/output/Model_ML_IRZ_3layers"
const SOURCE_DIR = joinpath(pwd(), "model_inputs", "Model_ML_IRZ")
const cleanup = true

const path_to_pce = joinpath("/work/andrea.perin/ThermoOptiPlan/results/pce")

const x_extractor = 2_250.0
const y_extractor = 0.0

const Δz_bottom = (-1374.8, -1364.6)
const Δz_middle = (-1338.0, -1309.1)
const Δz_top = (-1265.4, -1240.2)

const Δz_extractor = [
                      Δz_bottom,
                      Δz_middle,
                      Δz_top
                     ]

const sourcefile = "MULTI_BW_line_IRZ.prj"
extrafiles = [
              "Multi_BW_line_IRZ_domain_ini.vtu",
              "Multi_BW_line_IRZ_domain.vtu",
              "Multi_BW_line_IRZ.geo",
              "Multi_BW_line_IRZ.msh",
              "Multi_BW_line_IRZ_physical_group_boundary_ini.vtu",
              "Multi_BW_line_IRZ_physical_group_boundary.vtu",
              "Multi_BW_line_IRZ_physical_group_Clay_1.vtu",
              "Multi_BW_line_IRZ_physical_group_Clay_2.vtu",
              "Multi_BW_line_IRZ_physical_group_Inj_line_2.vtu",
              "Multi_BW_line_IRZ_physical_group_Inj_line_3.vtu",
              "Multi_BW_line_IRZ_physical_group_Inj_line_4.vtu",
              "Multi_BW_line_IRZ_physical_group_Pump_line_2.vtu",
              "Multi_BW_line_IRZ_physical_group_Pump_line_3.vtu",
              "Multi_BW_line_IRZ_physical_group_Pump_line_4.vtu",
              "Multi_BW_line_IRZ_physical_group_Sand_2.vtu",
              "Multi_BW_line_IRZ_physical_group_Sand_3.vtu",
              "Multi_BW_line_IRZ_physical_group_Sand_4.vtu",
              "Multi_BW_line_IRZ.pvd"
             ]

dist_thermal_conductivity_sandstone = Truncated(Normal(2.38, 0.25), 1.5675, 2.6125) # Approx 25% uncertainty window tapered towards the limits
dist_specific_heat_capacity_sandstone = Truncated(Normal(820, 80), 615, 1025) # Approx 25% uncertainty window tapered towards the limits
dist_density_sandstone = Truncated(Normal(2690, 100), 2400, 2800) # Approx 10% uncertainty window tapered towards the limits

dist_porosity_parameter_sandstone1 = Truncated(Normal(0.19, 0.05), 0.18, 0.23)
dist_porosity_parameter_sandstone2 = Truncated(Normal(0.20, 0.05), 0.18, 0.23)
dist_porosity_parameter_sandstone3 = Truncated(Normal(0.22, 0.05), 0.18, 0.23)

dist_kappa_sandstone1 = Truncated(Normal(1.56e-13, 2.868e-13), 1.912e-13, 2.868e-13)
dist_kappa_sandstone2 = Truncated(Normal(2.39e-13, 2.868e-13), 1.912e-13, 2.868e-13)
dist_kappa_sandstone3 = Truncated(Normal(4.97e-13, 2.868e-13), 1.912e-13, 2.868e-13)

thermal_conductivity_sandstone1 = RandomVariable(dist_thermal_conductivity_sandstone, :thermal_conductivity_sandstone1)
thermal_conductivity_sandstone2 = RandomVariable(dist_thermal_conductivity_sandstone, :thermal_conductivity_sandstone2)
thermal_conductivity_sandstone3 = RandomVariable(dist_thermal_conductivity_sandstone, :thermal_conductivity_sandstone3)

specific_heat_capacity_sandstone1 = RandomVariable(dist_specific_heat_capacity_sandstone, :specific_heat_capacity_sandstone1)
specific_heat_capacity_sandstone2 = RandomVariable(dist_specific_heat_capacity_sandstone, :specific_heat_capacity_sandstone2)
specific_heat_capacity_sandstone3 = RandomVariable(dist_specific_heat_capacity_sandstone, :specific_heat_capacity_sandstone3)

density_sandstone1 = RandomVariable(dist_density_sandstone, :density_sandstone1)
density_sandstone2 = RandomVariable(dist_density_sandstone, :density_sandstone2)
density_sandstone3 = RandomVariable(dist_density_sandstone, :density_sandstone3)

porosity_parameter1 = RandomVariable(dist_porosity_parameter_sandstone1, :sandstone_porosity_parameter_2)
porosity_parameter2 = RandomVariable(dist_porosity_parameter_sandstone2, :sandstone_porosity_parameter_3)
porosity_parameter3 = RandomVariable(dist_porosity_parameter_sandstone3, :sandstone_porosity_parameter_4)

kappa_sandstone1 = RandomVariable(dist_kappa_sandstone1, :kappa_Sandstone_2)
kappa_sandstone2 = RandomVariable(dist_kappa_sandstone2, :kappa_Sandstone_3)
kappa_sandstone3 = RandomVariable(dist_kappa_sandstone3, :kappa_Sandstone_4)

inputs = [thermal_conductivity_sandstone1, thermal_conductivity_sandstone2, thermal_conductivity_sandstone3, specific_heat_capacity_sandstone1, specific_heat_capacity_sandstone2, specific_heat_capacity_sandstone3, density_sandstone1, density_sandstone2, density_sandstone3, porosity_parameter1, porosity_parameter2, porosity_parameter3, kappa_sandstone1, kappa_sandstone2, kappa_sandstone3]

ogs = Solver(OGS_CMD,
             sourcefile;
             args="",
            )

extractor = Extractor(base -> begin
                          x = x_extractor
                          y = y_extractor
                          Δz = Δz_extractor

                          return extraction_temperatures_over_time(base, x, y, Δz)
                      end, :extraction_temperatures)

options = Dict(
               "job-name" => "pce_internal_scheduler",
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

ext = ExternalModel(
                    SOURCE_DIR,
                    sourcefile,
                    extractor,
                    ogs,
                    workdir=WORK_DIR,
                    extras=extrafiles,
                    cleanup=cleanup,
                    scheduler=slurm
                   )

function flow_percentage(kappa_bottom::Real, kappa_middle::Real, kappa_top::Real, Δz_bottom::Tuple, Δz_middle::Tuple, Δz_top::Tuple)
    kappas = [kappa_bottom, kappa_middle, kappa_top]
    thicknesses = [abs(a - b) for (a, b) in [Δz_bottom, Δz_middle, Δz_top]]
    flows = kappas .* thicknesses
    return flows ./ sum(flows)
end

flow_model = Model(df -> flow_percentage.(df.kappa_Sandstone_4, df.kappa_Sandstone_3, df.kappa_Sandstone_2, Ref(Δz_bottom), Ref(Δz_middle), Ref(Δz_top)), :flows)

function final_temperature(extraction_temperatures::Vector{Vector{Float64}}, flows::Vector{Float64})
    return [[sum(v[1:3] .* flows), v[4]] for v in extraction_temperatures]
end

final_T_model = Model(df -> final_temperature.(df.extraction_temperatures, df.flows), :final_T)

function crossing_year(final_temperature_data::Vector{Vector{Float64}})
    T = first.(final_temperature_data)
    t = last.(final_temperature_data)

    T0 = T[1]
    target = T0 - 1

    idx = findfirst(T .<= target)
    if idx === nothing || idx == 1
        return t[end]
    end
    # points around the crossing
    t1, t2 = t[idx-1], t[idx]
    T1, T2 = T[idx-1], T[idx]

    # linear interpolation
    return t1 + (target - T1) * (t2 - t1) / (T2 - T1)
end

crossing_year_model = Model(df -> crossing_year.(df.final_T), :crossing_year)

models = [ext, flow_model, final_T_model, crossing_year_model]

println("Preparing PCE with TRAIN_SAMPLES=$TRAIN_SAMPLES, degree=$PCE_DEGREE")

function choose_basis(rv)
    try
        d = rv.dist
        if d isa Uniform
            return LegendreBasis()
        else
            return HermiteBasis()
        end
    catch
        return HermiteBasis()
    end
end

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
