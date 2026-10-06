using UncertaintyQuantification

include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(@__DIR__, "extractor.jl"))

# Layer naming: the solid-property RVs are numbered 1/2/3, the OGS parameters 2/3/4.
#   sandstone1 / Sandstone_2 -> top    (Sand_2, medium 1)
#   sandstone2 / Sandstone_3 -> middle (Sand_3, medium 2)
#   sandstone3 / Sandstone_4 -> bottom (Sand_4, medium 3)

const OGS_CMD = "ogs"
const SOURCE_DIR = joinpath(MODELS_DIR, "multilayer_irz")

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
const extrafiles = [
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

# Random variables

# Normal with σ = 10% of the nominal value, truncated at ±20% around it
layer_distribution(nominal::Real) = Truncated(Normal(nominal, 0.1 * nominal), 0.8 * nominal, 1.2 * nominal)

dist_thermal_conductivity_sandstone = Truncated(Normal(2.38, 0.25), 0.75 * 2.38, 1.25 * 2.38) # Approx 25% uncertainty window tapered towards the limits
dist_specific_heat_capacity_sandstone = Truncated(Normal(820, 80), 615, 1025) # Approx 25% uncertainty window tapered towards the limits
dist_density_sandstone = Truncated(Normal(2690, 100), 2400, 2800) # Approx 10% uncertainty window tapered towards the limits

dist_porosity_parameter_sandstone1 = layer_distribution(0.19)
dist_porosity_parameter_sandstone2 = layer_distribution(0.20)
dist_porosity_parameter_sandstone3 = layer_distribution(0.22)

dist_kappa_sandstone1 = layer_distribution(1.56e-13)
dist_kappa_sandstone2 = layer_distribution(2.39e-13)
dist_kappa_sandstone3 = layer_distribution(4.97e-13)

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

# OGS model

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

function ogs_model(workdir::String; cleanup::Bool=false, scheduler=nothing)
    return ExternalModel(
        SOURCE_DIR,
        sourcefile,
        extractor,
        ogs;
        workdir=workdir,
        extras=extrafiles,
        cleanup=cleanup,
        scheduler=scheduler
    )
end

# Post-processing models

# Pumping rates per metre of well line, as in the .prj (Pump_4, Pump_3, Pump_2).
# The pump lines span exactly the Δz intervals, so the flow of each layer is rate * thickness.
const pump_rates = [7.966e-4, 3.72e-4, 2.43e-4] # bottom, middle, top

function flow_percentage(rates::Vector{<:Real}, Δz_bottom::Tuple, Δz_middle::Tuple, Δz_top::Tuple)
    thicknesses = [abs(a - b) for (a, b) in [Δz_bottom, Δz_middle, Δz_top]]
    flows = rates .* thicknesses
    return flows ./ sum(flows)
end

# Fixed shares (≈ 0.325, 0.430, 0.245), consistent with the rates OGS simulates.
# Future work: make the pumping rates κ-dependent in the .prj (fixed total, split by κ·h)
# and compute these shares from the sampled permeabilities instead.
const flow_shares = flow_percentage(pump_rates, Δz_bottom, Δz_middle, Δz_top)

flow_model = Model(df -> fill(flow_shares, size(df, 1)), :flows)

function final_temperature(extraction_temperatures::Vector{Vector{Float64}}, flows::Vector{Float64})
    return [[sum(v[1:3] .* flows), v[4]] for v in extraction_temperatures]
end

final_T_model = Model(df -> final_temperature.(df.extraction_temperatures, df.flows), :final_T)

# End of the simulation (t_end in the .prj) in years
const t_end_years = 2365200000 / 365 / 24 / 60 / 60

function crossing_year(final_temperature_data::Vector{Vector{Float64}})
    T = first.(final_temperature_data)
    t = last.(final_temperature_data)

    # incomplete (crashed) OGS run
    if t[end] < t_end_years - 1e-6
        return NaN
    end

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

postprocessing_models = [flow_model, final_T_model, crossing_year_model]

# PCE basis

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
