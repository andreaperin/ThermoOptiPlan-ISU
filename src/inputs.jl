# The 15 random inputs of the 3-layer model, one function per version of the distributions.
# Order and names match the dataset columns (see models/README.md for the layer naming).

using UncertaintyQuantification

"""
    inputs_v2()

Current per-layer distributions, used by run_samples.jl and for the v2 datasets (v2_128, v2_123).
"""
function inputs_v2()
    # Normal with σ = 10% of the nominal value, truncated at ±20% around it
    layer(nominal) = Truncated(Normal(nominal, 0.1 * nominal), 0.8 * nominal, 1.2 * nominal)

    λ = Truncated(Normal(2.38, 0.25), 0.75 * 2.38, 1.25 * 2.38) # ±25% window
    cp = Truncated(Normal(820, 80), 615, 1025)                 # ±25% window
    ρ = Truncated(Normal(2690, 100), 2400, 2800)               # ~±10% window

    return [
        RandomVariable(λ, :thermal_conductivity_sandstone1),
        RandomVariable(λ, :thermal_conductivity_sandstone2),
        RandomVariable(λ, :thermal_conductivity_sandstone3),
        RandomVariable(cp, :specific_heat_capacity_sandstone1),
        RandomVariable(cp, :specific_heat_capacity_sandstone2),
        RandomVariable(cp, :specific_heat_capacity_sandstone3),
        RandomVariable(ρ, :density_sandstone1),
        RandomVariable(ρ, :density_sandstone2),
        RandomVariable(ρ, :density_sandstone3),
        RandomVariable(layer(0.19), :sandstone_porosity_parameter_2),
        RandomVariable(layer(0.20), :sandstone_porosity_parameter_3),
        RandomVariable(layer(0.22), :sandstone_porosity_parameter_4),
        RandomVariable(layer(1.56e-13), :kappa_Sandstone_2),
        RandomVariable(layer(2.39e-13), :kappa_Sandstone_3),
        RandomVariable(layer(4.97e-13), :kappa_Sandstone_4),
    ]
end

"""
    inputs_v1()

Old distributions (before commit 6865b24), used for the v1 datasets (v1_969, v1_pce).
Porosity and permeability share one window for all layers, with σ much larger than the window,
so they are almost uniform; κ of the top and bottom layers has its nominal value outside the window.
"""
function inputs_v1()
    λ = Truncated(Normal(2.38, 0.25), 1.5675, 2.6125)
    cp = Truncated(Normal(820, 80), 615, 1025)
    ρ = Truncated(Normal(2690, 100), 2400, 2800)
    φ(nominal) = Truncated(Normal(nominal, 0.05), 0.18, 0.23)
    κ(nominal) = Truncated(Normal(nominal, 2.868e-13), 1.912e-13, 2.868e-13)

    return [
        RandomVariable(λ, :thermal_conductivity_sandstone1),
        RandomVariable(λ, :thermal_conductivity_sandstone2),
        RandomVariable(λ, :thermal_conductivity_sandstone3),
        RandomVariable(cp, :specific_heat_capacity_sandstone1),
        RandomVariable(cp, :specific_heat_capacity_sandstone2),
        RandomVariable(cp, :specific_heat_capacity_sandstone3),
        RandomVariable(ρ, :density_sandstone1),
        RandomVariable(ρ, :density_sandstone2),
        RandomVariable(ρ, :density_sandstone3),
        RandomVariable(φ(0.19), :sandstone_porosity_parameter_2),
        RandomVariable(φ(0.20), :sandstone_porosity_parameter_3),
        RandomVariable(φ(0.22), :sandstone_porosity_parameter_4),
        RandomVariable(κ(1.56e-13), :kappa_Sandstone_2),
        RandomVariable(κ(2.39e-13), :kappa_Sandstone_3),
        RandomVariable(κ(4.97e-13), :kappa_Sandstone_4),
    ]
end
