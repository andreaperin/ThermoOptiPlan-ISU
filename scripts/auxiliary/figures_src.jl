# Bar chart of the SRC² sensitivity indices computed by scripts/sensitivity/src.jl, v1 and v2 side by side.

using PGFPlotsX, LaTeXStrings

include(joinpath(@__DIR__, "..", "surrogates", "settings.jl")) # SENSITIVITY_DIR, FIGURES_DIR

# LaTeX label of an input column: symbol with the layer as subscript
function input_label(name::Symbol)
    s = string(name)
    layer = Dict('1' => "top", '2' => "mid", '3' => "bot")
    ogs_layer = Dict('2' => "top", '3' => "mid", '4' => "bot") # OGS parameters are numbered 2/3/4
    startswith(s, "thermal_conductivity") && return "\$\\lambda_{\\mathrm{$(layer[s[end]])}}\$"
    startswith(s, "specific_heat_capacity") && return "\$c_{p,\\mathrm{$(layer[s[end]])}}\$"
    startswith(s, "density") && return "\$\\rho_{\\mathrm{$(layer[s[end]])}}\$"
    startswith(s, "sandstone_porosity") && return "\$\\phi_{\\mathrm{$(ogs_layer[s[end]])}}\$"
    startswith(s, "kappa") && return "\$\\kappa_{\\mathrm{$(ogs_layer[s[end]])}}\$"
    return s
end

v1 = load(joinpath(SENSITIVITY_DIR, "src_v1_969.jld2"))
v2 = load(joinpath(SENSITIVITY_DIR, "src_v2_123.jld2"))

# inputs sorted by SRC² on v1, largest on the left
order = sortperm(v1["src2"]; rev = true)
k = collect(1:length(order))

p = @pgf Axis(
    {
        ybar,
        bar_width = "4pt",
        title = L"Sensitivity of the crossing year: $\mathrm{SRC}^2$ of a linear fit",
        ylabel = L"$\mathrm{SRC}^2$",
        xtick = k,
        xticklabels = input_label.(v1["inputs"][order]),
        ymin = 0,
        enlarge_x_limits = 0.04,
        width = "16cm",
        height = "8cm",
        grid = "major",
        legend_pos = "north east",
    },
    Plot({fill = "blue", draw = "none", area_legend}, Coordinates(k, v1["src2"][order])),
    LegendEntry("v1 969 (\$R^2 = $(round(v1["r2"], digits = 4))\$)"),
    Plot({fill = "orange", draw = "none", area_legend}, Coordinates(k, v2["src2"][order])),
    LegendEntry("v2 123 (\$R^2 = $(round(v2["r2"], digits = 4))\$)"),
)

mkpath(FIGURES_DIR)
file = joinpath(FIGURES_DIR, "sensitivity_src.pdf")
PGFPlotsX.save(file, p)
println("saved ", file)
