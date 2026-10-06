using UncertaintyQuantification
using JLD2
using LaTeXStrings
using PGFPlotsX
using DataFrames

# List of all rvs with their .jld2 name and their latexify version
function latexify(varname::String)
    # if varname == "thermal_conductivity_sandstone_3"
    #     return L"Th Cond"
    # elseif varname == "specific_heat_capacity_sandstone_3"
    #     return L"Spec Heat"
    # elseif varname == "density_sandstone_3"
    #     return L"Density"
    # elseif varname == "kappa_sandstone_3"
    #     return L"Permeability"
    # elseif varname == "sandstone_porosity_parameter_3"
    #     return L"Porosity"
    if varname == "thermal_conductivity_sandstone1"
        return L"\kappa s1"
    elseif varname == "thermal_conductivity_sandstone2"
        return L"\kappa s2"
    elseif varname == "thermal_conductivity_sandstone3"
        return L"\kappa s3"
    elseif varname == "specific_heat_capacity_sandstone1"
        return L"c 1"
    elseif varname == "specific_heat_capacity_sandstone2"
        return L"c 2"
    elseif varname == "specific_heat_capacity_sandstone3"
        return L"c 3"
    elseif varname == "density_sandstone1"
        return L"\rho 1"
    elseif varname == "density_sandstone2"
        return L"\rho 2"
    elseif varname == "density_sandstone3"
        return L"\rho 3"
    elseif varname == "sandstone_porosity_parameter_2"
        return L"\phi 1"
    elseif varname == "sandstone_porosity_parameter_3"
        return L"\phi 2"
    elseif varname == "sandstone_porosity_parameter_4"
        return L"\phi 3"
    elseif varname == "kappa_Sandstone_2"
        return L"k 1"
    elseif varname == "kappa_Sandstone_3"
        return L"k 2"
    elseif varname == "kappa_Sandstone_4"
        return L"k 3"
    end

end

function plot_sobols(
    df::DataFrame,
    title::LaTeXString;
    height::String="8cm",
    width::String="15cm",
    bar_width::Float64=0.2,
    bar_shift_amt::Float64=0.20,
    color_first_order::String="orange",
    color_total_effect::String="red",
    ymin::Real=0,
    ymax::Real=1
)

    class_labels = string.(df.Variables)
    class_labels_latex = latexify.(class_labels)

    group1_values = df.FirstOrder
    group2_values = df.TotalEffect
    n = length(class_labels)

    pgf = @pgf Axis(
        {
            ybar,
            bar_width = bar_width,
            xtick = 1:n,
            xticklabels = class_labels_latex,
            xticklabel_style = "{rotate=45, anchor=east}",  # Better for long labels
            enlargelimits = 0.05,
            ylabel = "",
            ymin = ymin,
            ymax = ymax,
            height = height,
            width = width,
            title = title,
            grid = "major"
        },
        # Group 1: shifted left
        Plot({fill = color_first_order, bar_shift = "-$(bar_shift_amt)cm"}, Coordinates([(i, group1_values[i]) for i in 1:n])),
        # Group 2: shifted right
        Plot({fill = color_total_effect, bar_shift = "$(bar_shift_amt)cm"}, Coordinates([(i, group2_values[i]) for i in 1:n])),
        Legend(["First Order", "Total Effect"])
    )
    return pgf
end