# Plots of the splines fitted by scripts/surrogates/spline.jl (run that first).
# Dataset and test split come from scripts/surrogates/settings.jl, as for the fit.

using UncertaintyQuantification
using PGFPlotsX, LaTeXStrings

include(joinpath(@__DIR__, "..", "surrogates", "settings.jl")) # DATASET, TEST_PERCENT, TAG

"""
    plot_spline(file; name)

Parity and error plots of the spline saved in `file` (with `test`, `ŷ`, `rmse`, `r2`).
`name` goes in the titles. The PDFs are saved in FIGURES_DIR as `<file name>_<plot>.pdf`.
"""
function plot_spline(file::String; name::String)
    res = load(file)
    spline, test, ŷ = res["spline"], res["test"], res["ŷ"]
    prefix = joinpath(FIGURES_DIR, splitext(basename(file))[1])
    mkpath(FIGURES_DIR)

    y = test[!, spline.output]
    err = ŷ .- y

    # 1. Parity plot: OGS value vs spline prediction
    ymin, ymax = minimum([ŷ; y]) - 0.2, maximum([ŷ; y]) + 0.2
    p_parity = @pgf Axis(
        {
            title = "spline $name - RMSE $(round(res["rmse"], digits = 3)) y, \$R^2\$ $(round(res["r2"], digits = 5))",
            xlabel = "crossing year OGS [y]",
            ylabel = "crossing year spline [y]",
            grid = "major",
            xmin = ymin, xmax = ymax, ymin = ymin, ymax = ymax,
        },
        Plot({only_marks, mark_size = "1pt", color = "blue"}, Coordinates(y, ŷ)),
        Plot({dashed, color = "black"}, Coordinates([(ymin, ymin), (ymax, ymax)])),
    )
    PGFPlotsX.save(prefix * "_parity.pdf", p_parity)

    # 2. Error (spline - OGS), test points sorted by OGS value (same x axis as the IPM band plot)
    order = sortperm(y)
    k = 1:length(y)
    p_err = @pgf Axis(
        {
            title = "spline $name: prediction minus OGS value",
            xlabel = "test sample (sorted by OGS value)",
            ylabel = "spline - OGS [y]",
            grid = "major",
            legend_pos = "outer north east",
        },
        Plot({only_marks, mark_size = "1pt", color = "blue"}, Coordinates(k, err[order])),
        LegendEntry("test points"),
        Plot({thick, color = "black"}, Coordinates([1, length(y)], [0, 0])),
        LegendEntry("OGS value"),
        Plot({dashed, color = "red"}, Coordinates([1, length(y)], [res["rmse"], res["rmse"]])),
        LegendEntry(L"\pm \mathrm{RMSE}"),
        Plot({dashed, color = "red", forget_plot}, Coordinates([1, length(y)], [-res["rmse"], -res["rmse"]])),
    )
    PGFPlotsX.save(prefix * "_error.pdf", p_err)

    return println("saved ", prefix, "_{parity,error}.pdf")
end

for file in filter(f -> startswith(f, "spline_") && endswith(f, "_$(TAG).jld2"), readdir(SURROGATES_DIR))
    basis = replace(file, "spline_" => "", "_$(TAG).jld2" => "")
    plot_spline(joinpath(SURROGATES_DIR, file); name = "$(replace(basis, "_" => " ")), $(replace(DATASET, "_" => " ")), test $(TEST_PERCENT)\\%")
end
