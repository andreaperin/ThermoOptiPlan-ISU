# Plots of the IPMs fitted by scripts/ipm_v1.jl (run that first, with TEST_PERCENT = 20 and/or 90).

using UncertaintyQuantification
using JLD2, DataFrames
using PGFPlotsX, LaTeXStrings

include(joinpath(@__DIR__, "..", "..", "src", "paths.jl"))

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")

"""
    plot_ipm(file, key; name=key)

Parity, band and reliability plots of the IPM saved in `file` under `key` (with `test`, `lb`, `ub`).
`name` goes in the titles. The PDFs are saved in FIGURES_DIR as `<file name>_<plot>.pdf`.
"""
function plot_ipm(file::String, key::String; name::String = key)
    res = load(file)
    ipm, test, lb, ub = res[key], res["test"], res["lb"], res["ub"]
    prefix = joinpath(FIGURES_DIR, splitext(basename(file))[1])
    mkpath(FIGURES_DIR)

    y = test[!, ipm.out]
    inside = lb .<= y .<= ub

    # 1. Parity plot: OGS value vs IPM interval, test points outside the interval in red
    mid = (lb .+ ub) ./ 2
    ymin, ymax = minimum([lb; y]) - 0.2, maximum([ub; y]) + 0.2
    p_parity = @pgf Axis(
        {
            title = "IPM $name - coverage $(round(100 * mean(inside), digits = 1))\\%",
            xlabel = "crossing year OGS [y]",
            ylabel = "crossing year IPM [y]",
            grid = "major",
            legend_pos = "north west",
            xmin = ymin, xmax = ymax, ymin = ymin, ymax = ymax,
        },
        Plot(
            {only_marks, mark_size = "1pt", color = "blue", "error bars/y dir" = "both", "error bars/y explicit"},
            Coordinates(y[inside], mid[inside]; yerrorplus = (ub .- mid)[inside], yerrorminus = (mid .- lb)[inside]),
        ),
        LegendEntry("inside"),
        Plot(
            {only_marks, mark_size = "1pt", color = "red", "error bars/y dir" = "both", "error bars/y explicit"},
            Coordinates(y[.!inside], mid[.!inside]; yerrorplus = (ub .- mid)[.!inside], yerrorminus = (mid .- lb)[.!inside]),
        ),
        LegendEntry("outside"),
        Plot({dashed, color = "black"}, Coordinates([(ymin, ymin), (ymax, ymax)])),
    )
    PGFPlotsX.save(prefix * "_parity.pdf", p_parity)

    # 2. IPM interval relative to the OGS value, test points sorted by OGS value:
    #    the OGS value is inside the interval where the band crosses 0
    order = sortperm(y)
    k = 1:length(y)
    p_band = @pgf Axis(
        {
            title = "IPM $name: interval minus OGS value",
            xlabel = "test sample (sorted by OGS value)",
            ylabel = "IPM bound - OGS [y]",
            grid = "major",
            legend_pos = "north west",
        },
        Plot({name_path = "lb", draw = "none", forget_plot}, Coordinates(k, (lb .- y)[order])),
        Plot({name_path = "ub", draw = "none", forget_plot}, Coordinates(k, (ub .- y)[order])),
        Plot({fill = "blue", fill_opacity = 0.3}, raw"fill between [of=lb and ub]"),
        LegendEntry("IPM interval"),
        Plot({thick, color = "black"}, Coordinates([1, length(y)], [0, 0])),
        LegendEntry("OGS value"),
        Plot({only_marks, mark = "x", color = "red"}, Coordinates(findall(.!inside[order]), zeros(count(.!inside)))),
        LegendEntry("outside"),
    )
    PGFPlotsX.save(prefix * "_band.pdf", p_band)

    # 3. Reliability: confidence 1-β that at least 1-ϵ of new data falls in the interval
    ϵs = range(0.01, 0.2; length = 100)
    p_rel = @pgf Axis(
        {
            title = "IPM $name (N = $(ipm.N), $(length(ipm.b)) terms)",
            xlabel = L"\epsilon",
            ylabel = L"confidence $1-\beta$",
            grid = "major",
            xtick = [0.0, 0.05, 0.1, 0.15, 0.2],
            xticklabels = ["0", "0.05", "0.1", "0.15", "0.2"],
            legend_pos = "south east",
            ymin = 0, ymax = 1.05,
        },
        Plot({thick, color = "blue"}, Coordinates(ϵs, 1 .- reliability.(Ref(ipm), ϵs))),
        LegendEntry(L"P(\mathrm{coverage} \geq 1-\epsilon)"),
        Plot({dashed, color = "red"}, Coordinates([1 - mean(inside), 1 - mean(inside)], [0, 1.05])),
        LegendEntry("test miss rate"),
    )
    PGFPlotsX.save(prefix * "_reliability.pdf", p_rel)

    return println("saved ", prefix, "_{parity,band,reliability}.pdf")
end

for test_percent in (20, 90), (key, name) in [("ipm_scaled", "scaled"), ("ipm_not_scaled", "not scaled")]
    file = joinpath(SURROGATES_DIR, "$(key)_v1_969_test$(test_percent).jld2")
    if !isfile(file)
        println("missing ", basename(file), ": run scripts/ipm_v1.jl with TEST_PERCENT = $test_percent")
        continue
    end
    plot_ipm(file, key; name = "$name, test $test_percent\\%")
end
