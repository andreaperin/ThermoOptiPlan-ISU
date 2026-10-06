include(joinpath(@__DIR__, "..", "src", "paths.jl"))
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))

# Sobols

pce128 = load(joinpath(SURROGATES_DIR, "2026_05_01_04_09_3layers_sobolsampling_128_deg_4.jld2"), "res")[1]

title = L"Sobol's Indices - Crossing Year - 128 samples - deg 4"
plot_sobol = plot_sobols(
    sobolindices(pce128),
    title,
    bar_width=0.3,
    ymin=0.0,
    ymax=0.4
)
path2save = "/home/perin/Documents/academic/slides/ThermoOptiPlan/5_presentation_INTERNAL/imgs"
PGFPlotsX.save(joinpath(path2save, "Sobols128.pdf"), plot_sobol)

pce256 = load(joinpath(SURROGATES_DIR, "2026_05_17_13_15_3layers_sobolsampling_256_deg_4.jld2"), "res")[1]

title = L"Sobol's Indices - Crossing Year - 256 samples - deg 4"
plot_sobol = plot_sobols(
    sobolindices(pce256),
    title,
    bar_width=0.3,
    ymin=0.0,
    ymax=0.4
)
path2save = "/home/perin/Documents/academic/slides/ThermoOptiPlan/5_presentation_INTERNAL/imgs"
PGFPlotsX.save(joinpath(path2save, "Sobols256.pdf"), plot_sobol)

pce512 = load(joinpath(SURROGATES_DIR, "2026_05_29_21_00_3layers_sobolsampling_512_deg_4.jld2"), "res")[1]

title = L"Sobol's Indices - Crossing Year - 512 samples - deg 4"
plot_sobol = plot_sobols(
    sobolindices(pce512),
    title,
    bar_width=0.3,
    ymin=0.0,
    ymax=0.4
)
path2save = "/home/perin/Documents/academic/slides/ThermoOptiPlan/5_presentation_INTERNAL/imgs"
PGFPlotsX.save(joinpath(path2save, "Sobols512.pdf"), plot_sobol)

# Tests

@load "/home/perin/Projects/ThermoOptiPlan/Tests/test_12.jld2"
true_values = df.crossing_year

df_test128 = df[:, 1:15]
evaluate!(pce128, df_test128)
pred128 = df_test128.crossing_year
minv128 = min(minimum(pred128), minimum(true_values))
maxv128 = max(maximum(pred128), maximum(true_values))
mse128 = mean((true_values .- pred128) .^ 2)
mse128 = trunc(mse128 * 100) / 100

p128 = @pgf Axis(
    {
        title = "PCE 128 Sample - Test Performance",
        xlabel = "pce - 128",
        ylabel = "OGS",
        grid = "major",
        legend_pos = "south east"
    },
    Plot(
        {
            only_marks,
            color = "blue",
        },
        Coordinates(zip(pred128, true_values))
    ),
    LegendEntry("MSE =" * string(mse128)),
    Plot(
        {
            dashed,
            color = "red"
        },
        Coordinates([(minv128, minv128), (maxv128, maxv128)])
    )
)
path2save = "/home/perin/Documents/academic/slides/ThermoOptiPlan/5_presentation_INTERNAL/imgs"
PGFPlotsX.save(joinpath(path2save, "test_performance128.pdf"), p128)


df_test256 = df[:, 1:15]
evaluate!(pce256, df_test256)
pred256 = df_test256.crossing_year
minv256 = min(minimum(pred256), minimum(true_values))
maxv256 = max(maximum(pred256), maximum(true_values))
mse256 = mean((true_values .- pred256) .^ 2)
mse256 = trunc(mse256 * 100) / 100
p256 = @pgf Axis(
    {
        title = "PCE 256 Sample - Test Performance",
        xlabel = "pce - 256",
        ylabel = "OGS",
        grid = "major",
        legend_pos = "south east"
    },
    Plot(
        {
            only_marks,
            color = "blue",
        },
        Coordinates(zip(pred256, true_values))
    ),
    LegendEntry("MSE =" * string(mse256)),
    Plot(
        {
            dashed,
            color = "red"
        },
        Coordinates([(minv256, minv256), (maxv256, maxv256)])
    )
)
path2save = "/home/perin/Documents/academic/slides/ThermoOptiPlan/5_presentation_INTERNAL/imgs"
PGFPlotsX.save(joinpath(path2save, "test_performance256.pdf"), p256)

df_test512 = df[:, 1:15]
evaluate!(pce512, df_test512)
pred512 = df_test512.crossing_year
minv512 = min(minimum(pred512), minimum(true_values))
maxv512 = max(maximum(pred512), maximum(true_values))
mse512 = mean((true_values .- pred512) .^ 2)
mse512 = trunc(mse512 * 100) / 100

p512 = @pgf Axis(
    {
        title = "PCE 512 Sample - Test Performance",
        xlabel = "pce - 512",
        ylabel = "OGS",
        grid = "major",
        legend_pos = "south east"
    },
    Plot(
        {
            only_marks,
            color = "blue",
        },
        Coordinates(zip(pred512, true_values))
    ),
    LegendEntry("MSE =" * string(mse512)),
    Plot(
        {
            dashed,
            color = "red"
        },
        Coordinates([(minv512, minv512), (maxv512, maxv512)])
    )
)
path2save = "/home/perin/Documents/academic/slides/ThermoOptiPlan/5_presentation_INTERNAL/imgs"
PGFPlotsX.save(joinpath(path2save, "test_performance512.pdf"), p512)
