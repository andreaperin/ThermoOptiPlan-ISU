using JLD2
using DataFrames

include(joinpath(@__DIR__, "..", "src", "model_3layers.jl"))

const input_names = [string(rv.name) for rv in inputs]

# Post-process the extraction temperatures with the current models and keep inputs + output.
# Incomplete OGS runs (crossing_year = NaN) are dropped.
function inputs_and_output(df::DataFrame)
    for m in postprocessing_models
        evaluate!(m, df)
    end
    df = df[.!isnan.(df.crossing_year), :]
    return df[:, [input_names; "crossing_year"]]
end

"""
    dataframe_from_pce(files)

DataFrame with the 15 inputs and `crossing_year` of every complete OGS run stored in the
PCE result files `files` (each holding `res = [pce, samples, mse]`).
"""
function dataframe_from_pce(files::Vector{String})
    df = reduce(vcat, [load(f, "res")[2][:, [input_names; "extraction_temperatures"]] for f in files])
    return inputs_and_output(df)
end

"""
    dataframe_from_ogs(run_dir)

Same as `dataframe_from_pce`, but read directly from the OGS output folders
`run_dir/sample-*`: inputs from the rendered .prj, temperatures from the .vtu files.
Samples without results are skipped.
"""
function dataframe_from_ogs(run_dir::String)
    rows = []
    for sample in filter(startswith("sample-"), readdir(run_dir))
        path = joinpath(run_dir, sample)
        any(endswith("000.vtu"), readdir(path)) || continue

        # every random input is written as <value(s)>number</value(s)> on the line of its placeholder
        template = readlines(joinpath(SOURCE_DIR, sourcefile))
        rendered = readlines(joinpath(path, sourcefile))
        values = Dict{Symbol,Float64}()
        for (t, r) in zip(template, rendered)
            m = match(r"\{\{\{ :(\w+) \}\}\}", t)
            m === nothing && continue
            values[Symbol(m.captures[1])] = parse(Float64, match(r"<values?>\s*([^<\s]+)\s*</values?>", r).captures[1])
        end

        temperatures = extraction_temperatures_over_time(path, x_extractor, y_extractor, Δz_extractor)
        push!(rows, (; values..., extraction_temperatures=temperatures))
    end
    return inputs_and_output(DataFrame(rows))
end
