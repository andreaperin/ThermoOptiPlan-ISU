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

# Run f(); while the (network) drive is unreachable, wait 1 min and retry, up to `minutes`
function wait_for_drive(f; minutes::Int=30)
    for attempt in 1:minutes
        try
            return f()
        catch err
            (err isa Base.IOError || err isa SystemError) || rethrow()
            attempt == minutes && error("Drive unreachable for $minutes min. Rerun the same call to resume.")
            @warn "Drive unreachable, retrying in 1 min ($attempt/$minutes)"
            sleep(60)
        end
    end
end

# One OGS sample folder -> (inputs..., extraction_temperatures), or nothing if it has no results
function read_ogs_sample(path::String)
    any(endswith("000.vtu"), readdir(path)) || return nothing

    # every random input is written as <value(s)>number</value(s)> on the line of its placeholder
    template = readlines(joinpath(SOURCE_DIR, sourcefile))
    rendered = readlines(joinpath(path, sourcefile))
    values = Dict{Symbol,Float64}()
    for (t, r) in zip(template, rendered)
        m = match(r"\{\{\{ :(\w+) \}\}\}", t)
        m === nothing && continue
        values[Symbol(m.captures[1])] = parse(Float64, match(r"<values?>\s*([^<\s]+)\s*</values?>", r).captures[1])
    end

    temperatures = extraction_temperatures_over_time_fast(path, x_extractor, y_extractor, Δz_extractor)
    return (; values..., extraction_temperatures=temperatures)
end

"""
    dataframe_from_ogs(run_dir; checkpoint)

Same as `dataframe_from_pce`, but read directly from the OGS output folders
`run_dir/sample-*`: inputs from the rendered .prj, temperatures from the .vtu files
(only the T array is read, so it is fast also on network drives).

Resumable: every finished sample is saved to `checkpoint` (default: a file in DATASETS_DIR
named after `run_dir`); calling it again skips those. If the drive is unreachable it waits
and retries. Samples without results or with unreadable .vtu files are skipped.
"""
function dataframe_from_ogs(run_dir::String;
        checkpoint::String=joinpath(DATASETS_DIR, "checkpoint_" * basename(normpath(run_dir)) * ".jld2"))
    done = isfile(checkpoint) ? load(checkpoint, "done") : Dict{String,Any}()
    isempty(done) || println("Resuming: ", length(done), " samples already in ", checkpoint)
    mkpath(dirname(checkpoint))

    samples = filter(startswith("sample-"), wait_for_drive(() -> readdir(run_dir)))
    for sample in samples
        haskey(done, sample) && continue
        row = try
            wait_for_drive(() -> read_ogs_sample(joinpath(run_dir, sample)))
        catch err
            err isa ErrorException && startswith(err.msg, "Drive unreachable") && rethrow()
            @warn "Skipping $sample: a .vtu file could not be read, even after 3 attempts"
            continue
        end
        row === nothing && continue
        done[sample] = row
        jldsave(checkpoint; done)
    end

    rows = [done[s] for s in sort(collect(keys(done)))]
    return inputs_and_output(DataFrame(rows))
end
