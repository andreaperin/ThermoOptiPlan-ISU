using Pkg

if Sys.islinux()
    osrelease = "/etc/os-release"
    data = read(osrelease, String)
    if occursin("Solus", data) || occursin("Ubuntu", data)
        ENV["PYTHON"] = normpath(joinpath(@__DIR__, "..", ".venv", "bin", "python"))
    elseif occursin("NixOS", data)
        ENV["PYTHON"] = "/home/lau/python_venv/bin/python"
    end
end

# Rebuild PyCall only if it was built for another python or its libpython is gone
function pycall_needs_build(python::String)
    depfile = joinpath(dirname(Base.find_package("PyCall")), "..", "deps", "deps.jl")
    isfile(depfile) || return true
    deps = read(depfile, String)
    built_python = match(r"const python = \"(.*)\"", deps)
    libpython = match(r"const libpython = \"(.*)\"", deps)
    (built_python === nothing || libpython === nothing) && return true
    return built_python.captures[1] != python || !isfile(libpython.captures[1])
end

if haskey(ENV, "PYTHON") && pycall_needs_build(ENV["PYTHON"])
    Pkg.build("PyCall")
end

using PyCall
using Statistics
pv = pyimport("pyvista")

function _average_extraction_temperature(outputfile::String, mask)
    mesh = pv.read(outputfile)
    T_data = mesh.point_data.get_array("T")

    Ts_extraction = T_data[mask]
    time_regex = r"t_(\d+(?:\.\d+)?)\.vtu"

    m = match(time_regex, outputfile)
    number_str = m.captures[1]
    number = parse(Float64, number_str)
    Δyear = number / 365 / 24 / 60 / 60
    return [mean(Ts_extraction), Δyear]
end

function extraction_temperatures_over_time(output_path::String, x::Float64, y::Float64, Δz::Vector{Tuple{Float64,Float64}}; tol_xy::Real=0.5)
    vtu_files = filter(f -> endswith(f, "000.vtu"), readdir(output_path))

    first_mesh = pv.read(joinpath(output_path, vtu_files[1]))
    points_coords = first_mesh.points

    # build masks for each layer
    masks = [
        map(row -> abs(row[1] - x) < tol_xy &&
                       abs(row[2] - y) < tol_xy &&
                       zmin ≤ row[3] ≤ zmax,
            eachrow(points_coords))
        for (zmin, zmax) in Δz
    ]
    results = map(vtu_file -> begin
            mesh = pv.read(joinpath(output_path, vtu_file))
            T_data = mesh.point_data.get_array("T")

            # compute mean T for each layer
            T_means = [mean(T_data[mask]) for mask in masks]

            # extract time
            time_regex = r"t_(\d+(?:\.\d+)?)\.vtu"
            m = match(time_regex, vtu_file)
            number = parse(Float64, m.captures[1])
            Δyear = number / 365 / 24 / 60 / 60

            return vcat(T_means, Δyear)
        end, vtu_files)

    return sort(results, by=x -> x[end])
end


# x = 2_250.0
# y = 0.0
# Δz_bottom = (-1374.8, -1364.6)
# Δz_middle = (-1338.0, -1309.1)
# Δz_top = (-1265.4, -1240.2)
# Δz = [
#     Δz_bottom,
#     Δz_middle,
#     Δz_top
# ]

# output_path = "/home/perin/Documents/projects/work/code/thermoptiplan_new/output/Model_ML_IRZ/2026-04-09-17-57-13/sample-1"

# extraction_temperatures = Vector{Vector{Any}}()
# push!(extraction_temperatures, extraction_temperatures_over_time(output_path, x, y, Δz))

# flows = Vector{Vector{Real}}([[0.3, 0.3, 0.4]])


# function final_temperature(extraction_temperatures::Vector{Vector{<:Real}}, flows::Vector{<:Real})
#     return [[sum(v[1:3] .* flows), v[4]] for v in extraction_temperatures]
# end