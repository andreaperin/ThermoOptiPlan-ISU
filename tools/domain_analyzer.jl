using Pkg

domain_vtu_file = "models/multilayer_irz/Multi_BW_line_IRZ_domain_ini.vtu"
pump2_vtu_file = "models/multilayer_irz/Multi_BW_line_IRZ_physical_group_Pump_line_2.vtu"
pump3_vtu_file = "models/multilayer_irz/Multi_BW_line_IRZ_physical_group_Pump_line_3.vtu"
pump4_vtu_file = "models/multilayer_irz/Multi_BW_line_IRZ_physical_group_Pump_line_4.vtu"

if Sys.islinux()
    osrelease = "/etc/os-release"
    data = read(osrelease, String)
    if occursin("Solus", data)
        ENV["PYTHON"] = normpath(joinpath(@__DIR__, "..", ".venv", "bin", "python"))
    elseif occursin("NixOS", data)
        ENV["PYTHON"] = "/home/lau/python_venv/bin/python"
    end
elseif Sys.iswindows()
    ENV["PYTHON"] = "C:/Users/laure/AppData/Local/Programs/Python/Python313/python.exe"
elseif Sys.isapple()
    ENV["PYTHON"] = "/Users/andreaperin_macos/Documents/Code/5_OpenGeoSys/.venv/bin/python"
else
    # ENV["PYTHON"] = "" # Default: try system Python or let PyCall auto-detect
end

Pkg.build("PyCall")

using PyCall
using Statistics
pv = pyimport("pyvista")

mesh = pv.read(domain_vtu_file)
bounds = mesh.bounds
xmin, xmax, ymin, ymax, zmin, zmax = bounds

println("MESH Bounds")
println("x: ", xmin, " → ", xmax)
println("y: ", ymin, " → ", ymax)
println("z: ", zmin, " → ", zmax)

println("Extraction Coordinates Bounds")
println("x: ", xmin, " → ", xmax)
println("y: ", ymin, " → ", ymax)
println("z: ", zmin, " → ", zmax)

pump2_mesh = pv.read(pump2_vtu_file)
points = Array(pump2_mesh.points)
x = points[:, 1]
y = points[:, 2]
z = points[:, 3]
x_well = mean(x)
y_well = mean(y)
println("PUMP2 location:")
println("x = ", x_well)
println("y = ", y_well)
println("z range: ", minimum(z), " → ", maximum(z))
pump3_mesh = pv.read(pump3_vtu_file)
points = Array(pump3_mesh.points)
x = points[:, 1]
y = points[:, 2]
z = points[:, 3]
x_well = mean(x)
y_well = mean(y)
println("PUMP3 location:")
println("x = ", x_well)
println("y = ", y_well)
println("z range: ", minimum(z), " → ", maximum(z))
pump4_mesh = pv.read(pump4_vtu_file)
points = Array(pump4_mesh.points)
x = points[:, 1]
y = points[:, 2]
z = points[:, 3]
x_well = mean(x)
y_well = mean(y)
println("PUMP4 location:")
println("x = ", x_well)
println("y = ", y_well)
println("z range: ", minimum(z), " → ", maximum(z))

material_ids = Array(mesh.cell_data.get_array("MaterialIDs"))
centers = mesh.cell_centers()
points = Array(centers.points)
x = points[:, 1]
y = points[:, 2]
z = points[:, 3]
function layer_bounds(ids, z, mat_id)
    z_layer = z[ids.==mat_id]
    return minimum(z_layer), maximum(z_layer)
end

println("SANSTONE location:")
for id in [1, 2, 3]
    zmin, zmax = layer_bounds(material_ids, z, id)
    println("Sandstone (ID=$id): ", zmin, " → ", zmax)
end

println("CLAY location:")
for id in [0, 4]
    zmin, zmax = layer_bounds(material_ids, z, id)
    println("Clay (ID=$id): ", zmin, " → ", zmax)
end