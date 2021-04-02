# test script for developing GMRES algorithm
# before running, make sure to navigate to the top directory in `MyPanel` and run
# `] activate .`

import MyPanel
pnl = MyPanel
using Statistics
using LinearAlgebra
using PyPlot

# GeometricTools module
gt = MyPanel.GeometricTools

module_path = normpath(joinpath(@__DIR__,"..")) # Path to this module
data_path = pnl.def_data_path         # Data path
airfoil_path = joinpath(data_path, "airfoils"); # Airfoil data path

save_path = joinpath(module_path, "../temps")
file_name = "sphere02"

# Parameters
nu = 1.443e-5                    # (m^2/s) kinematic viscosity
Re = 8800                        # Reynolds number V*d/nu
# R = 6.35e-3/2                    # (m) radius of sphere
# magVinf = 20                     # (m/s) freestream velocity   
R = 1
magVinf = Re*nu/(2*R)

P_min = [0.15, 0, 0]             # Lower bounds of (theta, phi, dummy)
P_max = [pi-P_min[1], 2*pi, 0]   # Upper bounds of (theta, phi, dummy)
NDIVS = 1*[15, 30, 0]            # Number of divisions (cells) of (theta, phi, dummy)
loop_dim = 2                     # Coordinate to loop (1==theta)

# Generates parametric (theta, phi) grid 
grid = pnl.gt.Grid(P_min, P_max, NDIVS, loop_dim)

# Transforms the grid into a spherical cartesian space
my_transform(X) = pnl.gt.spherical3D(vcat(R, X[1:2]))
pnl.gt.transform!(grid, my_transform)

# Splits the quad cells into triangular cells
dimsplit = 1
triang_grid = pnl.gt.GridTriangleSurface(grid, dimsplit)

# Creates non lifting body
body = pnl.NonLiftingBody(triang_grid)

# Adds normal vector field

pnl.gt.add_field(body.grid, "normal", "vector",
                    [pnl.gt.get_normal(body.grid, i) for i in 1:body.ncells], "cell")

pnl.gt.add_field(body.grid, "cellindex", "scalar",
                    [i for i in 1:body.ncells], "cell")
pnl.gt.add_field(body.grid, "nodeindex", "scalar",
                    [i for i in 1:body.nnodes], "node")

# Freestream at every control point
Vinf = magVinf*[1.0,0,0]
Vinfs = [Vinf for i in 1:body.ncells]
pnl.solve(body, Vinfs)

# Adds surface velocity field
CPs = [pnl.get_controlpoint(body, i) for i in 1:body.ncells]
Vsurf = [Vinf for i in 1:size(CPs,1)]
for i in 1:body.ncells
    pnodes = gt.get_cellnodes(body.grid, i)
    pnl.PanelSolver.Vconstant_source(pnodes, pnl.get_fieldval(body, "sigma", i), CPs, Vsurf)
end
point_data = [Dict("field_name"=>"V", "field_type"=>"vector", "field_data"=>Vsurf)]
pnl.gt.add_field(body.grid, "V", "vector", Vsurf, "cell")

# Creates a fluid domain grid
fdom = gt.Grid(-3*R*ones(3), 3*R*ones(3), 2*[10,10,10])
targets = [gt.get_node(fdom, i) for i in 1:fdom.nnodes]
V = [Vinf for i in 1:fdom.nnodes]
for i in 1:body.ncells
    pnodes = gt.get_cellnodes(body.grid, i)
    pnl.PanelSolver.Vconstant_source(pnodes, pnl.get_fieldval(body, "sigma", i), targets, V)
end

gt.add_field(fdom, "V", "vector", V, "node")

# Saves vtk and calls paraview
pnl.save(body, file_name; path=save_path)
gt.save(fdom, file_name*"_fdom"; path=save_path)
gt.generateVTK(file_name*"_CPs", CPs; point_data=point_data, path=save_path)
strn = "$(joinpath(save_path, file_name)).vtk;$(file_name)_fdom.vtk;$(file_name)_CPs.vtk"
run(`paraview --data=$strn`)
