# test script for developing GMRES algorithm
# before running, make sure to navigate to the top directory in `MyPanel` and run
# `] activate .`

import MyPanel
pnl = MyPanel
using Statistics
using LinearAlgebra
using PyPlot
using BenchmarkTools
using Infiltrator

# GeometricTools module
gt = MyPanel.GeometricTools

module_path = normpath(joinpath(@__DIR__,"..")) # Path to this module
data_path = pnl.def_data_path         # Data path
airfoil_path = joinpath(data_path, "airfoils"); # Airfoil data path

save_path = joinpath(module_path, "../temps")
file_name = "sphere02"

"""
    `benchmarking_gmres()`

Runs the sphere example and benchmarks the solve method from
    MyPanel using the specified algorithm.

ARGUMENTS:
`benchmarking::Bool` - either runs benchmarking if set to true,
    or if set to false runs the example and pulls up the
    visualization in Paraview.
<<<<<<< HEAD
`algorithm::String` - for now, either pnl.native or pnl.gmres
    will work.
=======
`algorithm::String`:
        - "all" runs all algorithms (currently native solve and gmres)
        - "native" runs native linear solve
        - "gmres" runs with GMRES
>>>>>>> 7f80c41c97ac5ff8c79f91f14c5974db7a6d8db8

Note: The benchmarking macro outputs a lot of useful information
    including min/max/median time and estimated memory used, but only if the
    macro is the last thing to be called. Maybe there's a way around
    this? In the meantime, it means we can't benchmark both algorithms
    right after the other with one command. It requires two commands,
    whether in the REPL or in another function. I also am not able to
    just hardcode the function to run both; only the latter will show
    up in the REPL.
"""
<<<<<<< HEAD
function benchmarking_sphere(;
                            benchmarking=true,
                            algorithm=pnl.gmres_bulky)
=======
function benchmarking_sphere(;
                            benchmarking=true,
                            algorithm="all",
                            panel_scale=3,
                            minimum=false)
>>>>>>> 7f80c41c97ac5ff8c79f91f14c5974db7a6d8db8

    # Parameters
    nu = 1.443e-5                    # (m^2/s) kinematic viscosity
    Re = 8800                        # Reynolds number V*d/nu
    # R = 6.35e-3/2                    # (m) radius of sphere
    # magVinf = 20                     # (m/s) freestream velocity
    R = 1
    magVinf = Re*nu/(2*R)

    P_min = [0.15, 0, 0]             # Lower bounds of (theta, phi, dummy)
    P_max = [pi-P_min[1], 2*pi, 0]   # Upper bounds of (theta, phi, dummy)
    NDIVS = panel_scale*[5, 10, 0]            # Number of divisions (cells) of (theta, phi, dummy)
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
    body = algorithm == pnl.gmres_agile ? pnl.NonLiftingBody(triang_grid, false) : pnl.NonLiftingBody(triang_grid, true)

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

    times = []
    memory = []
    if benchmarking # run benchmarking and ignore vtks
        if minimum
            if lowercase(algorithm)=="all"

                # println("Native linear solve:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.native, verbose=false)
                push!(times, minimum(t).time)
                push!(memory, minimum(t).memory)

                # println("GMRES solver:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.gmres, verbose=false)
                push!(times, minimum(t).time)
                push!(memory, minimum(t).memory)

            elseif lowercase(algorithm)=="native"

                # println("Native linear solve:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.native, verbose=false)
                push!(times, minimum(t).time)
                push!(memory, minimum(t).memory)

            elseif lowercase(algorithm)=="gmres"

                # println("GMRES solver:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.gmres, verbose=false)
                push!(times, minimum(t).time)
                push!(memory, minimum(t).memory)
            end
        else #median time instead
            if lowercase(algorithm)=="all"

                # println("Native linear solve:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.native, verbose=false)
                push!(times, median(t).time)
                push!(memory, median(t).memory)

                # println("GMRES solver:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.gmres, verbose=false)
                push!(times, median(t).time)
                push!(memory, median(t).memory)

            elseif lowercase(algorithm)=="native"

                # println("Native linear solve:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.native, verbose=false)
                push!(times, median(t).time)
                push!(memory, median(t).memory)

            elseif lowercase(algorithm)=="gmres"

                # println("GMRES solver:")
                t = @benchmark pnl.solve($body, $Vinfs; algorithm=pnl.gmres, verbose=false)
                push!(times, median(t).time)
                push!(memory, median(t).memory)
            end
        end

        return times, memory

    else # solve and run paraview
        pnl.solve(body, Vinfs; algorithm=algorithm, verbose=false)

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
    end
end

function panel_sweep(;minimum=false)

    scales = 1:7
    numpanels = 5*10*scales .^ 2 # [50; 200; 450; 800; 1250; 1800; 2450; 3200; 4050; 5000; etc.]
    numalgorithms = 2 # change to 3 if conjugate gradient is implemented

    times = Array{Float64}(undef, length(scales), numalgorithms)
    memory = Array{Float64}(undef, length(scales), numalgorithms)
    for scale in scales
        println("$(numpanels[scale]) panels:")
        times[scale, :], memory[scale,:] = benchmarking_sphere(panel_scale=scale, minimum=minimum)
    end

    # Plot results - time and memory
    figure()
    plot(numpanels, times[:,1]*1e-9, label="Native solver")
    plot(numpanels, times[:,2]*1e-9, label="GMRES")
    # plot(numpanels, times[:,3], label="Conjugate gradient")
    xlabel("Number of panels")
    if minimum
        ylabel("Minimum time in benchmarking [ms]")
    else
        ylabel("Median time in benchmarking [ms]")
    end
    legend()

    figure()
    plot(numpanels, memory[:,1] * 1e-6, label="Native solver")
    plot(numpanels, memory[:,2] * 1e-6, label="GMRES")
    # plot(numpanels, memory[:,3], label="Conjugate gradient")
    xlabel("Number of panels")
    ylabel("Memory in benchmarking [Mb]")
    legend()

    figure()
    plot(numpanels, times[:,1]*1e-9, label="Native solver")
    plot(numpanels, times[:,2]*1e-9, label="GMRES")
    # plot(numpanels, times[:,3], label="Conjugate gradient")
    yscale("log")
    xlabel("Number of panels")
    if minimum
        ylabel("Minimum time in benchmarking [ms]")
    else
        ylabel("Median time in benchmarking [ms]")
    end
    legend()

    figure()
    plot(numpanels, memory[:,1] * 1e-6, label="Native solver")
    plot(numpanels, memory[:,2] * 1e-6, label="GMRES")
    # plot(numpanels, memory[:,3], label="Conjugate gradient")
    yscale("log")
    xlabel("Number of panels")
    ylabel("Memory in benchmarking [Mb]")
    legend()
end
