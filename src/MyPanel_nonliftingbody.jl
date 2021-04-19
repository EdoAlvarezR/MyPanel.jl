#=##############################################################################
# DESCRIPTION
    Non-lifting paneled body types definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : AGPL-3.0
=###############################################################################


################################################################################
# NON-LIFTING BODY TYPE
################################################################################

# abstract type CompactG end

# function overload_getindex(grid::gt.GridTTriangleSurface)
#     Base.getindex(_G::CompactG) = 4

# end

"""
  `NonLiftingBody(grid::gt.GridTriangleSurface)`

Non-lifting paneled body that is solved using a constant source distribution.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
  * `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid

"""
struct NonLiftingBody <: AbstractBody

  # User inputs
  grid::gt.GridTriangleSurface              # Paneled geometry

  # Properties
  nnodes::Int64                             # Number of nodes
  ncells::Int64                             # Number of cells
  fields::Array{String, 1}                  # Available fields (solutions)
  Oaxis::Array{T1,2} where {T1<:RType}      # Coordinate system of original grid
  O::Array{T2,1} where {T2<:RType}          # Position of CS of original grid

  # Internal variables
  _G::Array{T3,2} where {T3<:RType}         # Geometric solution matrix
  _nodes::Array{T4,2} where {T4<:RType}
  _panels::Array{Array{Int64,1},1}
  _CPs::Array{Array{T5,1},1} where {T5<:RType}
  _normals::Array{Array{T6,1},1} where {T6<:RType}

  NonLiftingBody( grid,
                  nnodes=grid.nnodes, ncells=grid.ncells,
                    fields=Array{String,1}(),
                    Oaxis=Array(1.0I, 3, 3), O=zeros(3),
                  _G=_calc_Gsource(grid),
                  _nodes=grid.orggrid.nodes,                                   # Nodes
                  _panels=[gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  _CPs=[_get_controlpoint(grid, i) for i in 1:grid.ncells],    # CPs
                  _normals=[gt.get_normal(grid, i) for i in 1:grid.ncells]     # Normals
         ) = new( grid,
                  nnodes, ncells,
                    fields,
                    Oaxis, O,
                  _G, _nodes, _panels, _CPs, _normals
         )
end

"Function gets the influence on the ith control point of the jth panel of unit strength"
function Base.getindex(self::NonLiftingBody, i::Int, j::Int)
    N = self.ncells
    # Gij = 0.0
    # for j in 1:N # Iterates over columns (panels)
    Gij = PanelSolver.Vconstant_source_compact(
            [self._nodes[:,ind] for ind in self._panels[j]], # Nodes in j-th panel
            1.0,                                 # Unitary strength,
            self._CPs[i];                        # Target control point
            dot_with=self._normals[i]                     # Normal of every CP
        )
    # end
    return Gij
end

Base.size(A::NonLiftingBody, index::Int64) = A.ncells

function IS.LinearAlgebra.mul!(C::AbstractArray, A::NonLiftingBody, B::AbstractArray, alpha::Number, beta::Number)
    C .*= beta
    D = zeros(size(C))
    @inbounds for j = 1:length(B)
        @simd for i = 1:length(C)
            D[i] += A[i,j] * B[j]
        end
    end
    C .+= alpha .* D
end

function IS.LinearAlgebra.mul!(C::AbstractArray, A::NonLiftingBody, B::AbstractArray, alpha::Bool, beta::Bool)
    if !alpha && !beta
        return 0 .* C
    elseif !alpha && beta
        return C
    elseif alpha && !beta
        C .*= 0
        @inbounds for j = 1:length(B)
            @simd for i = 1:length(C)
                C[i] += A[i,j] * B[j]
            end
        end
        return C
    else
        return IS.LinearAlgebra.mul!(C, A, B, 1, 1)
    end
end


"Overload `Adivtype` to work with `::NonLiftingBody` type."
IS.Adivtype(A::NonLiftingBody,b) = typeof(one(eltype(b))/one(eltype(A[1,1])))

function solve(self::NonLiftingBody, Vinfs::Array{Array{T,1},1};
    algorithm::SolverAlgorithm = gmres,
    verbose=true) where {T<:RType}
  if size(Vinfs,1) != self.ncells
    error("Invalid Vinfs; expected size $(self.ncells), got $(size(Vinfs,1))")
  end

  lambda = [-dot(Vinfs[i], get_normal(self, i)) for i in 1:self.ncells]

  if Int(algorithm) == 1
    sigma = self._G\lambda
  elseif Int(algorithm) == 2
    sigma = IS.gmres(self._G, lambda; verbose=verbose)
  elseif Int(algorithm) == 3
    sigma = IS.gmres(self, lambda; verbose=verbose)
  elseif Int(algorithm) == 4
    sigma = IS.cg(self._G, lambda; verbose=verbose)
  else
    throw("algorithm not implemented.")
  end

  add_field(self, "Vinf", Vinfs)
  add_field(self, "sigma", sigma)
  _solvedflag(self, true)
end

##### INTERNAL FUNCTIONS  ######################################################
function _calc_Gsource(grid::gt.GridTriangleSurface)
  return PanelSolver.G_constant_source(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells],     # Normals
                )
end

"""
Returns the velocity induced by the body on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
function _Vind(self::NonLiftingBody, targets::Array{Array{T1,1},1},
                          out::Array{Array{T2,1},1}) where{T1<:RType, T2<:RType}
  # Iterates over panels
  for i in 1:self.ncells
    # Velocity of i-th  panel on every target
    PanelSolver.Vconstant_source(
                    gt.get_cellnodes(self.grid, i),    # Nodes in i-th panel
                    get_fieldval(self, "sigma", i; _check=false),  # Strength
                    targets,                           # Targets
                    out;                               # Outputs
                  )
  end
end
##### END OF NON-LIFTING BODY ##################################################


################################################################################
# NON-LIFTING DOUBLET BODY TYPE
################################################################################
"""
  `NonLiftingBodyDoublet(grid::gt.GridTriangleSurface)`

Non-lifting paneled body that is solved using a constant doublet distribution.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
  * `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid

"""
struct NonLiftingBodyDoublet <: AbstractBody

  # User inputs
  grid::gt.GridTriangleSurface              # Paneled geometry

  # Properties
  nnodes::Int64                             # Number of nodes
  ncells::Int64                             # Number of cells
  fields::Array{String, 1}                  # Available fields (solutions)
  Oaxis::Array{T1,2} where {T1<:RType}      # Coordinate system of original grid
  O::Array{T2,1} where {T2<:RType}          # Position of CS of original grid

  # Internal variables
  _G::Array{T3,2} where {T3<:RType}         # Geometric solution matrix

  NonLiftingBodyDoublet( grid,
                  nnodes=grid.nnodes, ncells=grid.ncells,
                    fields=Array{String,1}(),
                    Oaxis=eye(3), O=zeros(3),
                  _G=_calc_Gdoublet(grid)
         ) = new( grid,
                  nnodes, ncells,
                    fields,
                    Oaxis, O,
                  _G
         )
end


function solve(self::NonLiftingBodyDoublet, Vinfs::Array{Array{T,1},1};
                algorithm::SolverAlgorithm = gmres,
                verbose=true) where {T<:RType}
  if size(Vinfs,1) != self.ncells
    error("Invalid Vinfs; expected size $(self.ncells), got $(size(Vinfs,1))")
  end

  lambda = [-dot(Vinfs[i], get_normal(self, i)) for i in 1:self.ncells]

  if Int(algorithm) == 1
    mu = self._G\lambda
  elseif Int(algorithm) == 2
    mu = IS.gmres(self._G, lambda; verbose=verbose)
  elseif Int(algorithm) == 3
    mu = IS.gmres(self, lambda; verbose=verbose)
  elseif Int(algorithm) == 4
    mu = IS.cg(self._G, lambda; verbose=verbose)
  else
    throw("algorithm not implemented.")
  end

  add_field(self, "Vinf", Vinfs)
  add_field(self, "mu", mu)
  _solvedflag(self, true)
end




##### INTERNAL FUNCTIONS  ######################################################
function _calc_Gdoublet(grid::gt.GridTriangleSurface)
  return PanelSolver.G_constant_doublet(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells],     # Normals
                )
end

"""
Returns the velocity induced by the body on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
function _Vind(self::NonLiftingBodyDoublet, targets::Array{Array{T1,1},1},
                          out::Array{Array{T2,1},1}) where{T1<:RType, T2<:RType}
  # Iterates over panels
  for i in 1:self.ncells
    # Velocity of i-th  panel on every target
    PanelSolver.Vconstant_doublet(
                    gt.get_cellnodes(self.grid, i),    # Nodes in i-th panel
                    get_fieldval(self, "mu", i; _check=false),  # Strength
                    targets,                           # Targets
                    out;                               # Outputs
                  )
  end
end
##### END OF NON-LIFTING DOUBLET BODY ##########################################




################################################################################
# NON-LIFTING VORTEX-RING BODY TYPE
################################################################################
"""
  `NonLiftingBodyVRing(grid::gt.GridTriangleSurface)`

Non-lifting paneled body that is solved using vortex ring panels.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
  * `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid

"""
struct NonLiftingBodyVRing <: AbstractBody

  # User inputs
  grid::gt.GridTriangleSurface              # Paneled geometry

  # Properties
  nnodes::Int64                             # Number of nodes
  ncells::Int64                             # Number of cells
  fields::Array{String, 1}                  # Available fields (solutions)
  Oaxis::Array{T1,2} where {T1<:RType}      # Coordinate system of original grid
  O::Array{T2,1} where {T2<:RType}          # Position of CS of original grid

  # Internal variables
  _G::Array{T3,2} where {T3<:RType}         # Geometric solution matrix

  NonLiftingBodyVRing( grid,
                  nnodes=grid.nnodes, ncells=grid.ncells,
                    fields=Array{String,1}(),
                    Oaxis=eye(3), O=zeros(3),
                  _G=_calc_Gvring(grid)
         ) = new( grid,
                  nnodes, ncells,
                    fields,
                    Oaxis, O,
                  _G
         )
end


function solve(self::NonLiftingBodyVRing, Vinfs::Array{Array{T,1},1};
                algorithm::SolverAlgorithm = gmres,
                verbose=true) where {T<:RType}
  if size(Vinfs,1) != self.ncells
    error("Invalid Vinfs; expected size $(self.ncells), got $(size(Vinfs,1))")
  end

  lambda = [-dot(Vinfs[i], get_normal(self, i)) for i in 1:self.ncells]

  if Int(algorithm) == 1
    Gamma = self._G\lambda
  elseif Int(algorithm) == 2
    Gamma = IS.gmres(self._G, lambda; verbose=verbose)
  elseif Int(algorithm) == 3
    Gamma = IS.gmres(self, lambda; verbose=verbose)
  elseif Int(algorithm) == 4
    Gamma = IS.cg(self._G, lambda; verbose=verbose)
  else
    throw("algorithm not implemented.")
  end

  add_field(self, "Vinf", Vinfs)
  add_field(self, "Gamma", Gamma)
  _solvedflag(self, true)
end



##### INTERNAL FUNCTIONS  ######################################################
function _calc_Gvring(grid::gt.GridTriangleSurface)
  return PanelSolver.G_vortexring(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells],     # Normals
                )
end

"""
Returns the velocity induced by the body on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
function _Vind(self::NonLiftingBodyVRing, targets::Array{Array{T1,1},1},
                          out::Array{Array{T2,1},1}) where{T1<:RType, T2<:RType}
  # Iterates over panels
  for i in 1:self.ncells
    # Velocity of i-th  panel on every target
    PanelSolver.Vvortexring(
                    gt.get_cellnodes(self.grid, i),    # Nodes in i-th panel
                    get_fieldval(self, "Gamma", i; _check=false),  # Strength
                    targets,                           # Targets
                    out;                               # Outputs
                  )
  end
end
##### END OF NON-LIFTING VORTEX-RING BODY ######################################





################################################################################
# COMMON FUNCTIONS
################################################################################
"""
  `generate_loft_nonliftbody(args...; optargs...)`
Generates a lofted non-lifting body. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft_nonliftbody(args...; vortexring=false, dimsplit::Int64=2,
                                                                    optargs...)
  # Lofts the surface geometry
  grid = gt.generate_loft(args...; optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  if vortexring
    return NonLiftingBodyVRing(triang_grid)
  else
    return NonLiftingBody(triang_grid)
  end
end

"""
  `generate_revolution_nonliftbody(args...; optargs...)`
Generates a non-lifting body of a body of revolution. See documentation of
`GeometricTools.surface_revolution` for a description of the arguments of this
function.
"""
function generate_revolution_nonliftbody(args...; vortexring=false,
                                                  dimsplit::Int64=2,
                                                  loop_dim::Int64=2, optargs...)
  # Revolves the geometry
  grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  if vortexring
    return NonLiftingBodyVRing(triang_grid)
  else
    return NonLiftingBody(triang_grid)
  end
end
##### END OF COMMON FUNTIONS ###################################################
