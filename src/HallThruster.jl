module HallThruster

using StaticArrays
using CSV
using DataFrames
using OrdinaryDiffEq
using DiffEqBase
using LoopVectorization

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:end-1]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")

include("physicalconstants.jl")
include("gas.jl")
include("conservationlaws.jl")
include("fluid.jl")
include("thermodynamics.jl")
include("limiters.jl")
include("flux.jl")
include("ionization.jl")
include("geometry.jl")
include("boundaryconditions.jl")
include("simulation.jl")
include("sourceterms.jl")

#dere
# this change was made on main

end # module
