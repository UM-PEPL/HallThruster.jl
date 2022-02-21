module HallThruster

using StaticArrays
using CSV
using DataFrames
using OrdinaryDiffEq
using DiffEqBase
using LoopVectorization
using LinearAlgebra
using FileIO

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")

include("utility_functions.jl")
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
include("potential.jl")
include("postprocess.jl")
include("simulation.jl")
include("sourceterms.jl")
include("electrontransport.jl")
include("electronenergy.jl")
end # module
