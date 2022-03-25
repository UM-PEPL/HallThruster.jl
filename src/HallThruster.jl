module HallThruster

using StaticArrays
using CSV
using DataFrames
using OrdinaryDiffEq
using DiffEqBase
using LoopVectorization
using LinearAlgebra
using FileIO
using DiffEqCallbacks
using SparseArrays
using ForwardDiff
using Term
import Term.progress: track, ProgressBar, update, start, stop
import Term.progress: AbstractColumn, DescriptionColumn, BarColumn
import Term.style: apply_style  # to style strings
import Term.measure: Measure  # to define the column's size

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")

include("finite_differences.jl")
include("transition_functions.jl")
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
include("heavy_species.jl")
include("electronenergy.jl")
include("update_values.jl")
include("sourceterms.jl")
include("electrontransport.jl")
include("collisions.jl")
include("progressbar.jl")
include("configuration.jl")
include("utility_functions.jl")
include("simulation.jl")
include("walls.jl")

end # module
