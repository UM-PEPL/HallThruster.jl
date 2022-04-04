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
using PartialFunctions
import Term.progress: track, ProgressBar, update, start, stop
import Term.progress: AbstractColumn, DescriptionColumn, BarColumn
import Term.style: apply_style  # to style strings
import Term.measure: Measure  # to define the column's size
using MuladdMacro
using SpecialFunctions
using QuadGK

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")
const LANDMARK_RATES_FILE = joinpath(LANDMARK_FOLDER, "landmark_rates.csv")

include("numerics/finite_differences.jl")
include("transition_functions.jl")
include("physics/physicalconstants.jl")
include("physics/gas.jl")
include("physics/conservationlaws.jl")
include("physics/fluid.jl")
include("physics/thermodynamics.jl")
include("numerics/limiters.jl")
include("numerics/flux.jl")
include("reactions/reactions.jl")
include("reactions/ionization.jl")
include("reactions/excitation.jl")
include("geometry.jl")
include("physics/boundaryconditions.jl")
include("physics/potential.jl")
include("postprocess.jl")
include("physics/heavy_species.jl")
include("physics/electronenergy.jl")
include("update_values.jl")
include("physics/sourceterms.jl")
include("physics/electrontransport.jl")
include("physics/collisions.jl")
include("progressbar.jl")
include("physics/walls.jl")
include("configuration.jl")
include("utility_functions.jl")
include("simulation.jl")

end # module
