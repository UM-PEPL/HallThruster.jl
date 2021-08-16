module HallThruster

using StaticArrays

include("physicalconstants.jl")
include("gas.jl")
include("conservationlaws.jl")
include("fluid.jl")
include("thermodynamics.jl")
include("limiters.jl")
include("flux.jl")
include("ionization.jl")
include("solve.jl")

end # module
