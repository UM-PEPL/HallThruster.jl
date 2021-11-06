abstract type BoundaryCondition end

struct Dirichlet  <: BoundaryCondition 
    state::Vector{Float64}
end

struct Neumann <: BoundaryCondition end