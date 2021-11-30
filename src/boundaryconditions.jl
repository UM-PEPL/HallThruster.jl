abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition
    state::Vector{Float64}
end

struct Neumann <: BoundaryCondition end

function apply_bc!(U, bc::Dirichlet, left_or_right::Symbol)
    if left_or_right == :left
        @. @views U[:, begin] = bc.state
    elseif left_or_right == :right
        @. @views U[:, end] = bc.state
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc!(U, ::Neumann, left_or_right::Symbol)
    if left_or_right == :left
        @. @views U[:, begin] = U[:, begin + 1]
    elseif left_or_right == :right
        @. @views U[:, end] = U[:, end - 1]
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end