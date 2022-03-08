abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition
    state::Vector{Float64}
end

struct Neumann <: BoundaryCondition end

struct Dirichlet_ionbohm <: BoundaryCondition
    state::Vector{Float64}
end

struct Neumann_ionbohm <: BoundaryCondition end

struct Dirichlet_energy <: BoundaryCondition
    state::Float64
end

struct Dirichlet_energy_upd_ne <: BoundaryCondition
    int_energy::Float64
end

struct Neumann_energy <: BoundaryCondition
end

function left_boundary_state!(bc_state, U, params)
    (;propellant, Te_L, index, mdot_a, A_ch, fluids) = params
    mi = propellant.m

    un = fluids[1].conservation_laws.u
    bc_state[index.ρn] = mdot_a / A_ch / un

    ne = 0.0
    bohm_velocity = sqrt(e * 2/3 * Te_L / mi)
    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], 1]
        boundary_flux = U[index.ρiui[Z]]
        boundary_velocity = min(-sqrt(Z) * bohm_velocity, boundary_flux / boundary_density)
        boundary_density = boundary_flux / boundary_velocity

        bc_state[index.ρn[Z]] -= boundary_flux / un
        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux

        ne += Z * boundary_density / mi
    end

    bc_state[index.nϵ] = ne * Te_L
end

function right_boundary_state!(bc_state, U, params)
    (;propellant, Te_R, index) = params
    mi = propellant.m

    bc_state[index.ρn] = U[index.ρn, end]

    ne = 0.0
    bohm_velocity = sqrt(e * 2/3 * Te_R / mi)
    for Z in 1:params.config.ncharge
        boundary_density = U[index.ρi[Z], 1]
        boundary_flux = U[index.ρiui[Z]]
        boundary_velocity = max(sqrt(Z) * bohm_velocity, boundary_flux / boundary_density)
        boundary_density = boundary_flux / boundary_velocity

        bc_state[index.ρi[Z]] = boundary_density
        bc_state[index.ρiui[Z]] = boundary_flux

        ne += Z * boundary_density / mi
    end

    bc_state[index.nϵ] = ne * Te_R
end

function apply_bc!(U, bc::Dirichlet, left_or_right::Symbol, ϵ0::Float64, mᵢ::Float64)
    if left_or_right == :left
        @. @views U[:, begin] = bc.state
    elseif left_or_right == :right
        @. @views U[:, end] = bc.state
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc!(U, ::Neumann, left_or_right::Symbol, ϵ0::Float64, mᵢ::Float64)
    if left_or_right == :left
        @. @views U[:, begin] = U[:, begin + 1]
    elseif left_or_right == :right
        @. @views U[:, end] = U[:, end - 1]
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc!(U, bc::Dirichlet_ionbohm, left_or_right::Symbol, ϵ0::Float64, mᵢ::Float64)
    u_bohm = sqrt(2/3*e*ϵ0/mᵢ)

    if left_or_right == :left
        U[1, begin] = bc.state[1] - U[3, begin + 1]/150.0

        # Ion bohm condition, ui ≤ -u_bohm
        boundary_flux = U[3, begin+1]
        boundary_velocity = min(-u_bohm, boundary_flux / U[2, begin+1])
        boundary_density = boundary_flux / boundary_velocity
        U[2, begin] = boundary_density #3*10e17*HallThruster.Xenon.m #boundary_density #U[2, begin+1]
        U[3, begin] = U[3, begin+1]

    elseif left_or_right == :right

        U[1, end] = bc.state[1] + U[2, end-1]

        boundary_flux = U[3, end-1]
        boundary_velocity = max(U[3, end-1] / U[2, end-1], u_bohm)
        boundary_density = boundary_flux / boundary_velocity
        U[2, end] = boundary_density
        U[3, end] = boundary_density * boundary_velocity

    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc!(U, bc::Neumann_ionbohm, left_or_right::Symbol, ϵ0::Float64, mᵢ::Float64)
    u_bohm = sqrt(2/3*e*ϵ0/mᵢ)

    if left_or_right == :left
        @. @views U[1:2, begin] = U[1:2, begin+1]

        # Ion bohm condition, ui ≤ -u_bohm
        boundary_flux = U[3, begin+1]
        boundary_velocity = min(-u_bohm, boundary_flux / U[2, begin+1])
        boundary_density = boundary_flux / boundary_velocity
        U[2, begin] = boundary_density
        U[3, begin] = U[3, begin+1]

    elseif left_or_right == :right
        @. @views U[1:2, end] = U[1:2, end-1]

        # make sure ui[end] ≥ u_bohm
        U[3, end] = max(U[3, end-1], u_bohm * U[2, end])
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc_electron!(U, bc::Dirichlet_energy, left_or_right::Symbol, index::NamedTuple)
    if left_or_right == :left
        @views U[index.nϵ, begin] = bc.state
    elseif left_or_right == :right
        @views U[index.nϵ, end] = bc.state
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc_electron!(U, bc::Dirichlet_energy_upd_ne, left_or_right::Symbol, index::NamedTuple)
    if left_or_right == :left
        @views U[index.nϵ, begin] = bc.int_energy*U[2, begin]/HallThruster.Xenon.m
    elseif left_or_right == :right
        @views U[index.nϵ, end] = bc.int_energy*U[2, end]/HallThruster.Xenon.m
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end

function apply_bc_electron!(U, bc::Neumann_energy, left_or_right::Symbol, index::NamedTuple)
    if left_or_right == :left
        @views U[index.nϵ, begin] = U[index.nϵ, begin+1]
    elseif left_or_right == :right
        @views U[index.nϵ, end] = U[index.nϵ, end-1]
    else
        throw(ArgumentError("left_or_right must be either :left or :right"))
    end
end