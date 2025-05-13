abstract type AbstractFluid end

struct ContinuityFluid <: AbstractFluid 
    # Conservative variables
    density::Vector{Float64}

    # Time derivative of conservative variables
    dens_ddt::Vector{Float64}

    # Edge states
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}

    # Fluxes
    flux_dens::Vector{Float64}

    # Data
    wave_speed::Array{Float64, 0}
    max_timestep::Array{Float64, 0}
    species::Species
    sound_speed::Float64
    velocity::Float64

    function ContinuityFluid(species, vel, temp, num_cells)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables and time derivatives
            zeros(num_cells+2),
            zeros(num_cells+2),

            # Edge states
            zeros(num_cells+1), zeros(num_cells+1),

            # Fluxes
            zeros(num_cells+1),

            # Data
            fill(max(abs(vel+a), abs(vel-a))), fill(0.0), species, a, vel
        )
    end
end

struct IsothermalFluid <: AbstractFluid
    # Conservative variables
    density::Vector{Float64}
    momentum::Vector{Float64}

    # Time deriviatives
    dens_ddt::Vector{Float64}
    mom_ddt::Vector{Float64}

    # Edge states
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}
    mom_L::Vector{Float64}
    mom_R::Vector{Float64}    

    # Fluxes
    flux_dens::Vector{Float64}
    flux_mom::Vector{Float64}

    # Data
    wave_speed::Array{Float64, 0}
    max_timestep::Array{Float64, 0}
    species::Species
    sound_speed::Float64

    function IsothermalFluid(species, temp, num_cells)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables and time derivatives
            zeros(num_cells+2), zeros(num_cells+2),
            zeros(num_cells+2), zeros(num_cells+2),

            # Edge states
            zeros(num_cells+1), zeros(num_cells+1),
            zeros(num_cells+1), zeros(num_cells+1),

            # Fluxes
            zeros(num_cells+1), zeros(num_cells+1),

            # Data
            fill(a), fill(0.0), species, a
        )
    end
end

function allocate_fluids(propellant, ncharge, ncells, u_neutral, T_neutral, T_ion)
    continuity = [
        ContinuityFluid(propellant(0), u_neutral, T_neutral, ncells)
    ]

    isothermal= [
        IsothermalFluid(propellant(Z), T_ion, ncells) for Z in 1:ncharge 
    ]
    return (;continuity, isothermal)
end

function _to_state_vector!(U, continuity, isothermal)
    @. @views U[1, :] = continuity[1].density

    for (i, fluid) in eachindex(isothermal)
        @. @views U[2*i, :] = fluid.density
        @. @views U[2*i+1, :] = fluid.momentum
    end
end

function _from_state_vector!(continuity, isothermal, U)
    @. @views continuity[1].density = U[1, :]

    for (i, fluid) in enumerate(isothermal)
        @. @views fluid.density = U[2*i, :]
        @. @views fluid.momentum = U[2*i+1, :]
    end
end

@enum ConservationLawType begin
    _ContinuityOnly
    _IsothermalEuler
    _EulerEquations
end
Base.@kwdef struct Fluid
    species::Species
    type::ConservationLawType
    nvars::Int
    u::Union{Float64,Nothing} = nothing
    T::Union{Float64,Nothing} = nothing
    a::Union{Float64,Nothing} = nothing
end

@inline nvars(f::Fluid) = f.nvars

function Fluid(s; u = nothing, T = nothing)
    if isnothing(u) && isnothing(T)
        return EulerEquations(s)
    elseif isnothing(u)
        return Fluid(s, T)
    else
        return Fluid(s, u, T)
    end
end

function ContinuityOnly(s; u, T)
    a = √(s.element.γ * kB * T / s.element.m)
    return Fluid(species = s, type=_ContinuityOnly, nvars=1, u=Float64(u), T=Float64(T), a = Float64(a))
end
ContinuityOnly(s, u, T) = ContinuityOnly(s; u, T)
Fluid(s, u, T) = ContinuityOnly(s, u, T)

function IsothermalEuler(s; T)
    a = √(s.element.γ * kB * T / s.element.m)
    return Fluid(species = s, type=_IsothermalEuler, nvars=2, T=Float64(T), a = Float64(a))
end
IsothermalEuler(s, T) = IsothermalEuler(s; T)
Fluid(s, T) = IsothermalEuler(s, T)

function EulerEquations(s)
    return Fluid(species = s, type=_EulerEquations, nvars=3)
end

function ranges(fluids)
    fluid_ranges = fill(1:1, length(fluids))
    start_ind = 1
    last_ind = 1
    for (i, f) in enumerate(fluids)
        nf = nvars(f)
        last_ind = start_ind - 1 + nf
        fluid_ranges[i] = start_ind:last_ind
        start_ind = last_ind + 1
    end

    return fluid_ranges
end
