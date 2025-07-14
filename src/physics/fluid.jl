@enum ConservationLawType begin
    _ContinuityOnly
    _IsothermalEuler
    _EulerEquations
end

"""
$(TYPEDEF)

Struct containing necessary internal states and caches for solving the heavy species fluid equations and for interfacing with the electron solver.

# Fields
$(TYPEDFIELDS)
"""
struct FluidContainer
    """Mass density in kg/m^3"""
    density::Vector{Float64}
    """Momentum density in kg/m^2 s"""
    momentum::Vector{Float64}
    dens_ddt::Vector{Float64}
    mom_ddt::Vector{Float64}
    dens_cache::Vector{Float64}
    mom_cache::Vector{Float64}
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}
    mom_L::Vector{Float64}
    mom_R::Vector{Float64}
    flux_dens::Vector{Float64}
    flux_mom::Vector{Float64}
    """Maximum wave speed for this species"""
    wave_speed::Array{Float64, 0}
    """Maximum permissable timestep for this species"""
    max_timestep::Array{Float64, 0}
    """The `Species` whose properties are stored in this struct"""
    species::Species
    """The sound speed for this species"""
    sound_speed::Float64
    """For neutral species, the constant advection speed of this species"""
    const_velocity::Float64
    """The type of species (_ContinuityOnly or _IsothermalEuler)"""
    type::ConservationLawType

    function FluidContainer(type, species, num_cells; temp, vel = 0.0)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables, caches, and time derivatives
            zeros(num_cells + 2), zeros(num_cells + 2),
            zeros(num_cells + 2), zeros(num_cells + 2),
            zeros(num_cells + 2), zeros(num_cells + 2),

            # Edge states
            zeros(num_cells + 1), zeros(num_cells + 1),
            zeros(num_cells + 1), zeros(num_cells + 1),

            # Fluxes
            zeros(num_cells + 1), zeros(num_cells + 1),

            # Data
            fill(max(abs(vel + a), abs(vel - a))), fill(0.0), species, a, vel, type
        )
    end
end

function temperature(f::FluidContainer)
    R = f.species.element.R
    γ = f.species.element.γ
    return f.sound_speed^2 / (γ * R)
end

function allocate_fluids(p::Propellant, ncells)
    continuity = [
        FluidContainer(_ContinuityOnly, p.gas(0), ncells; vel = p.velocity_m_s, temp = p.temperature_K),
    ]

    isothermal = [
        FluidContainer(_IsothermalEuler, p.gas(Z), ncells; temp = p.ion_temperature_K) for Z in 1:p.max_charge
    ]
    return (; continuity, isothermal)
end
