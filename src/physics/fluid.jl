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
        R = R0 / species.element.M
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

"""
$(TYPEDEF)

All `FluidContainer`s in a simulation. `continuity` holds the ground-state neutral first,
then one fluid per tracked excited state; `isothermal` holds one per ion charge state.
"""
struct FluidContainerSet
    continuity::Vector{FluidContainer}
    isothermal::Vector{FluidContainer}
end

"""The ground-state neutral fluid (always the first continuity fluid)."""
ground_neutral(fluids::FluidContainerSet) = fluids.continuity[1]

"""All excited-state neutral fluids in the set."""
excited_fluids(fluids::FluidContainerSet) =
    [f for f in fluids.continuity if is_excited(f.species)]

"""
Excitation levels with per-channel rate data for `gas`; each gets its own neutral fluid.
Stub: the real implementation belongs in the reaction-loading code, scanning the rate
tables on disk. Empty preserves the lumped-excitation behavior.
"""
supported_excitation_levels(::Gas) = Int[]

function allocate_fluids(p::Propellant, ncells; excited_levels = supported_excitation_levels(p.gas))
    # excited neutrals advect with the ground-state velocity and temperature
    continuity = [
        FluidContainer(
            _ContinuityOnly, p.gas(0, n), ncells;
            vel = p.velocity_m_s, temp = p.temperature_K,
        )
        for n in [0; sort!(collect(excited_levels))]
    ]

    isothermal = [
        FluidContainer(_IsothermalEuler, p.gas(Z), ncells; temp = p.ion_temperature_K) for Z in p.allowed_charges
    ]
    return FluidContainerSet(continuity, isothermal)
end
