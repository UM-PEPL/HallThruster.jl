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
    """Cell-centered primitive velocity cache, refreshed before each derivative evaluation"""
    vel_prim::Vector{Float64}
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}
    vel_L::Vector{Float64}
    vel_R::Vector{Float64}
    flux_dens::Vector{Float64}
    flux_mom::Vector{Float64}
    """Edge-local wave speed used by continuity-only fluids"""
    wave_speed::Vector{Float64}
    """Maximum permissable timestep for this species"""
    max_timestep::Array{Float64, 0}
    """The `Species` whose properties are stored in this struct"""
    species::Species
    """The sound speed for this species"""
    sound_speed::Float64
    """The type of species (_ContinuityOnly or _IsothermalEuler)"""
    type::ConservationLawType

    function FluidContainer(type, species, grid; temp, vel = 0.0)
        num_cells = length(grid.cell_centers)
        num_edges = length(grid.edges)
        R = R0 / species.element.M
        γ = species.element.γ

        if type == _ContinuityOnly
            cell_velocity = vel.(grid.cell_centers)
            edge_velocity = vel.(grid.edges)
            edge_temperature = temp.(grid.edges)
            edge_sound_speed = @. sqrt(γ * R * edge_temperature)
            wave_speed = @. abs(edge_velocity) + edge_sound_speed
            sound_speed = maximum(edge_sound_speed)
            vel_L = copy(edge_velocity)
            vel_R = copy(edge_velocity)
        else
            cell_velocity = zeros(num_cells)
            wave_speed = zeros(num_edges)
            sound_speed = sqrt(γ * R * temp)
            vel_L = zeros(num_edges)
            vel_R = zeros(num_edges)
        end

        return new(
            # Conservative variables, caches, and time derivatives
            zeros(num_cells), zeros(num_cells),
            zeros(num_cells), zeros(num_cells),
            zeros(num_cells), zeros(num_cells),
            cell_velocity,

            # Edge states
            zeros(num_edges), zeros(num_edges),
            vel_L, vel_R,

            # Fluxes
            zeros(num_edges), zeros(num_edges),

            # Data
            wave_speed, fill(0.0), species, sound_speed, type
        )
    end
end

"""
$(TYPEDEF)

Collection of all `FluidContainer`s, grouped by conservation law type. `continuity` holds
the ground-state neutral first, then one fluid per excited state; `isothermal` holds one
fluid per ion charge state.
"""
struct FluidContainerSet
    continuity::Vector{FluidContainer}
    isothermal::Vector{FluidContainer}
end

"""
    ground_neutral(fluids::FluidContainerSet) -> FluidContainer
The ground-state neutral fluid (the first continuity fluid).
"""
ground_neutral(fluids::FluidContainerSet) = fluids.continuity[1]

"""
    excited_fluids(fluids::FluidContainerSet) -> Vector{FluidContainer}
All excited-state neutral fluids in the set.
"""
excited_fluids(fluids::FluidContainerSet) =
    [f for f in fluids.continuity if is_excited(f.species)]

function allocate_fluids(p::Propellant, grid; excited_levels = p.excited_levels)
    # Ground state first, then one fluid per excited level, all advecting with the neutral flow
    continuity = [
        FluidContainer(
            _ContinuityOnly, p.gas(0, excited_level), grid;
            vel = p.velocity_m_s, temp = p.temperature_K,
        )
            for excited_level in [0; sort!(collect(excited_levels))]
    ]

    isothermal = [
        FluidContainer(
            _IsothermalEuler, p.gas(Z, excited_level), grid;
            temp = p.ion_temperature_K,
        )
            for Z in p.allowed_charges
            for excited_level in [0; get(p.excited_ion_levels, Z, Int[])]
    ]
    return FluidContainerSet(continuity, isothermal)
end
