abstract type InitialCondition end

Base.@kwdef struct DefaultInitialization <: InitialCondition
    max_electron_temperature::Float64 = -1.0
    min_ion_density::Float64 = 2.0e17
    max_ion_density::Float64 = 1.0e18
end

#=============================================================================
 Serialization
==============================================================================#
Serialization.SType(::Type{T}) where {T <: InitialCondition} = Serialization.TaggedUnion()
Serialization.options(::Type{T}) where {T <: InitialCondition} = (; DefaultInitialization)

#=============================================================================
Definitions
==============================================================================#
initialize!(params, config) = initialize!(params, config, config.initial_condition)

function initialize!(_, _, model::InitialCondition)
    throw(ArgumentError("Function HallThruster.initialize!(params, model::$(typeof(model)) not yet implemented. For InitialCondition types other than DefaultInitialization(), this must be defined by the user!"))
end

function initialize_heavy_species_default!(params, config; kwargs...)
    (; anode_Tev, domain, discharge_voltage, propellants) = config

    ρn = inlet_neutral_density(config)
    un = config.neutral_velocity
    return initialize_heavy_species_default!(
        params, anode_Tev, domain, discharge_voltage,
        propellants[1].flow_rate_kg_s, ρn, un; kwargs...,
    )
end

function initialize_heavy_species_default!(
        params, anode_Tev, domain, discharge_voltage, anode_mass_flow_rate, ρn_0, un;
        min_ion_density = 2.0e17, max_ion_density = 1.0e18,
    )

    (; grid, cache, thruster, propellants) = params
    mi = propellants[1].gas.m
    ncharge = propellants[1].max_charge

    L_ch = thruster.geometry.channel_length
    z0 = domain[1]

    ni_center = L_ch / 2
    ni_width = L_ch / 3
    ni_min = min_ion_density
    ni_max = max_ion_density
    scaling_factor = sqrt(discharge_voltage / 300) * (anode_mass_flow_rate / 5.0e-6)
    ion_density_function = (z, Z) -> mi * scaling_factor *
        (
        ni_min +
            (ni_max - ni_min) *
            exp(-(((z - z0) - ni_center) / ni_width)^2)
    ) / Z^2

    bohm_velocity = Z -> -sqrt(Z * e * anode_Tev / mi)

    final_velocity = Z -> sqrt(2 * Z * e * discharge_voltage / mi)
    scale(Z) = 2 / 3 * (final_velocity(Z) - bohm_velocity(Z))
    ion_velocity_f1(z, Z) = bohm_velocity(Z) + scale(Z) * ((z - z0) / L_ch)^2
    function ion_velocity_f2(z, Z)
        return lerp(z, z0 + L_ch, domain[2], ion_velocity_f1(L_ch, Z), final_velocity(Z))
    end

    ion_velocity_function = (z, Z) -> if z - z0 < L_ch
        ion_velocity_f1(z, Z)
    else
        ion_velocity_f2(z, Z)
    end

    # add recombined neutrals
    for Z in 1:(ncharge)
        ρn_0 -= ion_velocity_function(0.0, Z) * ion_density_function(0.0, Z) / un
    end

    # Beam neutral density at outlet
    ρn_1 = 0.01 * ρn_0
    neutral_function = z -> smooth_if(z - z0, L_ch / 2, ρn_0, ρn_1, L_ch / 6)

    # Electron density
    number_density_function = z -> sum(
        Z * ion_density_function(z, Z) / mi
            for Z in 1:ncharge
    )

    # Fill the fluid containers
    @inbounds for fluid in params.fluid_containers.continuity
        @. fluid.density = neutral_function(grid.cell_centers)
    end

    @inbounds for fluid in params.fluid_containers.isothermal
        Z = fluid.species.Z
        @. fluid.density = ion_density_function(grid.cell_centers, Z)
        @. fluid.momentum = fluid.density * ion_velocity_function(grid.cell_centers, Z)
    end

    @. cache.ne = number_density_function(grid.cell_centers)

    return nothing
end

function initialize_electrons_default!(
        params, anode_Tev, cathode_Tev, domain,
        discharge_voltage; max_electron_temperature = -1.0,
    )
    (; grid, cache, min_Te, thruster) = params
    L_ch = thruster.geometry.channel_length
    z0 = domain[1]
    # Electron temperature
    Te_baseline = z -> lerp(z, domain[1], domain[2], anode_Tev, cathode_Tev)
    Te_max = max_electron_temperature > 0.0 ? max_electron_temperature :
        discharge_voltage / 10
    Te_width = L_ch / 3

    # Gaussian Te profile
    energy_function = z -> 3 / 2 * (
        Te_baseline(z) +
            (Te_max - min_Te) * exp(-(((z - z0) - L_ch) / Te_width)^2)
    )

    for (i, z) in enumerate(grid.cell_centers)
        cache.nϵ[i] = cache.ne[i] * energy_function(z)
        cache.Tev[i] = energy_function(z) / 1.5
    end

    return nothing
end

function initialize!(params, config, init::DefaultInitialization)
    (; max_electron_temperature, min_ion_density, max_ion_density) = init
    (; anode_Tev, cathode_Tev, domain, discharge_voltage, propellants) = config
    ρn = inlet_neutral_density(config)
    un = propellants[1].velocity_m_s
    initialize_heavy_species_default!(
        params, anode_Tev, domain, discharge_voltage,
        propellants[1].flow_rate_kg_s, ρn, un; max_ion_density, min_ion_density,
    )
    return initialize_electrons_default!(
        params, anode_Tev, cathode_Tev, domain, discharge_voltage; max_electron_temperature,
    )
end

"""
    $(TYPEDSIGNATURES)
Initialize fluid containers and other plasma variables form a restart
"""
function initialize_from_restart!(params, restart_file::String)
    restart = JSON3.read(read(restart_file))

    if haskey(restart, "output")
        restart = restart.output
    end

    if haskey(restart, "frames")
        frame = restart.frames[end]
    elseif haskey(restart, "average")
        frame = restart.average
    else
        throw(ArgumentError("Restart file $(restart_file) has no key `frames` or `average`."))
    end

    return initialize_from_restart!(params, frame)
end

function initialize_from_restart!(params, frame)
    (; grid, cache, propellants) = params
    mi = propellants[1].gas.m
    ncharge = propellants[1].max_charge
    ncharge_restart = length(frame.ni)

    # load ion properties, interpolated from restart grid to grid in params
    z = grid.cell_centers

    nn = LinearInterpolation(frame.z, frame.nn .* mi).(z)
    params.fluid_containers.continuity[1].density .= nn

    for Z in 1:min(ncharge, ncharge_restart)
        fluid = params.fluid_containers.isothermal[Z]
        fluid.density .= LinearInterpolation(frame.z, frame.ni[Z] .* mi).(z)
        fluid.momentum .= LinearInterpolation(frame.z, frame.niui[Z] .* mi).(z)
    end

    cache.ne .= LinearInterpolation(frame.z, frame.ne).(z)

    # load electron properties
    Te = LinearInterpolation(frame.z, frame.Tev).(z)
    phi = LinearInterpolation(frame.z, frame.potential).(z)
    E = LinearInterpolation(frame.z, frame.E).(z)

    @. cache.nϵ = 1.5 * cache.ne * Te
    @. cache.Tev = Te
    @. cache.∇ϕ = -E
    @. cache.ϕ = phi

    return nothing
end
