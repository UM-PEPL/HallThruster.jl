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

function initialize_gas!(propellant, fluids, params; max_ion_density, min_ion_density, anode_Tev, discharge_voltage)
    (; grid, thruster) = params
    mi = propellant.gas.m
    ncharge = propellant.max_charge
    flow_rate = propellant.flow_rate_kg_s
    un = propellant.velocity_m_s

    L_ch = thruster.geometry.channel_length
    z0 = grid.cell_centers[1]
    z1 = grid.cell_centers[end]

    ni_center = L_ch / 2
    ni_width = L_ch / 3
    ni_min = min_ion_density
    ni_max = max_ion_density

    # Scale density up for high flow rates and down for high voltages
    scaling_factor = sqrt(discharge_voltage / 300) * (flow_rate / 5.0e-6)
    ion_density_function(z, Z) = mi * scaling_factor * (
        ni_min + (ni_max - ni_min) * exp(-(((z - z0) - ni_center) / ni_width)^2)
    ) / Z^2

    # The ion velocity curve combines a quadratic part upstream and a linear part downstream
    bohm_velocity = Z -> -sqrt(Z * e * anode_Tev / mi)
    final_velocity = Z -> sqrt(2 * Z * e * discharge_voltage / mi)
    scale(Z) = 2 / 3 * (final_velocity(Z) - bohm_velocity(Z))
    ion_velocity_f1(z, Z) = bohm_velocity(Z) + scale(Z) * ((z - z0) / L_ch)^2
    ion_velocity_f2(z, Z) = lerp(z, z0 + L_ch, z1, ion_velocity_f1(L_ch, Z), final_velocity(Z))

    ion_velocity_function(z, Z) = if (z - z0) < L_ch
        ion_velocity_f1(z, Z)
    else
        ion_velocity_f2(z, Z)
    end

    # Neutral density at inlet
    ρn_0 = inlet_neutral_density(propellant, thruster.geometry.channel_area)
    # add recombined neutrals
    for Z in 1:ncharge
        ρn_0 -= ion_velocity_function(0.0, Z) * ion_density_function(0.0, Z) / un
    end

    # Neutral density at outlet
    ρn_1 = 0.01 * ρn_0
    # Tanh function steps between inlet and outlet densities
    neutral_function(z) = smooth_if(z - z0, L_ch / 2, ρn_0, ρn_1, L_ch / 6)

    # Fill the fluid containers
    @inbounds for fluid in fluids.continuity
        @. fluid.density = neutral_function(grid.cell_centers)
    end

    @inbounds for fluid in fluids.isothermal
        Z = fluid.species.Z
        @. fluid.density = ion_density_function(grid.cell_centers, Z)
        @. fluid.momentum = fluid.density * ion_velocity_function(grid.cell_centers, Z)
    end

    return
end

function initialize_heavy_species_default!(params; kwargs...)
    # Initialize each propellant species
    for (propellant, fluids) in zip(params.propellants, params.fluids_by_propellant)
        initialize_gas!(propellant, fluids, params; kwargs...)
    end

    # Compute the electron number density
    for fluid in params.fluid_array
        @. params.cache.ne += fluid.species.Z * fluid.density
    end

    return
end

function initialize_electrons_default!(params, anode_Tev, cathode_Tev, discharge_voltage; max_electron_temperature = -1.0)
    (; grid, cache, min_Te, thruster) = params
    L_ch = thruster.geometry.channel_length
    z0 = grid.cell_centers[1]
    z1 = grid.cell_centers[end]

    # Electron temperature
    Te_baseline(z) = lerp(z, z0, z1, anode_Tev, cathode_Tev)
    Te_max = max_electron_temperature > 0.0 ? max_electron_temperature : discharge_voltage / 10
    Te_width = L_ch / 3

    # Gaussian Te profile
    energy_function(z) = 3 / 2 * (Te_baseline(z) + (Te_max - min_Te) * exp(-(((z - z0) - L_ch) / Te_width)^2))

    for (i, z) in enumerate(grid.cell_centers)
        cache.nϵ[i] = cache.ne[i] * energy_function(z)
        cache.Tev[i] = energy_function(z) / 1.5
    end

    return
end

function initialize!(params, config, init::DefaultInitialization)
    (; max_electron_temperature, min_ion_density, max_ion_density) = init
    (; anode_Tev, cathode_Tev, discharge_voltage) = config

    initialize_heavy_species_default!(params; discharge_voltage, anode_Tev, max_ion_density, min_ion_density)
    initialize_electrons_default!(params, anode_Tev, cathode_Tev, discharge_voltage; max_electron_temperature)
    return
end

"""
    $(TYPEDSIGNATURES)
Initialize fluid containers and other plasma variables form a restart
"""
function initialize_from_restart!(params, restart_file::String)
    # TODO: multiple propellants
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
    # TODO: multiple propellants
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
