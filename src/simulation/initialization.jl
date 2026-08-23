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
    allowed_charges = propellant.allowed_charges
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

    ion_density_function(z, Z) = begin
        base = mi * scaling_factor * (
            ni_min + (ni_max - ni_min) * exp(-(((z - z0) - ni_center) / ni_width)^2)
        )

        #negative ions are assumed to have a much lower initial density
        if Z > 0
            base / Z^2
        elseif Z < 0
            1.0e-3 * base / abs(Z)^2
        end
    end
    # The ion velocity curve combines a quadratic part upstream and a linear part downstream
    # Negative Ions are assumed to initially have 0 average velocity to avoid initial instability
    anode_velocity(Z) = if Z > 0
        -sqrt(Z * e * anode_Tev / mi)
    elseif Z < 0
        0.0
    end

    cathode_velocity(Z) = if Z > 0
        sqrt(2 * Z * e * (discharge_voltage) / mi)
    elseif Z < 0
        0.0
    end

    scale(Z) = 2 / 3 * (cathode_velocity(Z) - anode_velocity(Z))
    ion_velocity_f1(z, Z) = anode_velocity(Z) + scale(Z) * ((z - z0) / L_ch)^2
    ion_velocity_f2(z, Z) = lerp(z, z0 + L_ch, z1, ion_velocity_f1(z0 + L_ch, Z), cathode_velocity(Z))

    ion_velocity_function(z, Z) = if (z - z0) < L_ch
        ion_velocity_f1(z, Z)
    else
        ion_velocity_f2(z, Z)
    end

    # Neutral density at inlet
    ρn_0 = inlet_neutral_density(propellant, thruster.geometry.channel_area)
    # add recombined neutrals
    for Z in allowed_charges
        ρn_0 -= ion_velocity_function(0.0, Z) * ion_density_function(0.0, Z) / un
    end

    # Neutral density at outlet
    ρn_1 = 0.01 * ρn_0
    # Tanh function steps between inlet and outlet densities
    neutral_function(z) = smooth_if(z - z0, L_ch / 2, ρn_0, ρn_1, L_ch / 6)

    # Only the ground-state neutral carries the neutral density at t=0; excited states
    # start empty and fill via excitation reactions.
    @inbounds for (i, fluid) in enumerate(fluids.continuity)
        if i == 1
            @. fluid.density = neutral_function(grid.cell_centers)
        else
            @. fluid.density = 1.0e-10 * neutral_function(grid.cell_centers)  # floor, ~0
        end
    end

    # Excited ions start near zero to avoid duplicating the initial charge density.
    @inbounds for fluid in fluids.isothermal
        Z = fluid.species.Z
        excitation_fraction = is_excited(fluid.species) ? 1.0e-10 : 1.0
        @. fluid.density = excitation_fraction * ion_density_function(grid.cell_centers, Z)
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
    params.cache.ne .= 0.0
    for fluid in params.fluid_array
        @. params.cache.ne += fluid.species.Z * fluid.density / fluid.species.element.m
    end

    return
end

function initialize_electrons_default!(params, anode_Tev, cathode_Tev, discharge_voltage; max_electron_temperature = -1.0)
    (; grid, cache, thruster) = params
    L_ch = thruster.geometry.channel_length
    z0 = grid.edges[1]
    z1 = grid.edges[end]

    # Electron temperature
    Te_baseline(z) = lerp(z, z0, z1, anode_Tev, cathode_Tev)
    base_Te = 0.5 * (anode_Tev + cathode_Tev)
    Te_max = max_electron_temperature > 0.0 ? max_electron_temperature : discharge_voltage / 10
    Te_width = L_ch / 3

    # Gaussian Te profile
    energy_function(z) = 1.5 * (Te_baseline(z) + (Te_max - base_Te) * exp(-(((z - z0) - L_ch) / Te_width)^2))

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
    update_pressure!(params.cache.pe, params.cache.nϵ, config.LANDMARK)
    update_pressure_gradient!(params.cache.∇pe, params.cache.pe, params.grid.cell_centers)

    return
end

"""
    $(TYPEDSIGNATURES)
Initialize fluid containers and other plasma variables form a restart
"""
function initialize_from_restart!(params, restart_file::String)
    restart = JSON.parsefile(restart_file)
    _validate_serialization_version(restart, "Restart file $(restart_file)")

    if haskey(restart, "output")
        restart = restart["output"]
    end

    if haskey(restart, "frames")
        frame = restart["frames"][end]
    elseif haskey(restart, "average")
        frame = restart["average"]
    else
        throw(ArgumentError("Restart file $(restart_file) has no key `frames` or `average`."))
    end

    return initialize_from_restart!(params, frame)
end

function initialize_from_restart!(params, frame)
    (; grid, cache, propellants, fluids_by_propellant, species_energies_eV) = params
    z = grid.cell_centers
    z_frame = _restart_field(frame, "z", "restart frame")

    _restart_collection(frame, "neutrals", "restart frame")
    _restart_collection(frame, "ions", "restart frame")
    excited_states = get(frame, "excited_states", nothing)
    if !isnothing(excited_states) && !(excited_states isa AbstractDict)
        throw(ArgumentError("Restart field `excited_states` must be a dictionary."))
    end

    # Restore every configured ground or excited fluid from the structured species output.
    for (propellant, fluids) in zip(propellants, fluids_by_propellant)
        gas_symbol = string(propellant.gas.formula)

        for fluid in fluids.continuity
            species = fluid.species
            state = if is_excited(species)
                isnothing(excited_states) ? nothing :
                    get(excited_states, string(species.symbol), nothing)
            else
                neutrals = frame["neutrals"]
                haskey(neutrals, gas_symbol) || throw(ArgumentError(
                    "Restart output has no ground-state $(gas_symbol) neutral."
                ))
                neutrals[gas_symbol]
            end
            if isnothing(state)
                fill!(fluid.density, 0.0)
                continue
            end
            _validate_restart_species!(
                state, species, species_energies_eV[species.symbol],
            )
            number_density = _restart_field(
                state, "n", "restart state $(species.symbol)", length(z_frame),
            )
            fluid.density .= LinearInterpolation(
                z_frame, number_density .* species.element.m
            ).(z)
        end

        for fluid in fluids.isothermal
            species = fluid.species
            state = if is_excited(species)
                isnothing(excited_states) ? nothing :
                    get(excited_states, string(species.symbol), nothing)
            else
                ions = frame["ions"]
                haskey(ions, gas_symbol) || throw(ArgumentError(
                    "Restart output has no ground-state $(gas_symbol) ions."
                ))
                ion_states = ions[gas_symbol]
                ion_states isa AbstractVector || throw(ArgumentError(
                    "Restart ground-state $(gas_symbol) ions must be an array."
                ))
                index = findfirst(ion_states) do ion
                    ion isa AbstractDict && get(ion, "Z", nothing) == species.Z
                end
                isnothing(index) && throw(ArgumentError(
                    "Restart output has no $(species.Z)-charged ground-state $(gas_symbol) ions."
                ))
                ion_states[index]
            end
            if isnothing(state)
                fill!(fluid.density, 0.0)
                fill!(fluid.momentum, 0.0)
                continue
            end
            _validate_restart_species!(
                state, species, species_energies_eV[species.symbol],
            )
            number_density = _restart_field(
                state, "n", "restart state $(species.symbol)", length(z_frame),
            )
            number_flux = _restart_field(
                state, "nu", "restart state $(species.symbol)", length(z_frame),
            )
            mass = species.element.m
            fluid.density .= LinearInterpolation(z_frame, number_density .* mass).(z)
            fluid.momentum .= LinearInterpolation(z_frame, number_flux .* mass).(z)
        end
    end

    ne = _restart_field(frame, "ne", "restart frame", length(z_frame))
    cache.ne .= LinearInterpolation(z_frame, ne).(z)

    # load electron properties
    Te = LinearInterpolation(
        z_frame, _restart_field(frame, "Tev", "restart frame", length(z_frame)),
    ).(z)
    phi = LinearInterpolation(
        z_frame, _restart_field(frame, "potential", "restart frame", length(z_frame)),
    ).(z)
    E = LinearInterpolation(
        z_frame, _restart_field(frame, "E", "restart frame", length(z_frame)),
    ).(z)

    @. cache.nϵ = 1.5 * cache.ne * Te
    @. cache.Tev = Te
    @. cache.∇ϕ = -E
    @. cache.ϕ = phi

    return nothing
end

function _restart_collection(frame, key, context)
    haskey(frame, key) || throw(ArgumentError("$(context) has no `$(key)` field."))
    collection = frame[key]
    collection isa AbstractDict || throw(ArgumentError(
        "Restart field `$(key)` must be a dictionary."
    ))
    return collection
end

function _restart_field(container, key, context, expected_length = nothing)
    haskey(container, key) || throw(ArgumentError("$(context) has no `$(key)` field."))
    values = container[key]
    values isa AbstractVector || throw(ArgumentError(
        "Field `$(key)` in $(context) must be an array."
    ))
    if !isnothing(expected_length) && length(values) != expected_length
        throw(ArgumentError(
            "Field `$(key)` in $(context) has $(length(values)) values; " *
                "expected $(expected_length)."
        ))
    end
    all(value -> value isa Real && isfinite(value), values) || throw(ArgumentError(
        "Field `$(key)` in $(context) must contain only finite numbers."
    ))
    return values
end

function _validate_restart_species!(state, species, expected_energy_eV)
    state isa AbstractDict || throw(ArgumentError(
        "Restart state $(species.symbol) must be a dictionary."
    ))
    for (key, expected) in (("Z", species.Z), ("excited_level", species.excited_level))
        if haskey(state, key) && state[key] != expected
            throw(ArgumentError(
                "Restart state $(species.symbol) has $(key)=$(state[key]); expected $(expected)."
            ))
        end
    end
    if haskey(state, "energy_eV")
        restart_energy_eV = state["energy_eV"]
        restart_energy_eV isa Real && isfinite(restart_energy_eV) || throw(ArgumentError(
            "Restart state $(species.symbol) has a non-numeric or non-finite energy."
        ))
        if !isapprox(
                restart_energy_eV, expected_energy_eV;
                atol = EXCITATION_ENERGY_MERGE_TOLERANCE_EV, rtol = 0,
            )
            throw(ArgumentError(
                "Restart state $(species.symbol) has energy $(restart_energy_eV) eV; " *
                    "the active chemistry uses $(expected_energy_eV) eV."
            ))
        end
    end
    return nothing
end
