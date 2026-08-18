@public time_average, thrust, discharge_current
@public ion_current, electron_current
@public voltage_eff, mass_eff, current_eff, divergence_eff, anode_eff

"""
    $(TYPEDSIGNATURES)
Average a `Solution` over time, starting at time `start_time`.
Return a `Solution` object with a single frame containing the averaged simulation properties
"""
function time_average(sol::Solution, start_time)
    start_time = convert_to_float64(start_time, units(:s))
    start_frame = findfirst(>=(start_time), sol.t)
    return time_average(sol, start_frame)
end

"""
    $(TYPEDSIGNATURES)
Average a `Solution` over time, starting at frame `start_frame`.
Return a `Solution` object with a single frame containing the averaged simulation properties
"""
function time_average(sol::Solution, start_frame::Integer = 1)
    avg_frame = deepcopy(sol.frames[end])
    fields = fieldnames(Frame)

    # Initialize avg to zero
    for f in fields
        avg = getfield(avg_frame, f)
        if f == :anom_variables
            for j in 1:num_anom_variables(sol.config.anom_model)
                avg[j] .= 0.0
            end
        elseif f == :neutrals
            for prop in sol.config.propellants
                symbol = prop.gas.short_name
                avg[symbol].n .= 0.0
                avg[symbol].nu .= 0.0
                avg[symbol].u .= 0.0
            end
        elseif f == :intensities
            for intensity in values(avg)
                intensity .= 0.0
            end
        elseif f == :ions
            for prop in sol.config.propellants
                symbol = prop.gas.short_name
                for ion in avg[symbol]
                    ion.n .= 0.0
                    ion.nu .= 0.0
                    ion.u .= 0.0
                end
            end
        else
            avg .= 0.0
        end
    end

    # Sum over all timesteps to get average
    tstamps = length(sol.t)
    dt = (tstamps - start_frame + 1)
    for i in start_frame:length(sol.t)
        for f in fields
            avg = getfield(avg_frame, f)
            field = getfield(sol.frames[i], f)
            if f == :anom_variables
                for j in 1:num_anom_variables(sol.config.anom_model)
                    avg[j] .+= field[j] ./ dt
                end
            elseif f == :neutrals
                for prop in sol.config.propellants
                    symbol = prop.gas.short_name
                    avg[symbol].n .+= field[symbol].n ./ dt
                    avg[symbol].nu .+= field[symbol].nu ./ dt
                    avg[symbol].u .+= field[symbol].u ./ dt
                end
            elseif f == :intensities
                for (wavelength, intensity) in avg
                    intensity .+= field[wavelength] ./ dt
                end
            elseif f == :ions
                for prop in sol.config.propellants
                    symbol = prop.gas.short_name
                    for (j, ion) in enumerate(avg[symbol])
                        ion.n .+= field[symbol][j].n ./ dt
                        ion.nu .+= field[symbol][j].nu ./ dt
                        ion.u .+= field[symbol][j].u ./ dt
                    end
                end
            else
                avg .+= field ./ dt
            end
        end
    end

    return Solution(
        sol.t[end:end],
        [avg_frame],
        sol.grid,
        sol.config,
        sol.simulation,
        sol.postprocess,
        sol.retcode,
        sol.error,
    )
end

"""
    $(TYPEDSIGNATURES)
Compute the thrust at a specific frame of a `Solution`.
"""
function thrust(sol::Solution, frame::Integer)
    f = sol.frames[frame]
    left_area = f.channel_area[begin]
    right_area = f.channel_area[end]
    thrust = 0.0

    for ions in values(f.ions)
        for ion in ions
            thrust += right_area * ion.m * ion.nu[end]^2 / ion.n[end]
            thrust -= left_area * ion.m * ion.nu[begin]^2 / ion.n[begin]
        end
    end

    # Multiply by sqrt of divergence efficiency to model loss of ions in radial direction
    if (sol.config.apply_thrust_divergence_correction)
        return thrust * sqrt(divergence_eff(sol, frame))
    else
        return thrust
    end
end

"""
    $(TYPEDSIGNATURES)
Compute the thrust at a each frame of a `Solution`.
"""
thrust(sol::Solution) = [thrust(sol, frame) for frame in eachindex(sol.frames)]

"""
    $(TYPEDSIGNATURES)
Compute the discharge current at a specific frame of a `Solution`.
"""
discharge_current(sol::Solution, frame::Integer) = sol.frames[frame].discharge_current[]

"""
    $(TYPEDSIGNATURES)
Compute the discharge current at a each frame of a `Solution`.
"""
discharge_current(sol::Solution) = [s.discharge_current[] for s in sol.frames]

"""
    $(TYPEDSIGNATURES)
Compute the anode efficiency at a specific frame of a `Solution`.
"""
function anode_eff(sol::Solution, frame::Integer)
    T = thrust(sol, frame)
    current = discharge_current(sol, frame)
    Vd = sol.config.discharge_voltage
    mdot_a = sum(prop.flow_rate_kg_s for prop in sol.config.propellants)
    anode_eff = 0.5 * T^2 / current / Vd / mdot_a
    return anode_eff
end

"""
    $(TYPEDSIGNATURES)
Compute the anode efficiency at each frame of a `Solution`.
"""
anode_eff(sol::Solution) = [anode_eff(sol, frame) for frame in eachindex(sol.frames)]

"""
    $(TYPEDSIGNATURES)
Compute the voltage/acceleration efficiency at a specific frame of a `Solution`.
"""
function voltage_eff(sol::Solution, frame::Integer)
    Vd = sol.config.discharge_voltage
    f = sol.frames[frame]
    ni_end = 0.0
    niV_end = 0.0

    for ion in values(f.ions)
        ui = ion[1].u[end]
        ni = ion[1].n[end]
        ni_end += ni
        niV_end += ni * (0.5 * ion[1].m * ui^2 / e)
    end

    V_accel = niV_end / ni_end
    return V_accel / Vd
end

"""
    $(TYPEDSIGNATURES)
Compute the voltage/acceleration efficiency at each frame of a `Solution`.
"""
voltage_eff(sol::Solution) = [voltage_eff(sol, frame) for frame in eachindex(sol.frames)]

"""
    $(TYPEDSIGNATURES)
Compute the divergence efficiency at a specific frame of a `Solution`.
"""
function divergence_eff(sol::Solution, frame::Integer)
    tanδ = sol.frames[frame].tan_div_angle[end]
    δ = atan(tanδ)
    return cos(δ)^2
end

"""
    $(TYPEDSIGNATURES)
Compute the divergence efficiency at each frame of a `Solution`.
"""
divergence_eff(sol::Solution) = [divergence_eff(sol, frame) for frame in eachindex(sol.frames)]

"""
    $(TYPEDSIGNATURES)
Compute the ion current at a specific frame of a `Solution`.
"""
function ion_current(sol::Solution, frame)
    f = sol.frames[frame]
    return f.ji[end] * f.channel_area[end]
end

"""
    $(TYPEDSIGNATURES)
Compute the ion current at each frame of a `Solution`.
"""
ion_current(sol::Solution) = [ion_current(sol, frame) for frame in eachindex(sol.frames)]

"""
    $(TYPEDSIGNATURES)
Compute the electron current at a specific frame of a `Solution`.
"""
electron_current(sol::Solution, frame) = discharge_current(sol, frame) - ion_current(sol, frame)

"""
    $(TYPEDSIGNATURES)
Compute the electron current at each frame of a `Solution`.
"""
function electron_current(sol::Solution)
    return [electron_current(sol, frame) for frame in eachindex(sol.frames)]
end

"""
    $(TYPEDSIGNATURES)
Compute the current/beam utilization efficiency at a specific frame of a `Solution`.
"""
current_eff(sol::Solution, frame) = ion_current(sol, frame) / discharge_current(sol, frame)

"""
    $(TYPEDSIGNATURES)
Compute the current/beam utilization efficiency at each frame of a `Solution`.
"""
current_eff(sol::Solution) = [current_eff(sol, frame) for frame in eachindex(sol.frames)]

"""
    $(TYPEDSIGNATURES)
Compute the mass utilization efficiency at a specific frame of a `Solution`.
"""
function mass_eff(sol::Solution, frame)
    mass_eff = 0.0
    f = sol.frames[frame]
    right_area = f.channel_area[end]
    mdot = sum(prop.flow_rate_kg_s for prop in sol.config.propellants)
    mdot_i = 0.0

    for ions in values(f.ions)
        for ion in ions
            mdot_i += ion.m * ion.nu[end] * right_area
        end
    end

    return mdot_i / mdot
end

"""
    $(TYPEDSIGNATURES)
Compute the mass utilization efficiency at each frame of a `Solution`.
"""
mass_eff(sol::Solution) = [mass_eff(sol, frame) for frame in eachindex(sol.frames)]

const _HC_EV_NM = 1239.8419843320026

_radiative_state_key(state) =
    (Symbol(state.species), Int(state.charge), Int(state.excited_level))

function _radiative_level_energies(chemistry, directories)
    edges = NamedTuple[]
    energies = Dict{Tuple{Symbol, Int, Int}, Float64}()
    for reaction in get(chemistry, "reactions", Any[])
        get(reaction, "type", "") == "electron_impact" || continue
        haskey(reaction, "equation") && haskey(reaction, "rate_coeff_file") || continue
        reactants, products = _parse_reaction_equation(reaction["equation"])
        path = find_file_in_dirs(reaction["rate_coeff_file"], directories)
        isnothing(path) && continue
        threshold = match(r":\s*([-+0-9.eE]+)", open(readline, path))
        isnothing(threshold) && continue
        delta_eV = parse(Float64, threshold.captures[1])
        for upper in keys(products), lower in keys(reactants)
            upper.species == "e" && continue
            lower.species == "e" && continue
            upper.species == lower.species || continue
            upper.charge == lower.charge || continue
            upper.excited_level > lower.excited_level || continue
            push!(edges, (; lower, upper, delta_eV))
            energies[(Symbol(lower.species), Int(lower.charge), 0)] = 0.0
        end
    end

    for _ in 1:(length(edges) + 1)
        changed = false
        for edge in edges
            lower = _radiative_state_key(edge.lower)
            upper = _radiative_state_key(edge.upper)
            if haskey(energies, lower) && !haskey(energies, upper)
                energies[upper] = energies[lower] + edge.delta_eV
                changed = true
            end
        end
        changed || break
    end
    return energies
end

function _radiative_transitions(chemistry, energies)
    transitions = NamedTuple[]
    for reaction in get(chemistry, "reactions", Any[])
        get(reaction, "type", "") == "de-excitation" || continue
        upper = _parse_species_term(reaction["target_species"])
        branches = Int.(reaction["branches"])
        half_lives = Float64.(reaction["half_lives"])
        length(branches) == length(half_lives) ||
            error("`branches` and `half_lives` must have equal length.")
        upper_energy = get(energies, _radiative_state_key(upper), NaN)
        isfinite(upper_energy) || error("Missing energy for $(upper).")

        for (lower_level, half_life) in zip(branches, half_lives)
            lower = RxnTerm(upper.species, upper.charge, lower_level)
            lower_energy = get(energies, _radiative_state_key(lower), NaN)
            isfinite(lower_energy) || error("Missing energy for $(lower).")
            photon_energy_eV = upper_energy - lower_energy
            photon_energy_eV > 0 || continue
            push!(transitions, (;
                upper,
                lower,
                upper_label = string(upper),
                lower_label = string(lower),
                rate_s = log(2.0) / half_life,
                photon_energy_eV,
                wavelength_nm = _HC_EV_NM / photon_energy_eV,
            ))
        end
    end
    return transitions
end

function _radiative_ion_indices(solution)
    indices = Dict{Tuple{Symbol, Int, Int}, Int}()
    for propellant in solution.config.propellants
        index = 0
        for charge in propellant.allowed_charges
            for level in [0; get(propellant.excited_ion_levels, charge, Int[])]
                index += 1
                indices[(propellant.gas.short_name, charge, level)] = index
            end
        end
    end
    return indices
end

function _radiative_state_density(frame, state, ion_indices)
    if state.charge == 0
        key = Symbol(_species_string(state.species, state.charge, state.excited_level))
        haskey(frame.neutrals, key) || error("State $(state) is not saved.")
        return frame.neutrals[key].n
    end
    key = _radiative_state_key(state)
    haskey(ion_indices, key) || error("State $(state) is not saved.")
    return frame.ions[Symbol(state.species)][ion_indices[key]].n
end

function populate_radiative_intensities!(solution)
    chemistry_file = solution.config.propellant_config
    isempty(chemistry_file) && return solution
    chemistry_file = abspath(chemistry_file)
    chemistry = TOML.parsefile(chemistry_file)
    directories = [dirname(chemistry_file); String.(solution.config.reaction_rate_directories)]
    energies = _radiative_level_energies(chemistry, directories)

    species = get(chemistry, "species", Any[])
    if !isempty(species)
        gas = Symbol(species[1]["symbol"])
        path = joinpath(dirname(chemistry_file), "$(gas)_level_energies.toml")
        if isfile(path)
            for (charge, values) in get(TOML.parsefile(path), "level_energies", Dict())
                for (level, energy_eV) in enumerate(values)
                    energies[(gas, parse(Int, charge), level - 1)] = Float64(energy_eV)
                end
            end
        end
    end

    transitions = _radiative_transitions(chemistry, energies)
    isempty(transitions) && return solution
    ion_indices = _radiative_ion_indices(solution)
    for frame in solution.frames
        empty!(frame.intensities)
        for transition in transitions
            intensity = get!(
                () -> zeros(length(solution.grid)),
                frame.intensities,
                transition.wavelength_nm,
            )
            intensity .+= transition.rate_s .* _radiative_state_density(
                frame, transition.upper, ion_indices,
            )
        end
    end
    return solution
end
