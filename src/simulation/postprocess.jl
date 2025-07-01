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
    # TODO: multiple props + containers
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
