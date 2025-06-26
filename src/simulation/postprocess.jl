@public Postprocess, time_average, thrust, discharge_current
@public ion_current, electron_current
@public voltage_eff, mass_eff, current_eff, divergence_eff, anode_eff

"""
$(TYPEDEF)
Contains postprocessing options for a given simulation.
When `run_simulation(config, sim_params; postprocess = Postprocess(...))` is called with a non-empty `output_file`, `HallThruster` will write the simulation results to a JSON file.
The results in the file will be transformed according to the fields.

## Fields
$(TYPEDFIELDS)
"""
@kwdef struct Postprocess
    """
    The file to which the output will be written. If empty, no output will be written.
    """
    output_file::String = ""
    """
    The time to begin averaging at. If less than zero, no averaged output will be written.
    """
    average_start_time::Float64 = -1
    """
    Whether time-resolved output will be saved.
    If true, each frame of the simulation will be written to the output file.
    """
    save_time_resolved::Bool = false
end

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
    avg_frames = deepcopy(sol.frames[end])
    fields = fieldnames(typeof(avg_frames))

    # Initialize avg to zero
    for f in fields
        if f == :anom_variables
            for j in 1:num_anom_variables(sol.config.anom_model)
                avg_frames[f][j] .= 0.0
            end
        else
            avg_frames[f] .= 0.0
        end
    end

    # Sum over all timesteps to get average
    tstamps = length(sol.t)
    Δt = (tstamps - start_frame + 1)
    for i in start_frame:length(sol.t)
        for f in fields
            if f == :anom_variables
                for j in 1:num_anom_variables(sol.config.anom_model)
                    avg_frames[f][j] .+= sol.frames[i][f][j] / Δt
                end
            else
                avg_frames[f] .+= sol.frames[i][f] / Δt
            end
        end
    end

    return Solution(
        sol.t[end:end],
        [avg_frames],
        sol.params,
        sol.config,
        sol.retcode,
        sol.error,
    )
end

"""
    $(TYPEDSIGNATURES)
Compute the thrust at a specific frame of a `Solution`.
"""
function thrust(sol::Solution, frame::Integer)
    mi = sol.params.propellants[1].gas.m
    f = sol.frames[frame]
    left_area = f.channel_area[begin]
    right_area = f.channel_area[end]
    thrust = 0.0
    #TODO: mutliple props
    for Z in 1:sol.config.propellants[1].max_charge
        thrust += right_area * mi * f.niui[Z, end]^2 / f.ni[Z, end]
        thrust -= left_area * mi * f.niui[Z, begin]^2 / f.ni[Z, begin]
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
discharge_current(sol::Solution, frame::Integer) = sol.frames[frame].Id[]

"""
    $(TYPEDSIGNATURES)
Compute the discharge current at a each frame of a `Solution`.
"""
discharge_current(sol::Solution) = [s.Id[] for s in sol.frames]

"""
    $(TYPEDSIGNATURES)
Compute the anode efficiency at a specific frame of a `Solution`.
"""
function anode_eff(sol::Solution, frame::Integer)
    T = thrust(sol, frame)
    current = discharge_current(sol, frame)
    Vd = sol.config.discharge_voltage
    mdot_a = sol.config.propellants[1].flow_rate_kg_s
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
    mi = sol.config.propellants[1].gas.m
    ui = sol.frames[frame].niui[1, end] / sol.frames[frame].ni[1, end]
    voltage_eff = 0.5 * mi * ui^2 / e / Vd
    return voltage_eff
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
    tanδ = sol.frames[frame].tanδ[end]
    δ = atan(tanδ)
    return cos(δ)^2
end

"""
    $(TYPEDSIGNATURES)
Compute the divergence efficiency at each frame of a `Solution`.
"""
divergence_eff(sol::Solution) = [
    divergence_eff(sol, frame)
        for frame in eachindex(sol.frames)
]

"""
    $(TYPEDSIGNATURES)
Compute the ion current at a specific frame of a `Solution`.
"""
function ion_current(sol::Solution, frame)
    Ii = 0.0
    f = sol.frames[frame]
    for Z in 1:(sol.config.propellants[1].max_charge)
        Ii += Z * e * f.niui[Z, end] * f.channel_area[end]
    end

    return Ii
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
electron_current(sol::Solution, frame) = discharge_current(sol, frame) -
    ion_current(sol, frame)
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
    prop = sol.config.propellants[1]
    mi = prop.gas.m
    mdot = prop.flow_rate_kg_s

    for Z in 1:prop.max_charge
        mass_eff += mi * f.niui[Z, end] * right_area / mdot
    end

    return mass_eff
end

"""
    $(TYPEDSIGNATURES)
Compute the mass utilization efficiency at each frame of a `Solution`.
"""
mass_eff(sol::Solution) = [mass_eff(sol, frame) for frame in eachindex(sol.frames)]
