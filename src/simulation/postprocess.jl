@public Postprocess, time_average, thrust, discharge_current
@public ion_current, electron_current
@public voltage_eff, mass_eff, current_eff, divergence_eff, anode_eff

"""
$(TYPEDEF)
Contains postprocessing options for a given simulation.

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
function time_average(sol::Solution, start_time::AbstractFloat)
    start_frame = findfirst(>=(start_time), sol.t)
    return time_average(sol, start_frame)
end

"""
    $(TYPEDSIGNATURES)
Average a `Solution` over time, starting at frame `start_frame`.
Return a `Solution` object with a single frame containing the averaged simulation properties
"""
function time_average(sol::Solution, start_frame::Integer = 1)
    avg = zeros(size(sol.u[1]))
    avg_savevals = deepcopy(sol.savevals[end])
    fields = fieldnames(typeof(avg_savevals))

    # Initialize avg to zero
    for f in fields
        if f == :anom_variables
            for j in 1:num_anom_variables(sol.config.anom_model)
                avg_savevals[f][j] .= 0.0
            end
        else
            avg_savevals[f] .= 0.0
        end
    end

    # Sum over all timesteps to get average
    tstamps = length(sol.t)
    Δt = (tstamps - start_frame + 1)
    for i in start_frame:length(sol.t)
        avg .+= sol.u[i] / Δt
        for f in fields
            if f == :anom_variables
                for j in 1:num_anom_variables(sol.config.anom_model)
                    avg_savevals[f][j] .+= sol.savevals[i][f][j] / Δt
                end
            else
                avg_savevals[f] .+= sol.savevals[i][f] / Δt
            end
        end
    end

    return Solution(
        sol.t[end:end],
        [avg],
        [avg_savevals],
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
    index = sol.params.index
    left_area = sol.savevals[frame].channel_area[1]
    right_area = sol.savevals[frame].channel_area[end]
    thrust = 0.0
    for Z in 1:(sol.config.ncharge)
        thrust += right_area * sol.u[frame][index.ρiui[Z], end]^2 /
                  sol.u[frame][index.ρi[Z], end]
        thrust -= left_area * sol.u[frame][index.ρiui[Z], 1]^2 /
                  sol.u[frame][index.ρi[Z], 1]
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
thrust(sol::Solution) = [thrust(sol, frame) for frame in eachindex(sol.savevals)]

"""
    $(TYPEDSIGNATURES)
Compute the discharge current at a specific frame of a `Solution`.
"""
discharge_current(sol::Solution, frame::Integer) = sol.savevals[frame].Id[]

"""
    $(TYPEDSIGNATURES)
Compute the discharge current at a each frame of a `Solution`.
"""
discharge_current(sol::Solution) = [s.Id[] for s in sol.savevals]

"""
    $(TYPEDSIGNATURES)
Compute the anode efficiency at a specific frame of a `Solution`.
"""
function anode_eff(sol::Solution, frame::Integer)
    T = thrust(sol, frame)
    current = discharge_current(sol, frame)
    Vd = sol.config.discharge_voltage
    mdot_a = sol.config.anode_mass_flow_rate
    anode_eff = 0.5 * T^2 / current / Vd / mdot_a
    return anode_eff
end

"""
    $(TYPEDSIGNATURES)
Compute the anode efficiency at each frame of a `Solution`.
"""
anode_eff(sol::Solution) = [anode_eff(sol, frame) for frame in eachindex(sol.savevals)]

"""
    $(TYPEDSIGNATURES)
Compute the voltage/acceleration efficiency at a specific frame of a `Solution`.
"""
function voltage_eff(sol::Solution, frame::Integer)
    Vd = sol.config.discharge_voltage
    mi = sol.config.propellant.m

    ui = sol.u[frame][sol.params.index.ρiui[1], end] /
         sol.u[frame][sol.params.index.ρi[1], end]
    voltage_eff = 0.5 * mi * ui^2 / e / Vd
    return voltage_eff
end

"""
    $(TYPEDSIGNATURES)
Compute the voltage/acceleration efficiency at each frame of a `Solution`.
"""
voltage_eff(sol::Solution) = [voltage_eff(sol, frame) for frame in eachindex(sol.savevals)]

"""
    $(TYPEDSIGNATURES)
Compute the divergence efficiency at a specific frame of a `Solution`.
"""
function divergence_eff(sol::Solution, frame::Integer)
    tanδ = sol.savevals[frame].tanδ[end]
    δ = atan(tanδ)
    return cos(δ)^2
end

"""
    $(TYPEDSIGNATURES)
Compute the divergence efficiency at each frame of a `Solution`.
"""
divergence_eff(sol::Solution) = [divergence_eff(sol, frame)
                                 for frame in eachindex(sol.savevals)]

"""
    $(TYPEDSIGNATURES)
Compute the ion current at a specific frame of a `Solution`.
"""
function ion_current(sol::Solution, frame)
    Ii = 0.0
    right_area = sol.savevals[frame].channel_area[end]
    mi = sol.config.propellant.m
    for Z in 1:(sol.config.ncharge)
        Ii += Z * e * sol.u[frame][sol.params.index.ρiui[Z], end] * right_area / mi
    end

    return Ii
end

"""
    $(TYPEDSIGNATURES)
Compute the ion current at each frame of a `Solution`.
"""
ion_current(sol::Solution) = [ion_current(sol, frame) for frame in eachindex(sol.savevals)]

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
    [electron_current(sol, frame) for frame in eachindex(sol.savevals)]
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
current_eff(sol::Solution) = [current_eff(sol, frame) for frame in eachindex(sol.savevals)]

"""
    $(TYPEDSIGNATURES)
Compute the mass utilization efficiency at a specific frame of a `Solution`.
"""
function mass_eff(sol::Solution, frame)
    mass_eff = 0.0
    right_area = sol.savevals[frame].channel_area[end]
    mdot = sol.config.anode_mass_flow_rate

    for Z in 1:(sol.config.ncharge)
        mass_eff += sol.u[frame][sol.params.index.ρiui[Z], end] * right_area / mdot
    end

    return mass_eff
end

"""
    $(TYPEDSIGNATURES)
Compute the mass utilization efficiency at each frame of a `Solution`.
"""
mass_eff(sol::Solution) = [mass_eff(sol, frame) for frame in eachindex(sol.savevals)]

function Base.getindex(sol::Solution, frame::Integer)
    return Solution(
        [sol.t[frame]],
        [sol.u[frame]],
        [sol.savevals[frame]],
        sol.params,
        sol.config,
        sol.retcode,
        sol.error,
    )
end

function Base.getindex(sol::Solution, field::Symbol, charge::Integer = 1)
    if charge > sol.config.ncharge && field in [:ni, :ui, :niui]
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. Maximum charge state in provided solution is $(sol.config.ncharge)."))
    end

    if field == :ni
        return [saved[:ni][charge, :] for saved in sol.savevals]
    elseif field == :ui
        return [saved[:ui][charge, :] for saved in sol.savevals]
    elseif field == :niui
        return [saved[:niui][charge, :] for saved in sol.savevals]
    elseif field == :B
        return [sol.params.cache.B]
    elseif field == :ωce || field == :cyclotron_freq || field == :omega_ce
        return [e * sol[:B, charge][1] / me]
    elseif field == :E
        return -sol[:∇ϕ]
    elseif field == :V || field == :potential
        return sol[:ϕ]
    elseif field == :nu_an
        return sol[:νan]
    elseif field == :nu_ei
        return sol[:νei]
    elseif field == :nu_en
        return sol[:νen]
    elseif field == :mobility
        return sol[:μ]
    elseif field == :thermal_conductivity
        return sol[:κ]
    else
        return [getproperty(saved, field) for saved in sol.savevals]
    end
end
