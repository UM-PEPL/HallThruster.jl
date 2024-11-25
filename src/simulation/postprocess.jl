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
            for j in 1:num_anom_variables(sol.params.config.anom_model)
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
                for j in 1:num_anom_variables(sol.params.config.anom_model)
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
    for Z in 1:(sol.params.config.ncharge)
        thrust += right_area * sol.u[frame][index.ρiui[Z], end]^2 /
                  sol.u[frame][index.ρi[Z], end]
        thrust -= left_area * sol.u[frame][index.ρiui[Z], 1]^2 /
                  sol.u[frame][index.ρi[Z], 1]
    end

    # Multiply by sqrt of divergence efficiency to model loss of ions in radial direction
    if (sol.params.config.apply_thrust_divergence_correction)
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
    Vd = sol.params.config.discharge_voltage
    mdot_a = sol.params.config.anode_mass_flow_rate
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
    Vd = sol.params.config.discharge_voltage
    mi = sol.params.config.propellant.m

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
    mi = sol.params.config.propellant.m
    for Z in 1:(sol.params.config.ncharge)
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
    mdot = sol.params.config.anode_mass_flow_rate

    for Z in 1:(sol.params.config.ncharge)
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
        sol.retcode,
        sol.error,
    )
end

function Base.getindex(sol::Solution, field::Symbol, charge::Integer = 1)
    if charge > sol.params.config.ncharge && field in [:ni, :ui, :niui]
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. Maximum charge state in provided solution is $(sol.params.config.ncharge)."))
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

function load_landmark_data(case; ncells = 100)
    fluid_1 = load_landmark_data(case, "fluid_1"; ncells)
    fluid_2 = load_landmark_data(case, "fluid_2"; ncells)
    hybrid = load_landmark_data(case, "hybrid"; ncells)
    return (fluid_1, fluid_2, hybrid)
end

function load_landmark_data(case, suffix; ncells = 100)
    folder = joinpath(PACKAGE_ROOT, "landmark", "case_$case")

    E_file = readdlm(joinpath(folder, "electric_field_$(suffix).csv"), ',')
    z_E, E = E_file[:, 1], E_file[:, 2]
    E_itp = LinearInterpolation(z_E, E)

    energy_file = readdlm(joinpath(folder, "energy_$(suffix).csv"), ',')
    z_ϵ, ϵ = energy_file[:, 1], energy_file[:, 2]
    ϵ_itp = LinearInterpolation(z_ϵ, ϵ)

    neutral_density_file = readdlm(joinpath(folder, "neutral_density_$(suffix).csv"), ',')
    z_nn, nn = neutral_density_file[:, 1], neutral_density_file[:, 2]
    nn_itp = LinearInterpolation(z_nn, nn)

    plasma_density_file = readdlm(joinpath(folder, "plasma_density_$(suffix).csv"), ',')
    z_ne, ne = plasma_density_file[:, 1], plasma_density_file[:, 2]
    ne_itp = LinearInterpolation(z_ne, ne)

    potential_file = readdlm(joinpath(folder, "potential_$(suffix).csv"), ',')
    z_ϕ, ϕ = potential_file[:, 1], potential_file[:, 2]
    ϕ_itp = LinearInterpolation(z_ϕ, ϕ)

    zs = range(0, 0.05, length = ncells)

    ui = fill(NaN, length(zs))
    ue = fill(NaN, length(zs))

    nes = ne_itp.(zs)
    nns = nn_itp.(zs)

    cache = (;
        ue = ue,
        Tev = 2 / 3 * ϵ_itp.(zs),
        pe = ϵ_itp.(zs),
        ne = ne_itp.(zs),
        ni = nes' |> collect,
        ui = ui' |> collect,
        niui = ui' |> collect,
        ∇ϕ = -E_itp.(zs),
        ϕ = ϕ_itp.(zs),
        nn = nns,
    )

    ionization_reactions = HallThruster.load_ionization_reactions(
        :Landmark, [Xenon(0), Xenon(1)],
    )

    mi = Xenon.m

    params = (;
        ncharge = 1,
        z_cell = zs,
        z_edge = zs,
        L_ch = 0.025,
        cache,
        mi,
        index = (;
            ρn = 1,
            ρi = [2],
            ρiui = [3],
            nϵ = 4,
        ),
        ionization_reactions,
        config = (;
            propellant = Xenon,
            ionization_model = :Landmark,
        ),
    )

    retcode = :LANDMARK

    u = zeros(4, length(zs))

    ρn = nn_itp.(zs) * mi
    ρi = ne_itp.(zs) * mi
    ρiui = ρi .* ui
    nϵ = ne_itp.(zs) .* ϵ_itp.(zs)

    u[1, :] = ρn
    u[2, :] = ρi
    u[3, :] = ρiui
    u[4, :] = nϵ

    return Solution([0.0], [u], [cache], params, retcode, "")
end
