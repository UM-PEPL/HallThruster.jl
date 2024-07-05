"""
    time_average(sol, tstampstart)
compute time-averaged solution, input Solution type and the frame at which averaging starts.
Returns a Solution object with a single frame.
"""
function time_average(sol::Solution, tstampstart = 1)
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
    Δt = (tstamps - tstampstart + 1)
    for i in tstampstart:length(sol.t)
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
        sol.retcode,
        sol.params
    )
end

"""
    compute_current(sol, location)
compute current at anode or cathode = outflow in
1D code.
"""
function compute_current(sol, location = "cathode")
    current = zeros(3, length(sol.t))
    if location == "cathode"
        loc = length(sol.savevals[1].ue) - 1
    elseif location == "anode"
        loc = 2
    else
        error("Type anode or cathode as location argument")
    end
    for i in 1:length(sol.t)
        (;ue, ne) = sol.savevals[i]
        Id = sol.savevals[i].Id[]
        Ie = -ne[loc] * ue[loc] * e * sol.savevals[i].channel_area[loc]
        current[1, i] = Id - Ie # Ion current
        current[2, i] = Ie      # Electron current
        current[3, i] = Id      # Discharge current
    end
    return current
end

function thrust(sol, frame)
    index = sol.params.index
    left_area = sol.savevals[frame].channel_area[1]
    right_area = sol.savevals[frame].channel_area[end]
    thrust = 0.0
    for Z in 1:sol.params.config.ncharge
        thrust += right_area * sol.u[frame][index.ρiui[Z], end]^2 / sol.u[frame][index.ρi[Z], end]
        thrust -= left_area * sol.u[frame][index.ρiui[Z], 1]^2 / sol.u[frame][index.ρi[Z], 1]
    end

    # Multiply by sqrt of divergence efficiency to model loss of ions in radial direction
    if (sol.params.config.apply_thrust_divergence_correction)
        return thrust * sqrt(divergence_eff(sol, frame))
    else
        return thrust
    end
end

discharge_current(sol, frame) = sol.savevals[frame].Id[]

function anode_eff(sol, frame)
    T = thrust(sol, frame)
    current = discharge_current(sol, frame)
    Vd = sol.params.config.discharge_voltage
    mdot_a = sol.params.config.anode_mass_flow_rate
    anode_eff = 0.5 * T^2 / current / Vd / mdot_a
    return anode_eff
end

function voltage_eff(sol, frame)
    Vd = sol.params.config.discharge_voltage
    mi = sol.params.config.propellant.m

    ui = sol.u[frame][sol.params.index.ρiui[1], end] / sol.u[frame][sol.params.index.ρi[1], end]
    voltage_eff = 0.5 * mi * ui^2 / e / Vd
    return voltage_eff
end

function divergence_eff(sol, frame)
    tanδ = sol.savevals[frame].tanδ[end]
    δ = atan(tanδ)
    return cos(δ)^2
end

function ion_current(sol, frame)
    Ii = 0.0
    right_area = sol.savevals[frame].channel_area[end]
    mi = sol.params.config.propellant.m
    for Z in 1:sol.params.config.ncharge
        Ii += Z * e * sol.u[frame][sol.params.index.ρiui[Z], end] * right_area / mi
    end

    return Ii
end

electron_current(sol, frame) = discharge_current(sol, frame) - ion_current(sol, frame)
current_eff(sol, frame) = ion_current(sol, frame) / discharge_current(sol, frame)

function mass_eff(sol, frame)
    mass_eff = 0.0
    right_area = sol.savevals[frame].channel_area[end]
    mdot = sol.params.config.anode_mass_flow_rate

    for Z in 1:sol.params.config.ncharge
        mass_eff += sol.u[frame][sol.params.index.ρiui[Z], end] * right_area / mdot
    end

    return mass_eff
end

for func in [:thrust, :discharge_current, :ion_current, :electron_current, :mass_eff, :voltage_eff, :anode_eff, :current_eff, :divergence_eff]
    eval(quote
        $(func)(sol) = [$(func)(sol, i) for i in eachindex(sol.t)]
    end)
end

function cut_solution(sol, tstampstart)
    sol_cut = Solution(sol.t[tstampstart:end], sol.u[tstampstart:end], sol.savevals[tstampstart:end], sol.retcode, sol.params)
    return sol_cut
end

function Base.getindex(sol::Solution, field::Symbol, charge::Int = 1)

    if charge > sol.params.ncharge && field in [:ni, :ui, :niui]
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
    elseif field == :ωce
        return [e * sol[:B, charge][1] / me]
    elseif field == :E
        return -sol[:∇ϕ]
    else
        return [getproperty(saved, field) for saved in sol.savevals]
    end
end

function load_landmark_data(case; ncells = 100)
    fluid_1 = load_landmark_data(case, "fluid_1"; ncells)
    fluid_2 = load_landmark_data(case, "fluid_2"; ncells)
    hybrid  = load_landmark_data(case, "hybrid"; ncells)
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
        Tev = 2/3 * ϵ_itp.(zs),
        pe = ϵ_itp.(zs),
        ne = ne_itp.(zs),
        ni = nes' |> collect,
        ui = ui' |> collect,
        niui = ui' |> collect,
        ∇ϕ = -E_itp.(zs),
        ϕ = ϕ_itp.(zs),
        nn = nns,
    )

    ionization_reactions = HallThruster.load_reactions(LandmarkIonizationLookup(), [Xenon(0), Xenon(1)]);

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
        )
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

    return Solution([0.0], [u], [cache], retcode, params)
end

function frame_dict(sol, frame)
    filler = zeros(length(sol.params.z_cell))
    ncharge = sol.params.config.ncharge
    Dict(
        "thrust" => thrust(sol, frame),
        "discharge_current" => discharge_current(sol, frame),
        "mass_eff" => mass_eff(sol, frame),
        "voltage_eff" => voltage_eff(sol, frame),
        "current_eff" => current_eff(sol, frame),
        "t" => sol.t[frame],
        "z" => sol.params.z_cell,
        "nn" => sol[:nn, 1][frame],
        "ni_1" => sol[:ni, 1][frame],
        "ni_2" => (ncharge > 1 ? sol[:ni, 2][frame] : filler),
        "ni_3" => (ncharge > 2 ? sol[:ni, 3][frame] : filler),
        "ne" => sol[:ne][frame],
        "ui_1" => sol[:ui, 1][frame],
        "ui_2" => (ncharge > 1 ? sol[:ui, 2][frame] : filler),
        "ui_3" => (ncharge > 2 ? sol[:ui, 3][frame] : filler),
        "niui_1" => sol[:niui, 1][frame],
        "niui_2" => (ncharge > 1 ? sol[:niui, 2][frame] : filler),
        "niui_3" => (ncharge > 2 ? sol[:niui, 3][frame] : filler),
        "ue" => sol[:ue][frame],
        "V" => sol[:ϕ][frame],
        "E" => sol[:E][frame],
        "Tev" => sol[:Tev][frame],
        "pe" => sol[:pe][frame],
        "grad_pe" => sol[:∇pe][frame],
        "nu_en" => sol[:νen][frame],
        "nu_ei" => sol[:νei][frame],
        "nu_anom" => sol[:νan][frame],
        "nu_class" => sol[:νc][frame],
        "mobility" => sol[:μ][frame],
    )
end

function write_to_json(filename, sol)
    num_frames = length(sol.t)
    frames = map(frame -> frame_dict(sol, frame), 1:num_frames)
    JSON3.write(filename, frames)
end
