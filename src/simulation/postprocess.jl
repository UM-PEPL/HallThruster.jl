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
        sol.t[end:end], [avg], [avg_savevals], sol.params, sol.retcode, sol.error
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
        (; ue, ne) = sol.savevals[i]
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

    ui = sol.u[frame][sol.params.index.ρiui[1], end] /
         sol.u[frame][sol.params.index.ρi[1], end]
    voltage_eff = 0.5 * mi * ui^2 / e / Vd
    return voltage_eff
end

function divergence_eff(sol, frame)
    tanδ = sol.savevals[frame].tan_div_angle[end]
    δ = atan(tanδ)
    return cos(δ)^2
end

function ion_current(sol, frame)
    Ii = 0.0
    right_area = sol.savevals[frame].channel_area[end]
    mi = sol.params.config.propellant.m
    for Z in 1:(sol.params.config.ncharge)
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

    for Z in 1:(sol.params.config.ncharge)
        mass_eff += sol.u[frame][sol.params.index.ρiui[Z], end] * right_area / mdot
    end

    return mass_eff
end

for func in [
    :thrust,
    :discharge_current,
    :ion_current,
    :electron_current,
    :mass_eff,
    :voltage_eff,
    :anode_eff,
    :current_eff,
    :divergence_eff
]
    eval(
        quote
        $(func)(sol) = [$(func)(sol, i) for i in eachindex(sol.t)]
    end,
    )
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

    zs = range(0, 0.05; length = ncells)

    ui = fill(NaN, length(zs))
    ue = fill(NaN, length(zs))

    nes = ne_itp.(zs)
    nns = nn_itp.(zs)
    nϵ = ne_itp.(zs) .* ϵ_itp.(zs)

    cache = (;
        ue = ue,
        Tev = 2 / 3 * ϵ_itp.(zs),
        pe = ϵ_itp.(zs),
        ne = ne_itp.(zs),
        ni = collect(nes'),
        ui = collect(ui'),
        niui = collect(ui'),
        electric_field = E_itp.(zs),
        potential = ϕ_itp.(zs),
        nn = nns,
        nϵ
    )

    ionization_reactions = HallThruster.load_ionization_reactions(
        :Lookup, [Xenon(0), Xenon(1)]
    )

    mi = Xenon.m

    params = (;
        ncharge = 1,
        z_cell = zs,
        z_edge = zs,
        L_ch = 0.025,
        cache,
        mi,
        index = (; ρn = 1, ρi = [2], ρiui = [3]),
        ionization_reactions,
        config = (; propellant = Xenon, ionization_model = :Lookup)
    )

    retcode = :LANDMARK

    u = zeros(4, length(zs))

    ρn = nn_itp.(zs) * mi
    ρi = ne_itp.(zs) * mi
    ρiui = ρi .* ui

    u[1, :] = ρn
    u[2, :] = ρi
    u[3, :] = ρiui

    return Solution([0.0], [u], [cache], params, retcode, "")
end
