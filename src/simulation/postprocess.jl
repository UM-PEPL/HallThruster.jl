struct Solution{T, U, P, S, D}
    t::T
    u::U
    savevals::S
    retcode::Symbol
    destats::D
    params::P
end

function Solution(sol::S, params::P, savevals::SV) where {S<:SciMLBase.AbstractODESolution, P, SV}
    return Solution(sol.t, sol.u, savevals, sol.retcode, sol.destats, params)
end

function Base.show(io::IO, mime::MIME"text/plain", sol::Solution)
    println(io, "Hall thruster solution with $(length(sol.u)) saved frames")
    println(io, "Retcode: $(string(sol.retcode))")
    print(io, "End time: $(sol.t[end]) seconds")
end

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
        avg_savevals[f] .= 0.0
    end

    # Sum over all timesteps to get average
    tstamps = length(sol.t)
    Δt = (tstamps - tstampstart + 1)
    for i in tstampstart:length(sol.t)
        avg .+= sol.u[i] / Δt
        for f in fields
            avg_savevals[f] .+= sol.savevals[i][f] / Δt
        end
    end

    return Solution(
        sol.t[end:end],
        [avg],
        [avg_savevals],
        sol.retcode,
        sol.destats,
        sol.params
    )
end

"""
    compute_current(sol, location)
compute current at anode or cathode = outflow in
1D code.
"""
function compute_current(sol, location = "cathode")
    index = sol.params.index
    current = zeros(3, length(sol.t))
    area = sol.params.A_ch
    mi = sol.params.config.propellant.m
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
        Ie = -ne[loc] * ue[loc] * e * area
        current[1, i] = Id - Ie # Ion current
        current[2, i] = Ie      # Electron current
        current[3, i] = Id      # Discharge current
    end
    return current
end

function compute_thrust(sol)
    index = sol.params.index
    thrust = zeros(length(sol.t))
    area = sol.params.A_ch
    for i in 1:length(sol.t)
        for Z in 1:sol.params.config.ncharge
            thrust[i] += area * sol.u[i][index.ρiui[Z], end]^2 / sol.u[i][index.ρi[Z], end]
            thrust[i] -= area * sol.u[i][index.ρiui[Z], 1]^2 / sol.u[i][index.ρi[Z], 1]
        end
    end
    return thrust
end

function compute_anode_eff(sol)
    thrust = compute_thrust(sol)
    current = reduce(vcat, sol[:Id])
    Vd = sol.params.config.discharge_voltage
    mdot_a = sol.params.config.anode_mass_flow_rate
    anode_eff = @. 0.5 * thrust^2 / current / Vd / mdot_a
    return anode_eff
end

function compute_voltage_eff(sol)
    Vd = sol.params.config.discharge_voltage
    nsave = length(sol.t)
    mi = sol.params.config.propellant.m

    ui = sol[:ui, 1]

    voltage_eff = [0.5 * mi * ui[i][end]^2 / e / Vd for i in 1:nsave]
    return voltage_eff
end

function compute_current_eff(sol)
    currents = compute_current(sol)
    Ii = currents[1, :]
    Id = currents[3, :]
    return @. Ii / Id
end

function compute_mass_eff(sol)
    mass_eff = [
        sum(sol.u[i][sol.params.index.ρiui[Z], end] for Z in 1:sol.params.config.ncharge) * sol.params.A_ch /
        sol.params.config.anode_mass_flow_rate
        for i in 1:length(sol.t)
    ]
    return mass_eff
end

function cut_solution(sol, tstampstart)
    sol_cut = Solution(sol.t[tstampstart:end], sol.u[tstampstart:end], sol.savevals[tstampstart:end], sol.retcode, sol.destats, sol.params)
    return sol_cut
end

function Base.getindex(sol::Solution, field::Symbol, charge::Int = 1)
    mi = sol.params.mi
    index = sol.params.index
    ncells = size(sol.u[1], 2)

    if charge > sol.params.ncharge && field in [:ni, :ui, :niui]
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. Maximum charge state in provided solution is $(sol.params.config.ncharge)."))
    end

    if field == :nn
        return [saved[:nn][charge, :] for saved in sol.savevals]
    elseif field == :ni
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

    zs = LinRange(0, 0.05, ncells)

    ui = fill(NaN, length(zs))
    ue = fill(NaN, length(zs))

    cache = (;
        ue = ue,
        Tev = 2/3 * ϵ_itp.(zs),
        pe = ϵ_itp.(zs),
        ne = ne_itp.(zs),
        ∇ϕ = -E_itp.(zs),
        ϕ = ϕ_itp.(zs),
        nn_tot = nn,
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
    )

    retcode = :LANDMARK
    destats = nothing

    u = zeros(4, length(zs))

    ρn = nn_itp.(zs) * mi
    ρi = ne_itp.(zs) * mi
    ρiui = ρi .* ui
    nϵ = ne_itp.(zs) .* ϵ_itp.(zs)

    u[1, :] = ρn
    u[2, :] = ρi
    u[3, :] = ρiui
    u[4, :] = nϵ

    return Solution([0.0], [u], [cache], retcode, destats, params)
end

function frame_dict(sol, frame)
    filler = zeros(length(sol.params.z_cell))
    ncharge = sol.params.config.ncharge
    Dict(
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
