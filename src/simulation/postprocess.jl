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
    (;Tev, ue, ϕ, ∇ϕ, ne, ϕ_cell, νan, νc, νei, νen, νiz, νex, νe, Id) = avg_savevals
    Tev .= 0.0
    ue .= 0.0
    ϕ .= 0.0
    ∇ϕ .= 0.0
    ne .= 0.0
    ϕ_cell .= 0.0
    νan .= 0.0
    νc .= 0.0
    νen .= 0.0
    νei .= 0.0

    tstamps = length(sol.t)
    Δt = (tstamps - tstampstart + 1)
    for i in tstampstart:length(sol.t)
        avg .+= sol.u[i] / Δt
        Tev .+= sol.savevals[i].Tev / Δt
        ue .+= sol.savevals[i].ue  / Δt
        ϕ .+= sol.savevals[i].ϕ / Δt
        ϕ_cell .+= sol.savevals[i].ϕ_cell / Δt
        ∇ϕ .+= sol.savevals[i].∇ϕ / Δt
        ne .+= sol.savevals[i].ne / Δt
        νan .+= sol.savevals[i].νan / Δt
        νc .+= sol.savevals[i].νc / Δt
        νen .+= sol.savevals[i].νen / Δt
        νei .+= sol.savevals[i].νei / Δt
        νe .+= sol.savevals[i].νe / Δt
        νiz .+= sol.savevals[i].νiz / Δt
        νex .+= sol.savevals[i].νex / Δt
        Id .+= sol.savevals[i].Id / Δt
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
        current[1, i] = sol.u[i][index.ρiui[1], loc]*HallThruster.e/mi*area
        current[2, i] = -ne[loc] * ue[loc]*HallThruster.e*area
        current[3, i] = current[1, i] + current[2, i]
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

function cut_solution(sol, tstampstart)
    sol_cut = Solution(sol.t[tstampstart:end], sol.u[tstampstart:end], sol.savevals[tstampstart:end], sol.retcode, sol.destats, sol.params)
    return sol_cut
end

function Base.getindex(sol::Solution, field::Symbol, charge::Int = 1)
    mi = sol.params.mi
    index = sol.params.index
    ncells = size(sol.u[1], 2)

    if charge > sol.params.ncharge
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. Maximum charge state in provided solution is $(sol.params.config.ncharge)."))
    end

    if field == :nn
        return [[u[index.ρn, i] / mi for i in 1:ncells] for u in sol.u]
    elseif field == :ni
        return [[u[index.ρi[charge], i] / mi for i in 1:ncells] for u in sol.u]
    elseif field == :ui
        return [[u[index.ρiui[charge], i] / u[index.ρi[charge], i] for i in 1:ncells] for u in sol.u]
    elseif field == :B
        return [sol.params.cache.B]
    elseif field == :ωce
        return [e * sol[:B, charge][1] / me]
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
        ϕ_cell = ϕ_itp.(zs),
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