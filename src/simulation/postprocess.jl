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

"""
    write_restart(path::AbstractString, sol)

Write a restart file to `path``.

This can be reloaded to resume a simulation. The filetype can be anything supported by FileIO, though JLD2 is preferred.
"""
function write_restart(path::AbstractString, sol)
    save(path, Dict(
        "t" => sol.t,
        "u" =>  sol.u,
        "savevals" => sol.savevals,
        "z_edge" => sol.params.z_edge,
        "z_cell" => sol.params.z_cell,
        "A_ch" => sol.params.A_ch,
        "B" => sol.params.cache.B,
        "index" => sol.params.index,
    ))
end

"""
    read_restart(path::AbstractString)

Load a restart file from `path`.

The filetype can be anything supported by FileIO, though JLD2 is preferred.
"""
function read_restart(path::AbstractString)
    dict = load(path)

    params = (;
        cache = (;B = dict["B"]),
        A_ch = dict["A_ch"],
        z_edge = dict["z_edge"],
        z_cell = dict["z_cell"],
        index = dict["index"]
    )

    return Solution(
        dict["t"], dict["u"], dict["savevals"], :Restart, nothing, params
    )
end

function Base.show(io::IO, mime::MIME"text/plain", sol::Solution)
    println(io, "Hall thruster solution with $(length(sol.u)) saved frames")
    println(io, "Retcode: $(string(sol.retcode))")
    print(io, "End time: $(sol.t[end]) seconds")
end

"""
    time_average(sol, tstampstart)
compute timeaveraged solution, input Solution type and the frame at which averaging starts. 
Returns a Solution object with a single frame.
"""
function time_average(sol, tstampstart = 1)
    avg = zeros(size(sol.u[1]))
    avg_savevals = deepcopy(sol.savevals[end])
    (;Tev, ue, ϕ, ∇ϕ, ne, ϕ_cell, νan, νc, νei, νen) = avg_savevals
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
        end
    end
    return thrust
end

function cut_solution(sol, tstampstart)
    sol_cut = Solution(sol.t[tstampstart:end], sol.u[tstampstart:end], sol.savevals[tstampstart:end], sol.retcode, sol.destats, sol.params)
    return sol_cut
end

function Base.getindex(sol::Solution, field::Symbol, charge::Int = 1)
    mi = sol.params.config.propellant.m
    index = sol.params.index
    ncells = size(sol.u[1], 2)

    if charge > sol.params.config.ncharge
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