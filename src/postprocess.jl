struct HallThrusterSolution{T, U, P, S}
    t::T
    u::U
    savevals::S
    retcode::Symbol
    destats::DiffEqBase.DEStats
    params::P
end

function HallThrusterSolution(sol::S, params::P, savevals::SV) where {S<:SciMLBase.AbstractODESolution, P, SV}
    return HallThrusterSolution(sol.t, sol.u, savevals, sol.retcode, sol.destats, params)
end

"""
    write_restart(path::AbstractString, sol)

Write a restart file to `path``.

This can be reloaded to resume a simulation. The filetype can be anything supported by FileIO, though JLD2 is preferred.
"""
function write_restart(path::AbstractString, sol)
    save(path, Dict(
        "u" =>  sol.u[end],
        "params" => sol.params
    ))
end

"""
    read_restart(path::AbstractString)

Load a restart file from `path`.

The filetype can be anything supported by FileIO, though JLD2 is preferred.
"""
function read_restart(path::AbstractString)
    dict = load(path)
    u, params = dict["u"], dict["params"]
    ncells = length(params.z_cell)-2
    grid = Grid1D(
        ncells,
        params.z_edge,
        params.z_cell,
        params.cell_volume
    )
    B = params.cache.B

    return u, grid, B
end

function Base.show(io::IO, mime::MIME"text/plain", sol::HallThrusterSolution)
    println(io, "Hall thruster solution with $(length(sol.u)) saved frames")
    println(io, "Retcode: $(string(sol.retcode))")
    print(io, "End time: $(sol.t[end]) seconds")
end


function timeaveraged(sol, tstampstart)
    avg = zeros(size(sol.u[1]))
    avg_savevals = deepcopy(sol.savevals[end])
    (;Tev, ue, ϕ, ∇ϕ, ne) = avg_savevals
    Tev .= 0.0
    ue .= 0.0
    ϕ .= 0.0
    ∇ϕ .= 0.0
    ne .= 0.0

    tstamps = length(sol.t)
    Δt = (tstamps - tstampstart + 1)
    for i in tstampstart:length(sol.t)
        avg .+= sol.u[i] / Δt
        Tev .+= sol.savevals[i].Tev / Δt
        ue .+= sol.savevals[i].ue  / Δt
        ϕ .+= sol.savevals[i].ϕ / Δt
        ∇ϕ .+= sol.savevals[i].∇ϕ / Δt
        ne .+= sol.savevals[i].ne / Δt
    end
    return avg, avg_savevals
end

function compute_current(sol)
    index = sol.params.index
    current = zeros(2, length(sol.t))
    area = sol.params.A_ch
    mi = sol.params.propellant.m
    for i in 1:length(sol.t)
        (;ue, ne) = sol.savevals[i]
        current[1, i] = sol.u[i][index.ρiui[1], end-1]*HallThruster.e/mi*area
        current[2, i] = -ne[end-1] * ue[end-1]*HallThruster.e*area
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