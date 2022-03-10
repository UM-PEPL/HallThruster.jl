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