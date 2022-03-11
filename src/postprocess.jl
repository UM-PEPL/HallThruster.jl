using JLD2

struct HallThrusterSolution{T, U, P}
    t::T
    u::U
    retcode::Symbol
    destats::DiffEqBase.DEStats
    params::P
end

function HallThrusterSolution(sol::S, params::P) where {S<:SciMLBase.AbstractODESolution, P}
    return HallThrusterSolution(sol.t, sol.u, sol.retcode, sol.destats, params)
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

function Base.getindex(sol::HallThrusterSolution, I1, I2)
    return sol.u[I1][I2, :]
end

function Base.getindex(sol::HallThrusterSolution, I1, field::String)
    index = sol.params.index
    mi = sol.params.fluids[1].species.element.m
    params = sol.params
    if field == "nn"
        return sol[I1, 1] / mi
    elseif field == "ni"
        return sol[I1, 2] / mi
    elseif field == "niui"
        return sol[I1, 3] / mi
    elseif field == "ui"
        return sol[I1, 3] ./ sol[I1, 2]
    elseif field == "ne"
        return sol[I1, index.ne]
    elseif field == "Te"
        return sol[I1, index.Tev]
    elseif field == "nϵ"
        return sol[I1, index.nϵ]
    elseif field == "pe"
        return sol[I1, index.pe]
    elseif field == "E"
        return -sol[I1, index.grad_ϕ]
    elseif field == "ϕ"
        return sol[I1, index.ϕ]
    elseif field == "B"
        return params.cache.B
    elseif field == "νan"
        return [get_v_an(z, B, params.L_ch) for (z, B) in zip(params.z_cell, params.cache.B)]
    elseif field == "νc"
        return @views [
                electron_collision_freq(Te, nn, ne, mi)
                for (Te, nn, ne) in zip(
                    sol.u[I1][index.Tev, :],
                    sol.u[I1][1, :]/mi,
                    sol.u[I1][2, :]/mi,
                )]
    elseif field == "z"
        return params.z_cell
    else
        throw(ArgumentError("Hall thruster has no field \"$(field)\""))
    end
end

Base.firstindex(s::HallThrusterSolution, args...) = Base.firstindex(s.u, args...)
Base.lastindex(s::HallThrusterSolution, args...) = Base.lastindex(s.u, args...)

#=
function extract_data(u::Matrix{T}, config)
    nvars, ncells = u
    d = Dict{String, Vector{T}}
    mi = config.propellant.m
    d["nn"] = u[1, :] / mi

    offset = 1
    if config.solve_neutral_velocity
        d["un"] = u[2, :] / mi
        offset += 1
    else
        d["un"] = config.un * ones(T,ncells)
    end

    for i in 1:config.ncharge
        d["ni($(i))"] = u[offset, :]
        d["ui($(i))"] = u[offset+1, :] ./ u[offset, :]
        if config.solve_ion_temperature
            d["Ti($(i))"] = u[offset+2, :] ./ u[offset, :]
            offset += 3
        else
            d["Ti($i))"] = config.ion_temperature * ones(T, ncells)
            offset += 2
        end
    end

    index = (;lf = lf, nϵ = lf+1, Tev = lf+2, ne = lf+3, pe = lf+4, ϕ = lf+5, grad_ϕ = lf+6, ue = lf+7)
    d["nϵ"] = u[offset, :]
    d["Te"] = u[offset+1, :]
    d["ne"] = u[offset+2, :]
    d["ϕ"] = u[offset+4, :]
    d["E"] = -u[offset+5, :]
    d["ue"] = u[offset+6, :]

    return d
end


Base.@kwdef struct HallThrusterConfig
    geometry
    ncells::Int
    restart::Bool
    restart_filename::String
    bfield
    propellant
    ncharge
    solve_neutral_velocity
    solve_ion_temperature
    anom_model
    anom_model_coeffs
    energy_equation
    ion_diffusion_coeff
    electron_pressure_coupled
    flux
    reconstruct
    limiter
    mass_flow_rate
    cathode_potential
    anode_potential
    cathode_Te
    anode_Te
end
=#