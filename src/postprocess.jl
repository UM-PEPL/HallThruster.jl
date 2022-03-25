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
    return avg, avg_savevals
end

function compute_current(sol)
    index = sol.params.index
    current = zeros(3, length(sol.t))
    area = sol.params.A_ch
    mi = sol.params.propellant.m
    for i in 1:length(sol.t)
        (;ue, ne) = sol.savevals[i]
        current[1, i] = sol.u[i][index.ρiui[1], end-1]*HallThruster.e/mi*area
        current[2, i] = -ne[end-1] * ue[end-1]*HallThruster.e*area
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

using DelimitedFiles

function load_hallis_output(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei,
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    return output[1:end-1, :]
end


function load_hallis_for_input()
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_01.out")
    ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, hallis.ϕ)
    grad_ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, -hallis.Ez)
    return ϕ_hallis, grad_ϕ_hallis
end
