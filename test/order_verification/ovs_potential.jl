module OVS_Potential

using Symbolics, HallThruster, Plots, LinearAlgebra

include("ovs_funcs.jl")

@variables x t

Dx = Differential(x)

L = 0.05

ϕ = sin_wave(x/L, amplitude = 250, phase = π/2, nwaves = 0.75, offset = 500.0)
ne = sin_wave(x/L, amplitude = 1e18, phase = π/4, nwaves = 0.5, offset = 1.1e18)
ui = sin_wave(x/L, amplitude = 13000, phase = π/4, nwaves = 0.83, offset = 10000)
μ = sin_wave(x/L, amplitude = 1e4, phase = π/2, nwaves = 4, offset = 1.1e4)
ϵ = sin_wave(x/L, amplitude = 40, phase = 3π/4, nwaves = 1.0, offset = 43)
ρiui = ne * ui * HallThruster.Xenon.m
pe = ne * ϵ

ϕ_func = eval(build_function(ϕ, [x]))
ne_func = eval(build_function(ne, [x]))
μ_func = eval(build_function(μ, [x]))
ρiui_func = eval(build_function(ρiui, [x]))
pe_func = eval(build_function(pe, [x]))

potential_eq = Dx(μ * ne * Dx(ϕ) - ne * ui - μ * Dx(pe))
source_potential = eval(build_function(expand_derivatives(potential_eq), [x]))

function verify_potential(ncells, plot_results = false)

    index = (;ρiui = [1])

    grid = HallThruster.generate_grid(HallThruster.SPT_100.geometry, ncells, (0.0, 0.05))

    z_cell = grid.cell_centers
    z_edge = grid.edges
    ncells = length(z_cell)
    nedges = length(z_edge)

    μ = μ_func.(z_cell)
    ne = ne_func.(z_cell)
    pe = pe_func.(z_cell)
    U = ρiui_func.(z_cell)' |> collect
    ϕ = zeros(nedges)

    ϕ_exact = ϕ_func.(z_edge)

    A = Tridiagonal(ones(nedges-1), ones(nedges), ones(nedges-1))
    b = zeros(nedges)

    ϕ_L = ϕ_exact[1]
    ϕ_R = ϕ_exact[end]

    source_func = (U, params, i) -> source_potential(params.z_edge[i])

    config = (propellant = HallThruster.Xenon, ncharge = 1, source_potential = source_func, LANDMARK = true)
    cache = (;A, b, μ, ϕ, pe, ne)
    params = (;z_cell, z_edge, index, ϕ_L, ϕ_R, cache, config)

    HallThruster.solve_potential_edge!(U, params)

    results = (;z = z_edge, exact = ϕ_exact, sim = ϕ)

    if plot_results
        p = plot(results.z, results.exact, label = "exact")
        plot!(p, results.z, results.sim, label = "sim")
        display(p)
    end

    return (results,)
end

∇pe = Dx(pe)
∇ϕ = Dx(ϕ)
ue = μ * (∇ϕ - ∇pe / ne)

∇pe_func = eval(build_function(expand_derivatives(∇pe), [x]))
∇ϕ_func = eval(build_function(expand_derivatives(∇ϕ), [x]))
ue_func = eval(build_function(expand_derivatives(ue), [x]))

function verify_gradients(ncells, plot_results = false)

    grid = HallThruster.generate_grid(HallThruster.SPT_100.geometry, ncells, (0.0, 0.05))

    z_cell, z_edge = grid.cell_centers, grid.edges
    ncells = length(z_cell)

    μ   = μ_func.(z_cell)
    ne  = ne_func.(z_cell)
    pe  = pe_func.(z_cell)
    ϕ   = ϕ_func.(z_edge)
    ϕ_cell = zeros(ncells)
    ∇ϕ  = zeros(ncells)
    ∇pe = zeros(ncells)

    ∇ϕ_exact  = ∇ϕ_func.(z_cell)
    ∇pe_exact = ∇pe_func.(z_cell)

    cache = (;μ, ne, pe, ϕ, ϕ_cell, ∇ϕ, ∇pe)
    params = (;z_cell, cache, z_edge, config = (;LANDMARK = true))

    HallThruster.compute_gradients!(∇ϕ, ∇pe, params)

    result_∇ϕ = (;z = z_cell, exact = ∇ϕ_exact, sim = ∇ϕ)
    result_∇pe = (;z = z_cell, exact = ∇pe_exact, sim = ∇pe)

    if plot_results
        p1 = plot(result_∇ϕ.z, result_∇ϕ.exact, label = "exact",   title = "Potential gradient")
        plot!(p1, result_∇ϕ.z, result_∇ϕ.sim, label = "sim")
        p2 = plot(result_∇pe.z, result_∇pe.exact, label = "", title = "Pressure gradient")
        plot!(p2, result_∇pe.z, result_∇pe.sim, label = "")
        p3 = plot(result_ue.z, result_ue.exact, label = "",   title = "Electron velocity")
        plot!(p3, result_ue.z, result_ue.sim, label = "")

        plot(p1, p2, p3, layout = (3, 1)) |> display
    end

    return (result_∇ϕ, result_∇pe)
end


end