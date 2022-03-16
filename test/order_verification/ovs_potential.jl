module OVS_Potential

using Symbolics, HallThruster, Plots, LinearAlgebra

@variables x t

Dx = Differential(x)

L = 0.05

function sin_wave(var; amplitude, phase, nwaves, offset = 0.0)
    return amplitude * sin(2π * nwaves * var + phase) + offset
end

ϕ = sin_wave(x/L, amplitude = 300, phase = π/2, nwaves = 0.25)
ne = sin_wave(x/L, amplitude = 1e18, phase = π/4, nwaves = 0.5, offset = 1.1e18)
ui = sin_wave(x/L, amplitude = 13000, phase = π/4, nwaves = 0.75, offset = 10000)
μ = sin_wave(x/L, amplitude = 1e3, phase = π/2, nwaves = 1.2, offset = 1.1e4)
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

    grid = HallThruster.generate_grid(HallThruster.SPT_100, ncells)

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

    A = Tridiagonal(ones(ncells-2), ones(ncells-1), ones(ncells-2))
    b = zeros(ncells) #for potential equation

    ϕ_L = ϕ_exact[1]
    ϕ_R = ϕ_exact[end]

    source_func = (U, params, i) -> source_potential(params.z_edge[i])

    config = (propellant = HallThruster.Xenon, ncharge = 1, source_potential = source_func)
    cache = (;A, b, μ, ϕ, pe, ne)
    params = (;z_cell, z_edge, index, ϕ_L, ϕ_R, cache, config)

    HallThruster.solve_potential_edge!(U, params)

    results = (;z = z_cell, exact = ϕ_exact, sim = ϕ)

    if plot_results
        p = plot(z_cell, results.exact, label = "exact")
        plot!(p, z_cell, results.sim, label = "sim")
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

    grid = HallThruster.generate_grid(HallThruster.SPT_100, ncells)

    z_cell = grid.cell_centers
    ncells = length(z_cell)

    μ   = μ_func.(z_cell)
    ne  = ne_func.(z_cell)
    pe  = pe_func.(z_cell)
    ϕ   = ϕ_func.(z_cell)
    ∇ϕ  = zeros(ncells)
    ∇pe = zeros(ncells)
    ue  = zeros(ncells)

    ∇ϕ_exact  = ∇ϕ_func.(z_cell)
    ∇pe_exact = ∇pe_func.(z_cell)
    ue_exact  = ue_func.(z_cell)

    cache = (;μ, ne, pe, ϕ, ∇ϕ, ∇pe, ue)
    params = (;z_cell, cache)
    U = nothing

    HallThruster.compute_gradients!(∇ϕ, ∇pe, ue, U, params)

    result_∇ϕ = (;z = z_cell, exact = ∇ϕ_exact, sim = ∇ϕ)
    result_∇pe = (;z = z_cell, exact = ∇pe_exact, sim = ∇pe)
    result_ue = (;z = z_cell, exact = ue_exact, sim = ue)

    if plot_results
        p1 = plot(result_∇ϕ.z, result_∇ϕ.exact, label = "exact",   title = "Potential gradient")
        plot!(p1, result_∇ϕ.z, result_∇ϕ.sim, label = "sim")
        p2 = plot(result_∇pe.z, result_∇pe.exact, label = "", title = "Pressure gradient")
        plot!(p2, result_∇pe.z, result_∇pe.sim, label = "")
        p3 = plot(result_ue.z, result_ue.exact, label = "",   title = "Electron velocity")
        plot!(p3, result_ue.z, result_ue.sim, label = "")

        plot(p1, p2, p3, layout = (3, 1)) |> display
    end

    return (result_∇ϕ, result_∇pe, result_ue)
end


end