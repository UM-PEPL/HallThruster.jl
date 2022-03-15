module OVS_Potential

using Symbolics, HallThruster, Plots

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

function verify_potential(ncells, plot = false)

    index = (;ρiui = [1])

    grid = HallThruster.generate_grid(HallThruster.SPT_100, ncells)

    z_cell = grid.cell_centers
    z_edge = grid.edges
    nedges = length(z_edge)

    μ = μ_func.(z_cell)
    ne = ne_func.(z_cell)
    pe = pe_func.(z_cell)
    U = ρiui_func.(z_cell)' |> collect
    ϕ = zeros(nedges)

    ϕ_exact = ϕ_func.(z_edge)

    A = Tridiagonal(ones(nedges-1), ones(nedges), ones(nedges-1))
    b = zeros(nedges) #for potential equation

    ϕ_L = ϕ_exact[1]
    ϕ_R = ϕ_exact[end]

    source_func = (U, params, i) -> source_potential(params.z_edge[i])

    config = (propellant = HallThruster.Xenon, ncharge = 1, source_potential = source_func)
    cache = (;A, b, μ, ϕ, pe, ne)
    params = (;z_cell, z_edge, index, ϕ_L, ϕ_R, cache, config)

    @time HallThruster.solve_potential_edge!(U, params)

    results = (;z_edge, exact = ϕ_exact, sim = params.cache.ϕ)

    if plot
        p = plot(z_edge, results.exact, label = "exact")
        plot!(p, z_edge, results.sim, label = "sim")
        display(p)
    end

    return results
end

end