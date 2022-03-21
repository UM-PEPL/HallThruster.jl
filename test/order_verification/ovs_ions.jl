module OVS_Ions

include("ovs_funcs.jl")

using Symbolics, HallThruster, Plots, LinearAlgebra, OrdinaryDiffEq

@variables x t

Dt = Differential(t)
Dx = Differential(x)

const k_ionization = HallThruster.ionization_fits_Xe(1)[1].rate_coeff
const un = 150
const mi = HallThruster.Xenon.m
const e = HallThruster.e
const Ti = 1000
const L = 0.025

ϕ = sin_wave(x/L, amplitude = 300, phase = π/2, nwaves = 0.25)
ne = sin_wave(x/L, amplitude = 1e16, phase = π/4, nwaves = 0.5, offset = 1.1e16)
nn = sin_wave(x/L, amplitude = 5e18, phase = pi/3, nwaves = 2.0, offset = 6e18)
ui = sin_wave(x/L, amplitude = 13000, phase = π/4, nwaves = 0.75, offset = 10000)
μ = sin_wave(x/L, amplitude = 1e4, phase = π/2, nwaves = 1.2, offset = 1.1e4)
ϵ = sin_wave(x/L, amplitude = 20, phase = 1.3*π/2, nwaves = 1.1, offset = 25)
nϵ = ne * ϵ
∇ϕ = Dx(ϕ)
∇pe = Dx(nϵ)
ρiui = ne * ui * HallThruster.Xenon.m
ρn = nn * HallThruster.Xenon.m
ue = μ * (∇ϕ - ∇pe/ne)
ρi = ne * HallThruster.Xenon.m
p = ne * HallThruster.kB * Ti
∇p = Dx(p)

ϕ_func = eval(build_function(ϕ, [x]))
ne_func = eval(build_function(ne, [x]))
μ_func = eval(build_function(μ, [x]))
ρiui_func = eval(build_function(ρiui, [x]))
nϵ_func = eval(build_function(nϵ, [x]))
ue_func = eval(build_function(expand_derivatives(ue), [x]))
∇ϕ_func = eval(build_function(expand_derivatives(∇ϕ), [x]))
ρn_func = eval(build_function(ρn, [x]))
ρi_func = eval(build_function(ρi, [x]))
p_func = eval(build_function(p, [x]))

continuity_neutrals = Dt(ρn) + Dx(ρn*un) + ne * nn * k_ionization(ϵ)
continuity_ions = Dt(ρi) + Dx(ρiui) - ne * nn * k_ionization(ϵ)

# Test that conservative and non-conservative forms are both satisfied, as well as both coupled and uncoupled
momentum_conservative_uncoupled = Dt(ρiui) + Dx(ρiui^2 / ρi + p) + e * ne * ∇ϕ
momentum_conservative_coupled = Dt(ρiui) + Dx(ρiui^2 / ρi + p + nϵ) + e * ue / μ
momentum_nonconservative_uncoupled = ρi * (Dt(ui) + ui * Dx(ui) + nn * ui * k_ionization(ϵ)) + ∇p + e * ne * ∇ϕ
momentum_nonconservative_coupled = ρi * (Dt(ui) + ui * Dx(ui) + nn * ui * k_ionization(ϵ)) + ∇p + ∇pe + e * ue / μ

source_ρn = eval(build_function(expand_derivatives(continuity_neutrals), [x]))
source_ρi = eval(build_function(expand_derivatives(continuity_ions), [x]))
source_ρiui_conservative_uncoupled = eval(build_function(expand_derivatives(momentum_conservative_uncoupled), [x]))
source_ρiui_conservative_coupled = eval(build_function(expand_derivatives(momentum_conservative_coupled), [x]))
source_ρiui_nonconservative_uncoupled = eval(build_function(expand_derivatives(momentum_nonconservative_uncoupled), [x]))
source_ρiui_nonconservative_coupled = eval(build_function(expand_derivatives(momentum_nonconservative_coupled), [x]))

function solve_ions(ncells, source_momentum, scheme, plot_results = false)

    grid = HallThruster.generate_grid(HallThruster.SPT_100, ncells)

    propellant = HallThruster.Xenon

    config = (;
        source_neutrals = (U, p, i) -> source_ρn(p.z_cell[i]),
        source_ion_continuity = (
            (U, p, i) -> source_ρi(p.z_cell[i])
        ),
        source_ion_momentum = (
            (U, p, i) -> source_momentum(p.z_cell[i])
        ),
        propellant,
        ncharge = 1,
        electron_pressure_coupled = coupled,
        min_electron_temperature = 1.0,
    )

    z_edge = grid.edges
    z_cell = grid.cell_centers

    ue = ue_func.(z_cell)
    μ = μ_func.(z_cell)
    nϵ = nϵ_func.(z_cell)
    ρn_exact = ρn_func.(z_cell)
    ρi_exact = ρi_func.(z_cell)
    ρiui_exact = ρiui_func.(z_cell)

    cache = (;ue, μ)

    reactions = HallThruster.ionization_fits_Xe()
    U = zeros(4, ncells)
    U[index.ρn, :] .= ρn_func(0.0)
    U[index.ρi[1], :] .= ρi_func(0.0)
    U[index.ρiui[1], :] .= ρiui_func(0.0)
    U[index.nϵ, :] .= nϵ

    species, fluids, fluid_ranges, species_range_dict = configure_fluids(config)
    index = configure_index(fluid_ranges)

    params = (;
        index,
        config,
        cache,
        scheme,
        fluids,
        species_range_dict,
        reactions,
        z_edge,
        z_cell,
    )

    tspan = (0, 1e-3)
    dt = 1e-8

    alg = Tsit5()

    prob = ODEProblem{true}(HallThruster.update_heavy_species!, U, tspan, params)
	sol = solve(
        prob, alg; saveat=tspan, adaptive=true, dt=dt, maxiters = 1_000_000,
    )

    ρn_sim = sol.u[end][index.ρn, :]
    ρi_sim = sol.u[end][index.ρi[1], :]
    ρiui_sim = sol.u[end][index.ρiui[1], :]

    if plot_results

        p1 = plot(z_cell, ρn_sim, label = "sim", title = "Neutral density")
        plot!(p1, z_cell, ρn_exact, label = "exact")

        p2 = plot(z_cell, ρi_sim, label = "sim", title = "Ion density")
        plot!(p2, z_cell, ρi_exact, label = "exact")

        p3 = plot(z_cell, ρiui_sim, label = "sim", title = "Ion momentum")
        plot!(p3, z_cell, ρiui_exact, label = "exact")

        p = plot(p1, p2, p3, layout = (1, 3), size = (1800, 600))
        display(p)
    end

    return (
        ρn = (z_cell, ρn_sim, ρn_exact),
        ρi = (z_cell, ρi_sim, ρi_exact),
        ρiui = (z_cell, ρiui_sim, ρiui_exact),
    )
end

end