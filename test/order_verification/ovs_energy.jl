module OVS_Energy

include("ovs_funcs.jl")

using Symbolics, LinearAlgebra
using HallThruster: HallThruster as het

struct R <: het.Reaction end

@variables x t

Dt = Differential(t)
Dx = Differential(x)

L = 0.05

ϕ = sin_wave(x / L, amplitude = 300, phase = π / 2, nwaves = 0.25)
ne = sin_wave(x / L, amplitude = 2.0e15, phase = π / 4, nwaves = 0.5, offset = 1.1e16)
nn = sin_wave(x / L, amplitude = 5.0e18, phase = pi / 3, nwaves = 2.0, offset = 6.0e18)
ui = sin_wave(x / L, amplitude = 13000, phase = π / 4, nwaves = 0.75, offset = 10000)
μ = sin_wave(x / L, amplitude = 1.0e4, phase = π / 2, nwaves = 1.2, offset = 1.1e4)
ϵ = sin_wave(x / L, amplitude = 20, phase = 1.3 * π / 2, nwaves = 1.1, offset = 30)
∇ϕ = Dx(ϕ)
niui = ne * ui
nϵ = ne * ϵ
ue = μ * (∇ϕ - Dx(nϵ) / ne)
κ = 10 / 9 * μ * nϵ

ϕ_func = eval(build_function(ϕ, [x]))
ne_func = eval(build_function(ne, [x]))
μ_func = eval(build_function(μ, [x]))
niui_func = eval(build_function(niui, [x]))
nϵ_func = eval(build_function(nϵ, [x]))
κ_func = eval(build_function(κ, [x]))
ue_func = eval(build_function(expand_derivatives(ue), [x]))
∇ϕ_func = eval(build_function(expand_derivatives(∇ϕ), [x]))
nn_func = eval(build_function(nn, [x]))
ϵ_func = eval(build_function(ϵ, [x]))

k(ϵ) = 8.32 * het.ovs_rate_coeff_ex(ϵ)
W(ϵ) = 1.0e7 * ϵ * exp(-20 / ϵ)
energy_eq = Dt(nϵ) + Dx(5 / 3 * nϵ * ue - 10 / 9 * μ * nϵ * Dx(nϵ / ne)) +
    ne * (-ue * Dx(ϕ) + nn * k(ϵ) + W(ϵ))
source_energy = eval(build_function(expand_derivatives(energy_eq), [x]))

function solve_energy!(params, config, max_steps, dt, rtol = sqrt(eps(Float64)))
    t = 0.0
    nϵ_old = copy(params.cache.nϵ)
    residual = Inf
    iter = 0
    res0 = 0.0
    while iter < max_steps && abs(residual / res0) > rtol
        het.update_electron_energy!(
            params, config.wall_loss_model, config.source_energy, dt,
        )
        params.cache.νiz .= 0.0
        params.cache.νex .= 0.0
        params.cache.inelastic_losses .= 0.0
        config.conductivity_model(params.cache.κ, params) # update thermal conductivity
        residual = Lp_norm(params.cache.nϵ .- nϵ_old, 2)
        if iter == 1
            res0 = residual
        end
        nϵ_old .= params.cache.nϵ
        t += dt
        iter += 1
    end

    return params
end

function verify_energy(ncells; niters = 20000)
    grid = het.generate_grid(het.UnevenGrid(ncells), het.SPT_100.geometry, (0.0, 0.05))
    z_cell = grid.cell_centers
    ncells = length(z_cell)

    Te_L = nϵ_func(z_cell[1]) / ne_func(z_cell[1])
    Te_R = nϵ_func(z_cell[end]) / ne_func(z_cell[end])

    config = (;
        domain = (0.0, 1.0),
        discharge_voltage = 0.0,
        anode_mass_flow_rate = 0.0,
        anode_Tev = 2 / 3 * Te_L,
        cathode_Tev = 2 / 3 * Te_R,
        ncharge = 1,
        source_energy = (params, i) -> source_energy(params.grid.cell_centers[i]),
        transition_length = 0.0,
        LANDMARK = true,
        propellant = het.Xenon,
        ionization_model = :OVS, excitation_model = :OVS,
        wall_loss_model = het.ConstantSheathPotential(20.0, 1.0, 1.0),
        thruster = het.SPT_100,
        anode_boundary_condition = :dirichlet,
        conductivity_model = het.LANDMARK_conductivity(),
        electron_plume_loss_scale = 1.0,
        implicit_energy = 1.0,
        anom_model = het.NoAnom(),
    )

    cfg_implicit = het.Config(; config..., implicit_energy = 1.0)
    cfg_cn = het.Config(; config..., implicit_energy = 1.0)

    species = [het.Xenon(0), het.Xenon(1)]
    species_range_dict = Dict([:Xe => 1, Symbol("Xe+") => 0])

    ionization_reactions = het.load_ionization_reactions(config.ionization_model, species)
    ionization_reactant_indices = het.reactant_indices(
        ionization_reactions, species_range_dict,
    )
    ionization_product_indices = het.product_indices(
        ionization_reactions, species_range_dict,
    )

    excitation_reactions = het.load_excitation_reactions(config.excitation_model, species)
    excitation_reactant_indices = het.reactant_indices(
        excitation_reactions, species_range_dict,
    )
    # fill cache values
    _, cache = het.allocate_arrays(grid, cfg_implicit)
    @. cache.μ = μ_func(z_cell)
    @. cache.κ = κ_func(z_cell)
    @. cache.ne = ne_func(z_cell)
    @. cache.ue = ue_func(z_cell)
    @. cache.∇ϕ = ∇ϕ_func(z_cell)
    @. cache.nn = nn_func(z_cell)
    @. cache.niui = niui_func(z_cell)'
    @. cache.Tev = 2 / 3 * ϵ_func(z_cell)
    @. cache.channel_area = 1.0
    @. cache.ni = cache.ne'

    nϵ_exact = nϵ_func.(z_cell)
    @. cache.pe = copy(nϵ_exact)
    @. cache.nϵ = copy(nϵ_exact)

    params_base = (;
        grid,
        min_Te = 0.01 * min(Te_L, Te_R),
        cache = deepcopy(cache),
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
    )

    params_implicit = (; params_base..., het.params_from_config(cfg_implicit)...)
    params_cn = (; params_base..., het.params_from_config(cfg_cn)...)

    dt = 8 / maximum(abs.(cache.ue)) * (z_cell[2] - z_cell[1])

    # Test backward euler implicit solve
    solve_energy!(params_implicit, cfg_implicit, niters, dt)
    results_implicit = (; z = z_cell, exact = nϵ_exact, sim = params_implicit.cache.nϵ[:])

    # Test crank-nicholson implicit solve
    solve_energy!(params_cn, cfg_cn, niters, dt)
    results_cn = (; z = z_cell, exact = nϵ_exact, sim = params_cn.cache.nϵ[:])

    return (results_implicit, results_cn)
end

end
