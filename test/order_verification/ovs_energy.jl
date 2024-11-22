module OVS_Energy

include("ovs_funcs.jl")

using Symbolics, HallThruster, LinearAlgebra

struct R <: HallThruster.Reaction end

@variables x t

Dt = Differential(t)
Dx = Differential(x)

L = 0.05

ϕ = sin_wave(x / L, amplitude = 300, phase = π / 2, nwaves = 0.25)
ne = sin_wave(x / L, amplitude = 2e15, phase = π / 4, nwaves = 0.5, offset = 1.1e16)
nn = sin_wave(x / L, amplitude = 5e18, phase = pi / 3, nwaves = 2.0, offset = 6e18)
ui = sin_wave(x / L, amplitude = 13000, phase = π / 4, nwaves = 0.75, offset = 10000)
μ = sin_wave(x / L, amplitude = 1e4, phase = π / 2, nwaves = 1.2, offset = 1.1e4)
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

k(ϵ) = 8.32 * HallThruster.rate_coeff(OVS_Excitation(), R(), ϵ)
W(ϵ) = 1e7 * ϵ * exp(-20 / ϵ)
energy_eq = Dt(nϵ) + Dx(5 / 3 * nϵ * ue - 10 / 9 * μ * nϵ * Dx(nϵ / ne)) +
            ne * (-ue * Dx(ϕ) + nn * k(ϵ) + W(ϵ))
source_energy = eval(build_function(expand_derivatives(energy_eq), [x]))

function solve_energy!(params, max_steps, dt, rtol = sqrt(eps(Float64)))
    t = 0.0
    nϵ_old = copy(params.cache.nϵ)
    residual = Inf
    iter = 0
    res0 = 0.0
    while iter < max_steps && abs(residual / res0) > rtol
        HallThruster.update_electron_energy!(params, dt)
        params.cache.νiz .= 0.0
        params.cache.νex .= 0.0
        params.cache.inelastic_losses .= 0.0
        params.config.conductivity_model(params.cache.κ, params) # update thermal conductivity
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
    grid = HallThruster.generate_grid(
        HallThruster.SPT_100.geometry, (0.0, 0.05), UnevenGrid(ncells),)
    z_cell = grid.cell_centers
    z_edge = grid.edges
    dz_cell = grid.dz_cell
    dz_edge = grid.dz_edge
    ncells = length(z_cell)

    # fill cache values
    _, cache = HallThruster.allocate_arrays(
        grid, (; ncharge = 1, anom_model = HallThruster.NoAnom()),)
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

    Te_L = nϵ_exact[1] / cache.ne[1]
    Te_R = nϵ_exact[end] / cache.ne[end]
    @. cache.nϵ = Te_L * cache.ne

    config = (;
        anode_Te = 2 / 3 * Te_L,
        cathode_Te = 2 / 3 * Te_R,
        ncharge = 1,
        source_energy = (params, i) -> source_energy(params.z_cell[i]),
        min_electron_temperature = 0.1 * min(Te_L, Te_R),
        transition_length = 0.0,
        LANDMARK = true,
        propellant = Xenon,
        ionization_model = OVS_Ionization(), excitation_model = OVS_Excitation(),
        wall_loss_model = HallThruster.ConstantSheathPotential(20.0, 1.0, 1.0),
        geometry = (; channel_length = 0.025),
        anode_boundary_condition = :dirichlet,
        conductivity_model = HallThruster.LANDMARK_conductivity(),
    )

    species = [HallThruster.Xenon(0), HallThruster.Xenon(1)]
    species_range_dict = Dict([:Xe => 1, Symbol("Xe+") => 0])

    ionization_reactions = HallThruster._load_reactions(config.ionization_model, species)
    ionization_reactant_indices = HallThruster.reactant_indices(
        ionization_reactions, species_range_dict,)
    ionization_product_indices = HallThruster.product_indices(
        ionization_reactions, species_range_dict,)

    excitation_reactions = HallThruster._load_reactions(config.excitation_model, species)
    excitation_reactant_indices = HallThruster.reactant_indices(
        excitation_reactions, species_range_dict,)

    dt = 8 / maximum(abs.(cache.ue)) * (z_cell[2] - z_cell[1])
    params_base = (;
        dt,
        z_cell, z_edge,
        min_Te = 0.1 * min(Te_L, Te_R),
        cache = deepcopy(cache),
        L_ch = config.geometry.channel_length,
        ionization_reactions,
        ionization_reactant_indices,
        ionization_product_indices,
        excitation_reactions,
        excitation_reactant_indices,
        Δz_cell = dz_cell, Δz_edge = dz_edge,
        ncells,
        grid,
    )

    # Test backward euler implicit solve
    params_implicit = (; params_base..., config = (; config..., implicit_energy = 1.0))
    solve_energy!(params_implicit, niters, dt)
    results_implicit = (; z = z_cell, exact = nϵ_exact, sim = params_implicit.cache.nϵ[:])

    # Test crank-nicholson implicit solve
    params_cn = (; params_base..., config = (; config..., implicit_energy = 0.5))
    solve_energy!(params_cn, niters, dt)
    results_cn = (; z = z_cell, exact = nϵ_exact, sim = params_cn.cache.nϵ[:])

    return (results_implicit, results_cn)
end

end
