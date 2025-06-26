module OVS_Energy

include("ovs_funcs.jl")

using Symbolics, LinearAlgebra
using HallThruster: HallThruster as het
using Accessors

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

function verify_energy(ncells; implicit_energy = 1.0, niters = 20000)
    domain = (0, L)
    simparams = het.SimParams(grid = het.UnevenGrid(ncells))

    Te_L = nϵ_func(domain[1]) / ne_func(domain[1])
    Te_R = nϵ_func(domain[end]) / ne_func(domain[end])

    config = het.Config(;
        domain = domain,
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
        implicit_energy = implicit_energy,
        anom_model = het.NoAnom(),
    )

    _, params = het.setup_simulation(config, simparams)
    @reset params.min_Te = 0.01 * min(Te_L, Te_R)
    z_cell = params.grid.cell_centers

    # fill cache values
    # Need to reallocate cache because setup_simulation does a bunch of initialization that messeṡ this up somehow
    _, cache = het.allocate_arrays(params.grid, config)
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

    dt = 8 / maximum(abs.(cache.ue)) * (z_cell[2] - z_cell[1])
    @reset params.cache = cache

    solve_energy!(params, config, niters, dt)

    return (; z = z_cell, exact = nϵ_exact, sim = params.cache.nϵ[:])
end

function verify_energy_all(ncells; niters = 20000)
    results_implicit = verify_energy(ncells; implicit_energy = 1.0, niters)
    results_cn = verify_energy(ncells; implicit_energy = 0.5, niters)
    return (results_implicit, results_cn)
end

end
