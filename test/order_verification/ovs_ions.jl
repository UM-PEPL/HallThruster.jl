module OVS_Ions

include("ovs_funcs.jl")

using Symbolics, HallThruster, LinearAlgebra

@variables x t

Dt = Differential(t)
Dx = Differential(x)

const k_ionization = OVS_rate_coeff_iz
const un = 1000
const mi = HallThruster.Xenon.m
const e = HallThruster.e
const Ti = 300
const L = 0.05

ϕ = 0.0 * sin_wave(x/L, amplitude = 300, phase = π/2, nwaves = 0.5)
ne = sin_wave(x/L, amplitude = 1e13, phase = π/2, nwaves = 2, offset = 6e13)
nn = sin_wave(x/L, amplitude = 2e18, phase = π/2, nwaves = 1, offset = 6e18)
ui = sin_wave(x/L, amplitude = 2000, phase = -π/2, nwaves = 0.5, offset = 3000)
μ = sin_wave(x/L, amplitude = 1e2, phase = 3π/2, nwaves = 0.6, offset = 1.1e2)
ϵ = sin_wave(x/L, amplitude = 3, phase = -π/2, nwaves = 1, offset = 6)
nϵ = ne * ϵ
pe = e * nϵ
∇ϕ = Dx(ϕ)
∇pe = Dx(pe)
ρi = ne * mi
ρiui = ρi * ui
ρn = nn * mi
ue = μ * (∇ϕ - Dx(nϵ)/ne)
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
ui_func = eval(build_function(ui, [x]))
ϵ_func = eval(build_function(ϵ, [x]))

continuity_neutrals = Dt(ρn) + un * Dx(ρn) + ne * ρn * k_ionization(ϵ)
continuity_ions = Dt(ρi) + Dx(ρiui) - ne * ρn * k_ionization(ϵ)

# Test that conservative and non-conservative forms are both satisfied, as well as both coupled and uncoupled
momentum_conservative_uncoupled = Dt(ρiui) + Dx(ρi * ui^2 + p) + e * ne * ∇ϕ
momentum_conservative_coupled = Dt(ρiui) + Dx(ρi * ui^2 + p + pe) + e * ne * ue / μ
momentum_nonconservative_uncoupled = ρi * (Dt(ui) + ui * Dx(ui) + nn * ui * k_ionization(ϵ)) + ∇p + e * ne * ∇ϕ
momentum_nonconservative_coupled = ρi * (Dt(ui) + ui * Dx(ui) + nn * ui * k_ionization(ϵ)) + ∇p + ∇pe + e * ne * ue / μ

source_ρn = eval(build_function(expand_derivatives(continuity_neutrals), [x]))
source_ρi = eval(build_function(expand_derivatives(continuity_ions), [x]))
source_ρiui_conservative_uncoupled = eval(build_function(expand_derivatives(momentum_conservative_uncoupled), [x]))
source_ρiui_conservative_coupled   = eval(build_function(expand_derivatives(momentum_conservative_coupled), [x]))
source_ρiui_nonconservative_uncoupled = eval(build_function(expand_derivatives(momentum_nonconservative_uncoupled), [x]))
source_ρiui_nonconservative_coupled   = eval(build_function(expand_derivatives(momentum_nonconservative_coupled), [x]))

function solve_ions(ncells, scheme, plot_results = true; t_end = 1e-4, coupled = true, conservative = true)

    grid = HallThruster.generate_grid(HallThruster.SPT_100.geometry, ncells, (0.0, 0.05))

    Δz_cell, Δz_edge = HallThruster.grid_spacing(grid)
    propellant = HallThruster.Xenon

    source_momentum = if conservative && coupled
        source_ρiui_conservative_coupled
    elseif conservative && !coupled
        source_ρiui_conservative_uncoupled
    elseif !conservative && coupled
        source_ρiui_nonconservative_uncoupled
    elseif !nonconservative  && !coupled
        source_ρiui_nonconservative_uncoupled
    end

    thruster = HallThruster.SPT_100

    A_ch = HallThruster.channel_area(thruster)
    anode_mass_flow_rate = un * (ρn_func(0.0) + ui_func(0.0) / un *ρi_func(0.0))* A_ch

    config = (;
        thruster,
        source_neutrals = ((U, p, i) -> source_ρn(p.z_cell[i]),),
        source_ion_continuity = (
            (U, p, i) -> source_ρi(p.z_cell[i]),
        ),
        source_ion_momentum = (
            (U, p, i) -> source_momentum(p.z_cell[i]),
        ),
        propellant,
        ncharge = 1,
        electron_pressure_coupled = true,
        min_electron_temperature = 1.0,
        neutral_velocity = un,
        neutral_temperature = 300.0,
        ion_temperature = Ti,
        solve_ion_energy = false,
        min_number_density = 1e6,
        anode_sheath = false,
        anode_mass_flow_rate,
        scheme,
        ionization_model = OVS_Ionization(),
        LANDMARK = true,
        ion_wall_losses = false,
        anode_boundary_condition = :dirichlet,
        solve_background_neutrals = false,
    )

    species = [HallThruster.Xenon(0), HallThruster.Xenon(1)]

    ionization_reactions = HallThruster._load_reactions(config.ionization_model, species)

    z_edge = grid.edges
    z_cell = grid.cell_centers
    #z_cell = LinRange(-0.05/ncells/2, 0.05 + 0.05/ncells/2, ncells+2)

    nedges = length(z_edge)

    ue = ue_func.(z_cell)
    μ = μ_func.(z_cell)
    nϵ = nϵ_func.(z_cell)
    ρn_exact = ρn_func.(z_cell)
    ρi_exact = ρi_func.(z_cell)
    ui_exact = ui_func.(z_cell)
    ∇ϕ = ∇ϕ_func.(z_cell)

    fluids, fluid_ranges, species, species_range_dict = HallThruster.configure_fluids(config)
    index = HallThruster.configure_index(fluids, fluid_ranges)

    F = zeros(4, nedges)
    UL = zeros(4, nedges)
    UR = zeros(4, nedges)
    λ_global = zeros(2)

    channel_area = A_ch .* ones(ncells+2)
    dA_dz = zeros(ncells+2)

    cache = (;ue, μ, F, UL, UR, ∇ϕ, λ_global, channel_area, dA_dz)

    U = zeros(4, ncells+2)
    z_end = z_cell[end]
    z_start = z_cell[1]
    line(v0, v1, z) = v0 + (v1 - v0) * (z - z_start) / (z_end - z_start)
    U[index.ρn[1], :] = [line(ρn_func(z_start), ρn_func(z_end), z) for z in z_cell]
    U[index.ρi[1], :] = [line(ρi_func(z_start), ρi_func(z_end), z) for z in z_cell]
    U[index.ρiui[1], :] = U[index.ρi[1], :] * ui_func(0.0)
    U[index.nϵ, :] = nϵ

    params = (;
        index,
        config,
        cache,
        fluids,
        species_range_dict,
        z_edge,
        z_cell,
        Te_L = 2/3 * ϵ_func(z_start),
        Te_R = 2/3 * ϵ_func(z_end),
        A_ch,
        ionization_reactions,
        ionization_reactant_indices = [index.ρn],
        ionization_product_indices = [index.ρi[1]],
        max_timestep = [Inf],
        num_neutral_fluids = 1,
        background_neutral_density = 0.0,
        background_neutral_velocity = 1.0,
        Δz_cell, Δz_edge
    )

    amax = maximum(ui_exact .+ sqrt.(2/3 * e * ϵ_func.(z_cell) / mi))

    tspan = (0, t_end)
    dt = 0.2 * (z_cell[end] - z_cell[1]) / ncells / amax

    dU = copy(U)

    t = 0.0
    while t < tspan[2]
        @views U[:, end] = U[:, end-1]
        HallThruster.update_heavy_species!(dU, U, params, t; apply_boundary_conditions = false)
        for i in eachindex(U)
            U[i] += dt * dU[i]
        end
        t += dt
    end

    sol = (;t = [t], u = [U])

    ρn_sim = sol.u[end][index.ρn[1], :]
    ρi_sim = sol.u[end][index.ρi[1], :]
    ρiui_sim = sol.u[end][index.ρiui[1], :]
    ui_sim = ρiui_sim ./ ρi_sim

    return (
        ρn = (;z_cell, sim = ρn_sim, exact = ρn_exact),
        ρi = (;z_cell, sim = ρi_sim, exact = ρi_exact),
        ui = (;z_cell, sim = ui_sim, exact = ui_exact),
    )
end

end
