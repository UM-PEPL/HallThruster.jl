module OVS_Ions

include("ovs_funcs.jl")

using HallThruster: HallThruster as het
using LinearAlgebra
using Symbolics
using Accessors

struct R <: het.Reaction end

@variables x t

Dt = Differential(t)
Dx = Differential(x)

k_ionization(ϵ) = 0.0

const un = 1000
const mi = het.Xenon.m
const e = het.e
const Ti = 300
const L = 0.05

ϕ = sin_wave(x / L, amplitude = 5, phase = π / 2, nwaves = 0.5)
ne = sin_wave(x / L, amplitude = 1.0e13, phase = π / 2, nwaves = 2, offset = 6.0e13)
nn = sin_wave(x / L, amplitude = 2.0e18, phase = π / 2, nwaves = 1, offset = 6.0e18)
ui = sin_wave(x / L, amplitude = 2000, phase = -π / 2, nwaves = 0.5, offset = 3000)
ϵ = sin_wave(x / L, amplitude = 3, phase = -π / 2, nwaves = 1, offset = 6)
nϵ = ne * ϵ
∇ϕ = Dx(ϕ)
ρi = ne * mi
ρiui = ρi * ui
ρn = nn * mi
p = ne * het.kB * Ti

ϕ_func = eval(build_function(ϕ, [x]))
ne_func = eval(build_function(ne, [x]))
ρiui_func = eval(build_function(ρiui, [x]))
nϵ_func = eval(build_function(nϵ, [x]))
∇ϕ_func = eval(build_function(expand_derivatives(∇ϕ), [x]))
ρn_func = eval(build_function(ρn, [x]))
ρi_func = eval(build_function(ρi, [x]))
p_func = eval(build_function(p, [x]))
ui_func = eval(build_function(ui, [x]))
ϵ_func = eval(build_function(ϵ, [x]))

continuity_neutrals = Dt(ρn) + un * Dx(ρn) + ne * ρn * k_ionization(ϵ)
continuity_ions = Dt(ρi) + Dx(ρiui) - ne * ρn * k_ionization(ϵ)
momentum_ions = Dt(ρiui) + Dx(ρi * ui^2 + p) + e * ne * ∇ϕ

source_ρn = eval(build_function(expand_derivatives(continuity_neutrals), [x]))
source_ρi = eval(build_function(expand_derivatives(continuity_ions), [x]))
source_ρiui = eval(build_function(expand_derivatives(momentum_ions), [x]))

function solve_ions(ncells, reconstruct; t_end = 1.0e-4)
    # Create config struct
    thruster = het.SPT_100
    A_ch = het.channel_area(thruster)
    anode_mass_flow_rate = un * (ρn_func(0.0) + ui_func(0.0) / un * ρi_func(0.0)) * A_ch
    domain = (0.0, 0.05)

    function source_heavy_species(fluid_containers, params)
        (; grid) = params
        (; continuity, isothermal) = fluid_containers

        for i in 2:(length(grid.cell_centers) - 1)
            continuity[1].dens_ddt[i] += source_ρn(grid.cell_centers[i])
            isothermal[1].dens_ddt[i] += source_ρi(grid.cell_centers[i])
            isothermal[1].mom_ddt[i] += source_ρiui(grid.cell_centers[i])
        end
        return
    end

    config = het.Config(;
        domain,
        discharge_voltage = 0.0,
        thruster,
        propellant = het.Xenon,
        ncharge = 1,
        neutral_velocity = un,
        neutral_temperature_K = 300.0,
        ion_temperature_K = Ti,
        anode_mass_flow_rate,
        reconstruct,
        ionization_model = :OVS,
        LANDMARK = true,
        conductivity_model = het.LANDMARK_conductivity(),
        ion_wall_losses = false,
        anode_boundary_condition = :dirichlet,
        anom_model = het.NoAnom(),
        source_heavy_species = source_heavy_species,
    )

    simparams = het.SimParams(
        grid = het.UnevenGrid(ncells),
        dt = 1.0e-9,
        duration = 1.0,
    )

    params = het.setup_simulation(config, simparams)
    # Need to overwrite min_Te b/c here Te can be lower than Te_L or Te_R
    z_cell = params.grid.cell_centers
    z_start = z_cell[1]
    z_end = z_cell[end]
    @reset params.min_Te = 0.01 * 2 / 3 * min(ϵ_func(z_start), ϵ_func(z_end))
    (; cache, grid) = params

    # Allocate arrays and fill variables
    ρn_exact = ρn_func.(z_cell)
    ρi_exact = ρi_func.(z_cell)
    ui_exact = ui_func.(z_cell)
    @. cache.nϵ = nϵ_func.(z_cell)
    @. cache.ϵ = ϵ_func.(z_cell)
    @. cache.∇ϕ = ∇ϕ_func.(z_cell)
    @. cache.channel_area = A_ch

    # Set up initial condition
    line(v0, v1, z) = v0 + (v1 - v0) * (z - z_start) / (z_end - z_start)
    ρn_0 = [line(ρn_func(z_start), ρn_func(z_end), z) for z in z_cell]
    ρi_0 = [line(ρi_func(z_start), ρi_func(z_end), z) for z in z_cell]
    ρiui_0 = ρi_0 * ui_func(0.0)
    fluids = params.fluid_containers
    fluids.continuity[1].density .= ρn_0
    fluids.isothermal[1].density .= ρi_0
    fluids.isothermal[1].momentum .= ρiui_0

    # Compute timestep
    amax = maximum(abs.(ui_exact) .+ sqrt.(2 / 3 * e * ϵ_func.(z_cell) / mi))
    cache.dt[] = 0.1 * minimum(grid.dz_cell) / amax

    t = 0.0
    while t < t_end
        # Neumann BC
        fluids.continuity[1].density[end] = fluids.continuity[1].density[end - 1]
        fluids.isothermal[1].density[end] = fluids.isothermal[1].density[end - 1]
        fluids.isothermal[1].momentum[end] = fluids.isothermal[1].momentum[end - 1]
        het.step_heavy_species!(fluids, params, source_heavy_species, cache.dt[])
        t += cache.dt[]
    end

    ρn_sim = fluids.continuity[1].density
    ρi_sim = fluids.isothermal[1].density
    ρiui_sim = fluids.isothermal[1].momentum
    ui_sim = ρiui_sim ./ ρi_sim

    return (
        ρn = (; z_cell, sim = ρn_sim, exact = ρn_exact),
        ρi = (; z_cell, sim = ρi_sim, exact = ρi_exact),
        ui = (; z_cell, sim = ui_sim, exact = ui_exact),
    )
end

end
