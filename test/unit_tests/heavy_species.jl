using HallThruster: HallThruster as het

@testset "Primitive velocity cache" begin
    # Primitive velocity is derived once per RK stage and stored with its fluid;
    # density and momentum remain the authoritative conservative state.
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    grid = het.Grid1D(range(0.0, 1.0; length = 4))
    ion = only(het.allocate_fluids(propellant, grid).isothermal)
    ion.density .= [0.0, 2.0, 4.0, 5.0, 0.0]
    ion.momentum .= [1.0, 6.0, -8.0, 0.0, -1.0]

    @test het.update_primitive_velocity!(ion) === nothing
    @test ion.vel_prim == [0.0, 3.0, -2.0, 0.0, 0.0]
end

@testset "Profiled neutral flux" begin
    grid = het.Grid1D([0.0, 0.04, 0.08])
    propellant = het.Propellant(
        het.Xenon, 5.0e-6;
        velocity_m_s = het.LinearInterpolation([0.0, 0.08], [100.0, 300.0]),
        temperature_K = het.LinearInterpolation([0.0, 0.08], [400.0, 800.0]),
    )
    neutral = only(het.allocate_fluids(propellant, grid).continuity)
    neutral.dens_L .= [1.0, 2.0, 3.0]
    neutral.dens_R .= [4.0, 5.0, 6.0]

    het.compute_fluxes_continuity!(neutral, grid)

    expected_flux = @. 0.5 * (
        neutral.vel_L * (neutral.dens_L + neutral.dens_R) -
            neutral.wave_speed * (neutral.dens_R - neutral.dens_L)
    )
    @test neutral.flux_dens ≈ expected_flux
    @test neutral.max_timestep[] ≈ minimum(grid.dz_edge ./ neutral.wave_speed)
end

@testset "Isothermal ion flux" begin
    grid = het.Grid1D([0.0, 0.4, 1.0])
    temperature = 900.0
    propellant = het.Propellant(
        het.Xenon, 0.0; max_charge = 1, ion_temperature_K = temperature,
    )
    ion = only(het.allocate_fluids(propellant, grid).isothermal)
    RT = het.R0 / het.Xenon.M * temperature

    # For p = ρRT, both the pressure flux and the eigenspeeds must use the
    # isothermal sound speed sqrt(RT), without an adiabatic gamma factor.
    @test ion.sound_speed ≈ sqrt(RT)

    ion.dens_L .= [1.0, 2.0, 3.0]
    ion.dens_R .= [1.5, 2.5, 3.5]
    ion.vel_L .= [-2.0, 1.0, 4.0]
    ion.vel_R .= [-1.0, 3.0, 2.0]
    het.compute_fluxes_isothermal!(ion, grid)

    smax = @. max(abs(ion.vel_L), abs(ion.vel_R)) + sqrt(RT)
    expected_density_flux = @. 0.5 * (
        ion.dens_L * ion.vel_L + ion.dens_R * ion.vel_R -
            smax * (ion.dens_R - ion.dens_L)
    )
    expected_momentum_flux = @. 0.5 * (
        ion.dens_L * (ion.vel_L^2 + RT) +
            ion.dens_R * (ion.vel_R^2 + RT) -
            smax * (ion.dens_R * ion.vel_R - ion.dens_L * ion.vel_L)
    )

    @test ion.flux_dens ≈ expected_density_flux
    @test ion.flux_mom ≈ expected_momentum_flux
    @test ion.max_timestep[] ≈ minimum(grid.dz_edge ./ smax)
end

@testset "Zero-density heavy species" begin
    @test het.primitive_velocity(0.0, 0.0) == 0.0
    @test isfinite(het.primitive_velocity(0.0, 0.0))
    @test het.primitive_velocity(6.0, 2.0) == 3.0

    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluids = het.allocate_fluids(propellant, grid)
    fluid_arr = [fluids.continuity[1], fluids.isothermal[1]]

    rxn = het.ElectronImpactReaction(0.0, het.Xenon(0), [het.Xenon(1)], ones(256))
    groups = het.build_electron_impact_groups([rxn], [1], [[2]], fluid_arr)
    num_grid_cells = ncells + 2
    fluid_arr[2].density .= propellant.gas.m * 1.0e18
    cache = (;
        inelastic_losses = zeros(num_grid_cells),
        νiz = zeros(num_grid_cells),
        νex_explicit = zeros(num_grid_cells),
        ϵ = zeros(num_grid_cells),
        ne = zeros(num_grid_cells),
        K = zeros(num_grid_cells),
        nϵ = fill(1.0e19, num_grid_cells),
        cell_cache_1 = zeros(num_grid_cells),
        reaction_loss_frequency = zeros(num_grid_cells),
        reaction_rate_indices = zeros(Int, num_grid_cells),
        reaction_rate_fractions = zeros(num_grid_cells),
        reaction_rate_index_limit = [254],
        dt_iz = [Inf],
    )

    het.apply_reaction_groups!(fluid_arr, groups, cache, false)

    @test cache.dt_iz[] == Inf
    @test all(isfinite, cache.νiz)
    @test all(isfinite, cache.inelastic_losses)
    @test all(iszero, cache.νiz)
    @test all(iszero, cache.νex_explicit)
    @test all(iszero, cache.inelastic_losses)
    @test all(iszero, fluid_arr[1].dens_ddt)
    @test all(iszero, fluid_arr[2].dens_ddt)

    # A populated reactant produces the expected grouped chemistry timestep.
    fluid_arr[1].density[2:(end - 1)] .= propellant.gas.m
    het.apply_reaction_groups!(fluid_arr, groups, cache, false)
    @test cache.dt_iz[] ≈ 1.0e-18
end

@testset "Empty heavy-species populations" begin
    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluids = het.allocate_fluids(propellant, grid)
    cache = het.allocate_arrays(length(first(fluids.continuity).density), 0)

    # User-provided initial conditions and restarts may contain an empty cell
    # before the normal density limiter runs; derived plasma fields must stay finite.
    het.update_heavy_species_cache!(fluids, cache, false)
    @test all(isfinite, cache.avg_neutral_vel)
    @test all(isfinite, cache.avg_ion_vel)
    @test all(isfinite, cache.m_eff)
    @test all(isfinite, cache.Z_eff)
    @test all(iszero, cache.avg_neutral_vel)
    @test all(iszero, cache.avg_ion_vel)
    @test all(==(propellant.gas.m), cache.m_eff)
    @test all(==(1.0), cache.Z_eff)

    # A configuration with no charged fluids cannot provide the ion properties
    # required by the electron and wall models, so fail with a useful error.
    neutral_only = (; continuity = fluids.continuity, isothermal = typeof(fluids.isothermal)())
    @test_throws ArgumentError het.update_heavy_species_cache!(
        neutral_only, cache, false,
    )
end

@testset "Ion acceleration timestep" begin
    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 3)
    fluid_grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluids = het.allocate_fluids(propellant, fluid_grid).isothermal
    grid = (; dz_cell = [1.0, 0.4, 0.7, 0.5, 1.0])
    cache = (; ∇ϕ = [0.0, -2.0, 4.0, -8.0, 0.0], dt_E = fill(0.0))

    for (j, fluid) in enumerate(fluids)
        fluid.density .= j
    end

    het.apply_ion_acceleration!(fluids, grid, cache)

    expected_dt = Inf
    for fluid in fluids
        qe_m = fluid.species.Z * het.e / fluid.species.element.m
        for i in 2:(ncells + 1)
            qE_m = -qe_m * cache.∇ϕ[i]
            expected_dt = min(expected_dt, abs(grid.dz_cell[i] / qE_m))
            @test fluid.mom_ddt[i] == qE_m * fluid.density[i]
        end
    end

    @test cache.dt_E[] == sqrt(expected_dt)
end
