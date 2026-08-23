using HallThruster: HallThruster as het

@testset "Primitive velocity cache" begin
    # Primitive velocity is derived once per RK stage and stored with its fluid;
    # density and momentum remain the authoritative conservative state.
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    ion = only(het.allocate_fluids(propellant, 3).isothermal)
    ion.density .= [0.0, 2.0, 4.0, 5.0, 0.0]
    ion.momentum .= [1.0, 6.0, -8.0, 0.0, -1.0]

    @test het.update_primitive_velocities!([ion]) === nothing
    @test ion.vel_prim == [0.0, 3.0, -2.0, 0.0, 0.0]
end

@testset "Zero-density heavy species" begin
    @test het.primitive_velocity(0.0, 0.0) == 0.0
    @test isfinite(het.primitive_velocity(0.0, 0.0))
    @test het.primitive_velocity(6.0, 2.0) == 3.0

    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    fluids = het.allocate_fluids(propellant, ncells)
    fluid_arr = [fluids.continuity[1], fluids.isothermal[1]]

    rxn = het.ElectronImpactReaction(0.0, het.Xenon(0), [het.Xenon(1)], ones(256))
    rxn_cache = (zeros(ncells + 2), zeros(ncells + 2))
    ne = fill(1.0e18, ncells + 2)
    energy = fill(10.0, ncells + 2)
    νiz = zeros(ncells + 2)
    νex_explicit = zeros(ncells + 2)
    inelastic_losses = zeros(ncells + 2)

    dt_max = het.apply_reaction!(
        fluid_arr,
        1,
        [2],
        rxn.product_coeffs,
        rxn_cache,
        ne,
        energy,
        rxn,
        νiz,
        νex_explicit,
        inelastic_losses,
        false,
    )

    @test dt_max == Inf
    @test all(isfinite, νiz)
    @test all(isfinite, inelastic_losses)
    @test all(iszero, νiz)
    @test all(iszero, νex_explicit)
    @test all(iszero, inelastic_losses)
    @test all(iszero, fluid_arr[1].dens_ddt)
    @test all(iszero, fluid_arr[2].dens_ddt)
    @test all(iszero, rxn_cache[1])
    @test all(iszero, rxn_cache[2])

    # Standalone callers without the summed-loss cache still receive their
    # per-reaction timestep limit.
    fluid_arr[1].density[2:(end - 1)] .= propellant.gas.m
    dt_max = het.apply_reaction!(
        fluid_arr,
        1,
        [2],
        rxn.product_coeffs,
        rxn_cache,
        ne,
        energy,
        rxn,
        νiz,
        νex_explicit,
        inelastic_losses,
        false,
    )
    @test dt_max ≈ 1.0e-18
end

@testset "Empty heavy-species populations" begin
    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    fluids = het.allocate_fluids(propellant, ncells)
    cache = het.allocate_arrays(length(first(fluids.continuity).density), [propellant], 0)

    # User-provided initial conditions and restarts may contain an empty cell
    # before the normal density limiter runs; derived plasma fields must stay finite.
    het.update_heavy_species_cache!(fluids, cache, nothing, false)
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
        neutral_only, cache, nothing, false,
    )
end

@testset "Ion acceleration timestep" begin
    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 3)
    fluids = het.allocate_fluids(propellant, ncells).isothermal
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
