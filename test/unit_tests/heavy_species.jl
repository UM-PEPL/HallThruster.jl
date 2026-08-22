using HallThruster: HallThruster as het

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
