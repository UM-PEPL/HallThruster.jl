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
        inelastic_losses,
        false,
    )

    @test dt_max == Inf
    @test all(isfinite, νiz)
    @test all(isfinite, inelastic_losses)
    @test all(iszero, νiz)
    @test all(iszero, inelastic_losses)
    @test all(iszero, fluid_arr[1].dens_ddt)
    @test all(iszero, fluid_arr[2].dens_ddt)
    @test all(iszero, rxn_cache[1])
    @test all(iszero, rxn_cache[2])

    # For first-order electron-impact destruction, ρ / ρdot = 1 / (k * ne).
    # The reaction kernel uses this frequency form to avoid a division per cell.
    fluid_arr[1].density[2:(end - 1)] .= 2.0
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
        inelastic_losses,
        false,
    )
    @test dt_max ≈ inv(ne[2])
    expected_νiz = 2.0 / fluid_arr[1].species.element.m
    @test all(νiz[2:(end - 1)] .≈ expected_νiz)
end
