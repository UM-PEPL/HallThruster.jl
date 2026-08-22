using HallThruster: HallThruster as het

@testset "Exact radiative cascade" begin
    # Verify a two-step 2 -> 1 -> 0 cascade against the analytic Bateman solution.
    # Starting with N atoms in level 2 and none in levels 1 or 0, decay rates
    # A21 and A10 give
    #   n2(t) = N exp(-A21 t)
    #   n1(t) = N A21 / (A10 - A21) * (exp(-A21 t) - exp(-A10 t))
    #   n0(t) = N - n1(t) - n2(t).
    # The integrated photon counts are N - n2 for 2 -> 1 and n0 for 1 -> 0.
    ncells = 3
    propellant = het.Propellant(het.Xenon, 1.0e-6; excited_levels = [1, 2])
    fluids = het.allocate_fluids(propellant, ncells).continuity
    fluid_array = collect(fluids)

    rate_21 = 2.0
    rate_10 = 3.0
    reactions = [
        het.DeExcitationReaction(het.Xenon(0, 2), [het.Xenon(0, 1)], [rate_21]),
        het.DeExcitationReaction(het.Xenon(0, 1), [het.Xenon(0)], [rate_10]),
    ]
    reactant_indices = het.reactant_indices(reactions, fluid_array)
    product_indices = het.product_indices(reactions, fluid_array)
    level_energies = Dict(
        het.Xenon(0).symbol => 0.0,
        het.Xenon(0, 1).symbol => 3.0,
        het.Xenon(0, 2).symbol => 5.0,
    )
    networks, emissions, transitions = het.build_radiative_networks(
        fluid_array, reactions, reactant_indices, product_indices, level_energies,
    )
    @test only(networks).transition_energies_eV == [2.0, 3.0]
    @test getfield.(transitions, :energy_eV) == [2.0, 3.0]
    @test getfield.(transitions, :frequency) == [rate_21, rate_10]

    initial_density = 10.0
    mass = propellant.gas.m
    for cell in 2:(ncells + 1)
        fluids[3].density[cell] = mass * initial_density
    end

    dt = 1.0
    het.apply_radiative_decay!(fluid_array, networks, emissions, dt)
    photon_output = het.photon_emissions(transitions, emissions, dt)

    expected_2 = initial_density * exp(-rate_21 * dt)
    expected_1 = initial_density * rate_21 / (rate_10 - rate_21) *
        (exp(-rate_21 * dt) - exp(-rate_10 * dt))
    expected_0 = initial_density - expected_1 - expected_2

    for cell in 2:(ncells + 1)
        populations = [fluid.density[cell] / mass for fluid in fluids]
        @test populations ≈ [expected_0, expected_1, expected_2]
        @test sum(populations) ≈ initial_density
        @test emissions[:, cell] ≈ [initial_density - expected_2, expected_0]
        @test getfield.(photon_output, :emission_rate)[1][cell] ≈
            (initial_density - expected_2) / dt
        @test getfield.(photon_output, :emission_rate)[2][cell] ≈ expected_0 / dt
    end
    @test getfield.(photon_output, :upper) == [
        het.Xenon(0, 2).symbol, het.Xenon(0, 1).symbol,
    ]
    @test getfield.(photon_output, :lower) == [
        het.Xenon(0, 1).symbol, het.Xenon(0).symbol,
    ]
    @test getfield.(photon_output, :frequency) == [rate_21, rate_10]

    # Ghost cells are boundary data and must not be advanced by the source update.
    @test all(iszero, emissions[:, [1, end]])
    @test all(iszero, fluid.density[1] for fluid in fluids)
    @test all(iszero, fluid.density[end] for fluid in fluids)
end

@testset "Radiative branching and ion momentum" begin
    # Exercise competing branches and a downstream cascade in a charged species.
    # The decay should preserve total ion population and axial momentum while
    # recording one photon for every transition that occurs.
    ncells = 1
    propellant = het.Propellant(
        het.Xenon, 1.0e-6;
        max_charge = 1,
        excited_ion_levels = Dict(1 => [1, 2]),
    )
    fluids = het.allocate_fluids(propellant, ncells).isothermal
    fluid_array = collect(fluids)
    reactions = [
        het.DeExcitationReaction(
            het.Xenon(1, 2), [het.Xenon(1, 1), het.Xenon(1)], [2.0, 1.0],
        ),
        het.DeExcitationReaction(het.Xenon(1, 1), [het.Xenon(1)], [4.0]),
    ]
    reactant_indices = het.reactant_indices(reactions, fluid_array)
    product_indices = het.product_indices(reactions, fluid_array)
    level_energies = Dict(
        het.Xenon(1).symbol => 0.0,
        het.Xenon(1, 1).symbol => 3.0,
        het.Xenon(1, 2).symbol => 5.0,
    )
    networks, emissions, transitions = het.build_radiative_networks(
        fluid_array, reactions, reactant_indices, product_indices, level_energies,
    )
    @test only(networks).transition_energies_eV == [2.0, 5.0, 3.0]
    @test getfield.(transitions, :energy_eV) == [2.0, 5.0, 3.0]

    mass = propellant.gas.m
    initial_density = 12.0
    velocity = 4.0e3
    fluids[3].density[2] = mass * initial_density
    fluids[3].momentum[2] = mass * initial_density * velocity

    het.apply_radiative_decay!(fluid_array, networks, emissions, 20.0)

    @test sum(fluid.density[2] for fluid in fluids) / mass ≈ initial_density
    @test sum(fluid.momentum[2] for fluid in fluids) / mass ≈
        initial_density * velocity
    @test emissions[:, 2] ≈ [8.0, 4.0, 8.0] atol = 1.0e-10
end

@testset "Summed electron-impact loss frequency" begin
    # Two reactions consume the same reactant, so the positivity/CFL constraint
    # must use the sum of their loss frequencies rather than either one alone.
    ncells = 3
    propellant = het.Propellant(
        het.Xenon, 1.0e-6; max_charge = 1, excited_levels = [1, 2],
    )
    fluid_set = het.allocate_fluids(propellant, ncells)
    fluids = [fluid_set.continuity; fluid_set.isothermal]
    mass = propellant.gas.m
    fluids[1].density .= mass * 2.0e18
    fluids[4].density .= mass * 1.0e18

    k1 = 2.0e-14
    k2 = 3.0e-14
    reactions = [
        het.ElectronImpactReaction(3.0, het.Xenon(0), [het.Xenon(0, 1)], fill(k1, 256)),
        het.ElectronImpactReaction(5.0, het.Xenon(0), [het.Xenon(0, 2)], fill(k2, 256)),
    ]
    reactant_indices = het.reactant_indices(reactions, fluids)
    product_indices = het.product_indices(reactions, fluids)
    rxns = zip(reactions, reactant_indices, product_indices)

    num_grid_cells = ncells + 2
    cache = (;
        inelastic_losses = zeros(num_grid_cells),
        νiz = zeros(num_grid_cells),
        νex_explicit = zeros(num_grid_cells),
        ϵ = fill(10.0, num_grid_cells),
        ne = zeros(num_grid_cells),
        K = zeros(num_grid_cells),
        nϵ = fill(1.0e19, num_grid_cells),
        cell_cache_1 = zeros(num_grid_cells),
        cell_cache_2 = zeros(num_grid_cells),
        dt_iz = [Inf],
    )
    loss_frequencies = zeros(length(fluids), num_grid_cells)

    het.apply_reactions!(fluids, rxns, cache, false, loss_frequencies)

    expected_frequency = (k1 + k2) * 1.0e18
    @test cache.dt_iz[] ≈ inv(expected_frequency)
    @test all(
        loss_frequencies[1, 2:(end - 1)] .≈ expected_frequency
    )
    @test all(iszero, loss_frequencies[2:end, :])
    expected_excitation_frequency = (k1 + k2) * 2.0e18
    expected_energy_loss = 1.0e18 * 2.0e18 * (3.0 * k1 + 5.0 * k2)
    @test all(cache.νex_explicit[2:(end - 1)] .≈ expected_excitation_frequency)
    @test all(cache.inelastic_losses[2:(end - 1)] .≈ expected_energy_loss)
end

@testset "Derived excited-state energies" begin
    # Recover an intermediate level through a reverse graph edge, demonstrating
    # that reaction files need not contain a direct ground-to-level transition.
    species = [het.Xenon(0), het.Xenon(0, 1), het.Xenon(0, 2)]
    reactions = [
        het.ElectronImpactReaction(5.0, species[1], [species[3]], zeros(256)),
        het.ElectronImpactReaction(2.0, species[2], [species[3]], zeros(256)),
    ]

    energies = het.derive_species_energies(species, reactions)

    @test energies[species[1].symbol] == 0.0
    @test energies[species[2].symbol] == 3.0
    @test energies[species[3].symbol] == 5.0

    # Independently derived paths may differ slightly because reaction-header
    # energies are rounded; differences below 0.01 eV should merge cleanly.
    rounded_reactions = [
        het.ElectronImpactReaction(10.0, species[1], [species[2]], zeros(256)),
        het.ElectronImpactReaction(10.9016, species[1], [species[3]], zeros(256)),
        het.ElectronImpactReaction(0.9015, species[2], [species[3]], zeros(256)),
    ]
    rounded_energies = het.derive_species_energies(species, rounded_reactions)
    @test rounded_energies[species[3].symbol] == 10.9016

    # Larger disagreement still indicates inconsistent reaction metadata.
    inconsistent_reactions = copy(rounded_reactions)
    inconsistent_reactions[3] =
        het.ElectronImpactReaction(0.92, species[2], [species[3]], zeros(256))
    @test_throws ErrorException het.derive_species_energies(
        species, inconsistent_reactions,
    )

    # An excited-state ionization threshold supplies the missing level energy
    # when the corresponding ground-state ionization threshold is also known.
    ionization_species = [het.Xenon(0), het.Xenon(0, 31), het.Xenon(1)]
    ionization_reactions = [
        het.ElectronImpactReaction(
            12.13, ionization_species[1], [ionization_species[3]], zeros(256),
        ),
        het.ElectronImpactReaction(
            1.2284, ionization_species[2], [ionization_species[3]], zeros(256),
        ),
    ]
    ionization_energies = het.derive_species_energies(
        ionization_species, ionization_reactions,
    )
    @test ionization_energies[ionization_species[2].symbol] ≈ 10.9016
    @test ionization_energies[ionization_species[3].symbol] == 12.13

    # A direct excitation path must agree with the energy inferred from the two
    # ionization thresholds.
    pushfirst!(
        ionization_reactions,
        het.ElectronImpactReaction(
            10.9015, ionization_species[1], [ionization_species[2]], zeros(256),
        ),
    )
    @test het.derive_species_energies(
        ionization_species, ionization_reactions,
    )[ionization_species[2].symbol] == 10.9015

    inconsistent_ionization = copy(ionization_reactions)
    inconsistent_ionization[3] = het.ElectronImpactReaction(
        1.20, ionization_species[2], [ionization_species[3]], zeros(256),
    )
    @test_throws ErrorException het.derive_species_energies(
        ionization_species, inconsistent_ionization,
    )
end
