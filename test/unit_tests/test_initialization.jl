@testset "Initialization" begin

    domain = (0.0, 0.08)
    thruster = HallThruster.SPT_100
    config = HallThruster.Config(;
        anom_model = HallThruster.NoAnom(),
        ncharge = 3,
        propellant = HallThruster.Xenon,
        thruster,
        domain,
        cathode_Te = 5.0,
        anode_Te = 3.0,
        neutral_velocity = 300.0,
        neutral_temperature = 100.0,
        ion_temperature = 300.0,
        initial_condition = HallThruster.DefaultInitialization(),
        anode_mass_flow_rate = 3e-6,
        discharge_voltage = 500.0
    )

    mi = config.propellant.m

    ncells = 100
    fluids, fluid_ranges, species, species_range_dict = HallThruster.configure_fluids(config)
    grid = HallThruster.generate_grid(config.thruster.geometry, ncells, domain)
    U, cache = HallThruster.allocate_arrays(grid, fluids)
    index = HallThruster.configure_index(fluids, fluid_ranges)

    params = (;
        index,
        cache,
        z_cell = grid.cell_centers,
        config,
    )

    HallThruster.initialize!(U, params)

    @test abs((U[index.ρn[1], 1] / U[index.ρn[1],end]) - 100) < 1

    ne = [HallThruster.electron_density(U, params, i) for i in 1:ncells+2]
    ϵ = U[index.nϵ, :] ./ ne

    ui = [
        U[index.ρiui[Z], :] ./ U[index.ρi[Z], :] for Z in 1:config.ncharge
    ]

    @test ui[1][1] ≈ -sqrt(HallThruster.e * config.anode_Te / mi)
    @test ui[2][1] ≈ -sqrt(2 * HallThruster.e * config.anode_Te / mi)
    @test ui[3][1] ≈ -sqrt(3 * HallThruster.e * config.anode_Te / mi)

    @test abs(2/3 * ϵ[1] - config.anode_Te) < 0.1
    @test abs(2/3 * ϵ[end] - config.cathode_Te) < 0.1

    struct NewInitialization <: HallThruster.InitialCondition end
    @test_throws ArgumentError HallThruster.initialize!(U, params, NewInitialization())
end

@testset "Configuration" begin

    common_opts = (;
        ncharge=3,
        discharge_voltage = 300u"V",
        anode_mass_flow_rate = 5u"mg/s",
        thruster = HallThruster.SPT_100,
        domain = (0.0u"cm", 5.0u"cm"),
    )

    config = HallThruster.Config(;
        solve_background_neutrals = false,
        background_pressure = 0.0u"Torr",
        background_neutral_temperature = 0.0u"K",
        common_opts...
    )

    fluids, fluid_ranges, species, species_range_dict = HallThruster.configure_fluids(config)

    @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]
    @test species == [Xenon(0), Xenon(1), Xenon(2), Xenon(3)]
    @test species_range_dict == Dict(
        Symbol("Xe") => [1:1],
        Symbol("Xe+") => [2:3],
        Symbol("Xe2+") => [4:5],
        Symbol("Xe3+") => [6:7],
    )

    @test fluids[1].conservation_laws == HallThruster.ContinuityOnly(config.neutral_velocity, config.neutral_temperature)
    @test fluids[2].conservation_laws == HallThruster.IsothermalEuler(config.ion_temperature)
    @test fluids[3].conservation_laws == HallThruster.IsothermalEuler(config.ion_temperature)
    @test fluids[4].conservation_laws == HallThruster.IsothermalEuler(config.ion_temperature)


    index = HallThruster.configure_index(fluids, fluid_ranges)
    @test keys(index) == (:ρn, :ρi, :ρiui, :nϵ, :lf)
    @test values(index) == ([1], [2, 4, 6], [3, 5, 7], 8, 7)

    # load collisions and reactions
    ionization_reactions = HallThruster._load_reactions(config.ionization_model, unique(species))
    ionization_reactant_indices = HallThruster.reactant_indices(ionization_reactions, species_range_dict)
    @test ionization_reactant_indices == [[1], [1], [1], [2], [2], [4]]

    ionization_product_indices = HallThruster.product_indices(ionization_reactions, species_range_dict)
    @test ionization_product_indices == [[2], [4], [6], [4], [6], [6]]

    excitation_reactions = HallThruster._load_reactions(config.excitation_model, unique(species))
    excitation_reactant_indices = HallThruster.reactant_indices(excitation_reactions, species_range_dict)
    @test excitation_reactant_indices == [[1]]

    # Test that initialization and configuration works properly when background neutrals are included

    pB = 5e-6u"Torr"
    TB = 120u"K"

    config_bg = HallThruster.Config(;
        solve_background_neutrals = true,
        background_pressure = pB,
        background_neutral_temperature = TB,
        common_opts...
    )

    fluids, fluid_ranges, species, species_range_dict = HallThruster.configure_fluids(config_bg)
    @test fluid_ranges == [1:1, 2:2, 3:4, 5:6, 7:8]
    @test species == [Xenon(0), Xenon(0), Xenon(1), Xenon(2), Xenon(3)]
    @test species_range_dict == Dict(
        Symbol("Xe") => [1:1, 2:2],
        Symbol("Xe+") => [3:4],
        Symbol("Xe2+") => [5:6],
        Symbol("Xe3+") => [7:8],
    )

    @test fluids[1].conservation_laws == HallThruster.ContinuityOnly(config.neutral_velocity, config.neutral_temperature)
    @test fluids[2].conservation_laws == HallThruster.ContinuityOnly(-sqrt(ustrip(TB) * HallThruster.kB / Xenon.m), ustrip(TB))
    @test fluids[3].conservation_laws == HallThruster.IsothermalEuler(config.ion_temperature)
    @test fluids[4].conservation_laws == HallThruster.IsothermalEuler(config.ion_temperature)
    @test fluids[5].conservation_laws == HallThruster.IsothermalEuler(config.ion_temperature)

    index = HallThruster.configure_index(fluids, fluid_ranges)
    @test keys(index) == (:ρn, :ρi, :ρiui, :nϵ, :lf)
    @test values(index) == ([1, 2], [3, 5, 7], [4, 6, 8], 9, 8)

    # load collisions and reactions
    ionization_reactions = HallThruster._load_reactions(config.ionization_model, unique(species))
    ionization_reactant_indices = HallThruster.reactant_indices(ionization_reactions, species_range_dict)
    @test ionization_reactant_indices == [[1, 2], [1, 2], [1, 2], [3], [3], [5]]

    ionization_product_indices = HallThruster.product_indices(ionization_reactions, species_range_dict)
    @test ionization_product_indices == [[3], [5], [7], [5], [7], [7]]

    excitation_reactions = HallThruster._load_reactions(config.excitation_model, unique(species))
    excitation_reactant_indices = HallThruster.reactant_indices(excitation_reactions, species_range_dict)
    @test excitation_reactant_indices == [[1, 2]]
end
