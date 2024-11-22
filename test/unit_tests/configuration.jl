function test_configuration()
    @testset "Configuration" begin
        common_opts = (;
            ncharge = 3,
            discharge_voltage = 300u"V",
            anode_mass_flow_rate = 5u"mg/s",
            thruster = het.SPT_100,
            domain = (0.0u"cm", 5.0u"cm"),
        )

        config = het.Config(;
            background_pressure = 0.0u"Torr",
            background_neutral_temperature = 0.0u"K",
            common_opts...,
        )

        fluids, fluid_ranges, species, species_range_dict, is_velocity_index = het.configure_fluids(config)

        @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]
        @test species == [Xenon(0), Xenon(1), Xenon(2), Xenon(3)]
        @test species_range_dict == Dict(
            Symbol("Xe") => 1:1,
            Symbol("Xe+") => 2:3,
            Symbol("Xe2+") => 4:5,
            Symbol("Xe3+") => 6:7,
        )

        @test fluids[1] == het.ContinuityOnly(
            species[1], config.neutral_velocity, config.neutral_temperature,)
        @test fluids[2] == het.IsothermalEuler(species[2], config.ion_temperature)
        @test fluids[3] == het.IsothermalEuler(species[3], config.ion_temperature)
        @test fluids[4] == het.IsothermalEuler(species[4], config.ion_temperature)
        @test is_velocity_index == [false, false, true, false, true, false, true]

        index = het.configure_index(fluids, fluid_ranges)
        @test keys(index) == (:ρn, :ρi, :ρiui)
        @test values(index) == (1, [2, 4, 6], [3, 5, 7])

        # load collisions and reactions
        ionization_reactions = het._load_reactions(
            config.ionization_model, unique(species),)
        ionization_reactant_indices = het.reactant_indices(
            ionization_reactions, species_range_dict,)
        @test ionization_reactant_indices == [1, 1, 1, 2, 2, 4]

        ionization_product_indices = het.product_indices(
            ionization_reactions, species_range_dict,)
        @test ionization_product_indices == [2, 4, 6, 4, 6, 6]

        excitation_reactions = het._load_reactions(
            config.excitation_model, unique(species),)
        excitation_reactant_indices = het.reactant_indices(
            excitation_reactions, species_range_dict,)
        @test excitation_reactant_indices == [1]

        # Test that initialization and configuration works properly when background neutrals are included

        pB = 5e-6u"Torr"
        TB = 120u"K"

        config_bg = het.Config(;
            background_pressure = pB,
            background_neutral_temperature = TB,
            common_opts...,
        )

        fluids, fluid_ranges, species, species_range_dict = het.configure_fluids(config_bg)
        @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]
        @test species == [Xenon(0), Xenon(1), Xenon(2), Xenon(3)]
        @test species_range_dict == Dict(
            Symbol("Xe") => 1:1,
            Symbol("Xe+") => 2:3,
            Symbol("Xe2+") => 4:5,
            Symbol("Xe3+") => 6:7,
        )

        @test fluids[1] == het.ContinuityOnly(
            species[1], config.neutral_velocity, config.neutral_temperature,)
        @test fluids[2] == het.IsothermalEuler(species[2], config.ion_temperature)
        @test fluids[3] == het.IsothermalEuler(species[3], config.ion_temperature)
        @test fluids[4] == het.IsothermalEuler(species[4], config.ion_temperature)
        @test is_velocity_index == [false, false, true, false, true, false, true]

        index = het.configure_index(fluids, fluid_ranges)
        @test keys(index) == (:ρn, :ρi, :ρiui)
        @test values(index) == (1, [2, 4, 6], [3, 5, 7])

        # load collisions and reactions
        ionization_reactions = het._load_reactions(
            config.ionization_model, unique(species),)
        ionization_reactant_indices = het.reactant_indices(
            ionization_reactions, species_range_dict,)
        @test ionization_reactant_indices == [1, 1, 1, 2, 2, 4]

        ionization_product_indices = het.product_indices(
            ionization_reactions, species_range_dict,)
        @test ionization_product_indices == [2, 4, 6, 4, 6, 6]

        excitation_reactions = het._load_reactions(
            config.excitation_model, unique(species),)
        excitation_reactant_indices = het.reactant_indices(
            excitation_reactions, species_range_dict,)
        @test excitation_reactant_indices == [1]
    end
end
