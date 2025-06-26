using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

function test_config_serialization()
    @testset "Serialization" begin
        cfg = het.Config(;
            thruster = het.SPT_100,
            discharge_voltage = 300.0,
            domain = (0.0, 0.8),
            anode_mass_flow_rate = 5.0e-6,
        )

        d = het.serialize(cfg)
        for k in het.Serialization.exclude(het.Config)
            @test !haskey(d, string(k))
        end

        test_roundtrip(het.Config, cfg)

        OD = het.Serialization.OrderedDict
        d = OD(
            :thruster => het.serialize(het.SPT_100),
            :discharge_voltage => 800.0,
            :domain => (0.0, 0.4),
            :anode_mass_flow_rate => 9.0e-6,
            :wall_loss_model => OD(
                :type => "NoWallLosses"
            ),
        )
        test_roundtrip(het.Config, d)
    end
    return
end

function test_configuration()
    @testset "Configuration" begin
        common_opts = (;
            ncharge = 3,
            discharge_voltage = 300,
            anode_mass_flow_rate = 5.0e-6,
            thruster = het.SPT_100,
            domain = (0.0, 5.0e-2),
        )

        config = het.Config(;
            background_pressure_Torr = 0.0,
            background_temperature_K = 0.0,
            common_opts...,
        )

        fluids, fluid_ranges, species, species_range_dict, is_velocity_index = het.configure_fluids(config)

        @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]
        @test species == [het.Xenon(0), het.Xenon(1), het.Xenon(2), het.Xenon(3)]
        @test species_range_dict == Dict(
            Symbol("Xe") => 1:1,
            Symbol("Xe+") => 2:3,
            Symbol("Xe2+") => 4:5,
            Symbol("Xe3+") => 6:7,
        )

        prop = config.propellants[1]
        @test fluids[1] == het.ContinuityOnly(species[1], prop.velocity_m_s, prop.temperature_K)
        @test fluids[2] == het.IsothermalEuler(species[2], prop.ion_temperature_K)
        @test fluids[3] == het.IsothermalEuler(species[3], prop.ion_temperature_K)
        @test fluids[4] == het.IsothermalEuler(species[4], prop.ion_temperature_K)
        @test is_velocity_index == [false, false, true, false, true, false, true]

        index = het.configure_index(fluids, fluid_ranges)
        @test keys(index) == (:ρn, :ρi, :ρiui)
        @test values(index) == (1, [2, 4, 6], [3, 5, 7])

        # load collisions and reactions
        ionization_reactions = het.load_ionization_reactions(
            config.ionization_model, unique(species),
        )
        ionization_reactant_indices = het.reactant_indices(
            ionization_reactions, species_range_dict,
        )
        @test ionization_reactant_indices == [1, 1, 1, 2, 2, 4]

        ionization_product_indices = het.product_indices(
            ionization_reactions, species_range_dict,
        )
        @test ionization_product_indices == [2, 4, 6, 4, 6, 6]

        excitation_reactions = het.load_excitation_reactions(
            config.excitation_model, unique(species),
        )
        excitation_reactant_indices = het.reactant_indices(
            excitation_reactions, species_range_dict,
        )
        @test excitation_reactant_indices == [1]

        # Test that initialization and configuration works properly when background neutrals are included
        pB_Torr = 5.0e-6
        TB_K = 120.0

        config_bg = het.Config(;
            background_pressure_Torr = pB_Torr,
            background_temperature_K = TB_K,
            common_opts...,
        )

        fluids, fluid_ranges, species, species_range_dict = het.configure_fluids(config_bg)
        @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]
        @test species == [het.Xenon(0), het.Xenon(1), het.Xenon(2), het.Xenon(3)]
        @test species_range_dict == Dict(
            Symbol("Xe") => 1:1,
            Symbol("Xe+") => 2:3,
            Symbol("Xe2+") => 4:5,
            Symbol("Xe3+") => 6:7,
        )

        @test fluids[1] == het.ContinuityOnly(species[1], prop.velocity_m_s, prop.temperature_K)
        @test fluids[2] == het.IsothermalEuler(species[2], prop.ion_temperature_K)
        @test fluids[3] == het.IsothermalEuler(species[3], prop.ion_temperature_K)
        @test fluids[4] == het.IsothermalEuler(species[4], prop.ion_temperature_K)
        @test is_velocity_index == [false, false, true, false, true, false, true]

        index = het.configure_index(fluids, fluid_ranges)
        @test keys(index) == (:ρn, :ρi, :ρiui)
        @test values(index) == (1, [2, 4, 6], [3, 5, 7])

        # load collisions and reactions
        ionization_reactions = het.load_ionization_reactions(
            config.ionization_model, unique(species),
        )
        ionization_reactant_indices = het.reactant_indices(
            ionization_reactions, species_range_dict,
        )
        @test ionization_reactant_indices == [1, 1, 1, 2, 2, 4]

        ionization_product_indices = het.product_indices(
            ionization_reactions, species_range_dict,
        )
        @test ionization_product_indices == [2, 4, 6, 4, 6, 6]

        excitation_reactions = het.load_excitation_reactions(
            config.excitation_model, unique(species),
        )
        excitation_reactant_indices = het.reactant_indices(
            excitation_reactions, species_range_dict,
        )
        @test excitation_reactant_indices == [1]
    end
    return
end

test_config_serialization()
test_configuration()
