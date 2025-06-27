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
    anom_model = het.NoAnom()
    config = het.Config(;
        ncharge = 3,
        discharge_voltage = 300,
        anode_mass_flow_rate = 5.0e-6,
        thruster = het.SPT_100,
        domain = (0.0, 5.0e-2),
        anom_model = anom_model,
        initial_condition = het.DefaultInitialization(;
            max_electron_temperature = 10.0
        )
    )

    ncells = 100

    simparams = het.SimParams(
        grid = het.EvenGrid(ncells),
        duration = 1.0e-6,
    )

    params = het.setup_simulation(config, simparams)

    @testset "Configuration" begin
        species = Set([f.species for f in params.fluid_array])
        expected_species = Set([het.Xenon(0), het.Xenon(1), het.Xenon(2), het.Xenon(3)])
        @test species == expected_species

        neutral_fluid = params.fluid_containers.continuity[1]
        @test length(params.fluid_containers.continuity) == 1
        @test neutral_fluid.const_velocity == config.propellants[1].velocity_m_s
        @test het.temperature(neutral_fluid) ≈ config.propellants[1].temperature_K

        for ion_fluid in params.fluid_containers.isothermal
            @test het.temperature(ion_fluid) ≈ config.propellants[1].ion_temperature_K
        end

        # Check array sizes
        (;
            Aϵ, bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue,
            ∇pe, νen, νei, radial_loss_frequency, νew_momentum, ni, ui, niui, nn, ji,
        ) = params.cache

        for arr in (
                bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen,
                νei, radial_loss_frequency, νew_momentum, ji, nn,
            )
            @test size(arr) == (ncells + 2,)
        end

        for arr in (ni, ui, niui)
            @test size(arr) == (3, ncells + 2)
        end

        @test length(Aϵ.dl) == ncells + 1
        @test length(Aϵ.d) == ncells + 2
        @test length(Aϵ.du) == ncells + 1
    end

    return @testset "Anom initialization" begin
        initial_model = het.TwoZoneBohm(1 // 160, 1 / 16)
        v = anom_model(zeros(ncells + 2), params, config)
        init = initial_model(zeros(ncells + 2), params, config)

        @test all(init .== params.cache.νan)
        @test all(v .!= params.cache.νan)

        sol = het.run_from_setup(params, config)
        @test all(init .!= sol.params.cache.νan)
        @test all(v .== sol.params.cache.νan)
    end
end

test_config_serialization()
test_configuration()
