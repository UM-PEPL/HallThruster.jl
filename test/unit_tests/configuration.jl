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

        step_trough_bohm = het.StepTroughBohm1(0.05, 1.0, 0.8, 0.5, 0.1, 0.2, 0.5)
        serialized_model = het.serialize(step_trough_bohm)
        @test serialized_model["trough_floor"] == 0.1
        test_subtype(het.AnomalousTransportModel, step_trough_bohm)

        profiled_cfg = het.Config(;
            thruster = het.SPT_100,
            discharge_voltage = 300.0,
            domain = (0.0, 0.08),
            propellants = [
                het.Propellant(
                    het.Xenon, 5.0e-6;
                    velocity_m_s = het.LinearInterpolation(
                        [0.0, 0.08], [150.0, 300.0],
                    ),
                    temperature_K = het.LinearInterpolation(
                        [0.0, 0.08], [500.0, 800.0],
                    ),
                ),
            ],
        )
        profiled_dict = het.serialize(profiled_cfg)
        @test profiled_dict["propellants"][1]["velocity_m_s"] == Dict(
            "xs" => [0.0, 0.08], "ys" => [150.0, 300.0],
        )
        test_roundtrip(het.Config, profiled_cfg)

        OD = het.Serialization.OrderedDict
        d = OD(
            "thruster" => het.serialize(het.SPT_100),
            "discharge_voltage" => 800.0,
            "domain" => (0.0, 0.4),
            "anode_mass_flow_rate" => 9.0e-6,
            "wall_loss_model" => OD(
                "type" => "NoWallLosses"
            ),
        )
        test_roundtrip(het.Config, d)
    end
    return
end

function test_configuration()
    anom_model = het.TwoZoneBohm(1 // 100, 1 // 10)
    config = het.Config(;
        ncharge = 3,
        discharge_voltage = 200,
        anode_mass_flow_rate = 5.0e-6,
        thruster = het.SPT_100,
        domain = (0.0, 5.0e-2),
        anom_model = anom_model,
        initial_condition = het.DefaultInitialization(;
            max_electron_temperature = 10.0
        )
    )

    legacy_profile = het.LinearInterpolation([0.0, 0.05], [150.0, 250.0])
    legacy_profile_config = het.Config(;
        discharge_voltage = 200.0,
        anode_mass_flow_rate = 5.0e-6,
        thruster = het.SPT_100,
        domain = (0.0, 0.05),
        neutral_velocity = legacy_profile,
    )
    legacy_config_velocity = only(legacy_profile_config.propellants).velocity_m_s
    @test legacy_config_velocity.xs == legacy_profile.xs
    @test legacy_config_velocity.ys == legacy_profile.ys

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
        @test neutral_fluid.vel_prim ==
            config.propellants[1].velocity_m_s.(params.grid.cell_centers)

        # Check array sizes
        (;
            Aϵ, bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue,
            ∇pe, νen, νei, radial_loss_frequency, νew_momentum, nn, ji,
        ) = params.cache

        for arr in (
                bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen,
                νei, radial_loss_frequency, νew_momentum, ji, nn,
            )
            @test size(arr) == (ncells + 2,)
        end

        @test length(Aϵ.dl) == ncells + 1
        @test length(Aϵ.d) == ncells + 2
        @test length(Aϵ.du) == ncells + 1
    end

    @testset "Anom initialization" begin
        # When we first initialize a simulation, the anom transport is set to this TwoZoneBohm
        # This is because we require some value of the anomalous collision frequency to
        # properly initialize the other plasma properties.
        initial_model = het.TwoZoneBohm(1 // 160, 1 / 16)
        v = anom_model(zeros(ncells + 2), params)

        # Only include interior cells here
        inds = 2:(ncells + 1)

        init = initial_model(zeros(ncells + 2), params)
        @test all(init[inds] .== params.cache.νan[inds])
        @test all(v[inds] .!= params.cache.νan[inds])

        # After one timestep, the anom. transport model we choose kicks in
        sol = het.run_from_setup(params, config)
        @test all(init[inds] .!= sol.frames[end].nu_an[inds])
        @test all(v[inds] .== sol.frames[end].nu_an[inds])
    end

    return
end

function test_multiple_propellants()

    Xe = het.Propellant(het.Xenon, flow_rate_kg_s = 4.0e-6, max_charge = 3)
    Kr = het.Propellant(het.Krypton, flow_rate_kg_s = 1.0e-6, max_charge = 2)
    ncells = 50
    simparams = het.SimParams(grid = het.UnevenGrid(ncells), duration = 1.0e-3, num_save = 1000, verbose = false)

    @testset "Setup" begin

        config = het.Config(
            propellants = [Xe, Kr],
            domain = (0.0, 0.08),
            discharge_voltage = 300.0,
            thruster = het.SPT_100,
            background_pressure_Torr = 1.0e-5,
        )

        params = het.setup_simulation(config, simparams)

        @test length(config.propellants) == 2

        expected_species = [
            het.Xenon(0);
            [het.Xenon(Z) for Z in Xe.allowed_charges]...;
            het.Krypton(0);
            [het.Krypton(Z) for Z in Kr.allowed_charges]...;
        ]

        species = [f.species for f in params.fluid_array]
        @test expected_species == species

        # Should have electron-neutral collisions for both species
        @test length(params.electron_neutral_collisions) == 2

        # Should have excitation reactions for both species
        # Excitation reactant indices should be the indices of the neutral species
        @test length(params.excitation_reactions) == 2
        @test params.excitation_reactant_indices == [findfirst(==(het.Xenon(0)), species), findfirst(==(het.Krypton(0)), species)]

        # Number of ionization reactions should be the Nth triangular number, where N is the maximum charge
        triangular(n::T) where {T <: Integer} = T(n * (n + 1) // 2)

        @test sum(length(group.channels) for group in params.reaction_groups) ==
            triangular(length(Kr.allowed_charges)) + triangular(length(Xe.allowed_charges))
        for group in params.reaction_groups
            @test all(
                channel.reaction.reactant == species[group.reactant_index]
                    for channel in group.channels
            )
            for channel in group.channels
                for (product, product_index) in zip(
                        channel.reaction.products, channel.product_indices,
                    )
                    @test product == species[product_index]
                end
            end
        end

        # Neutral ingestion densities
        ndot_xe_back, ndot_kr_back = params.ingestion_flow_rates ./ [Xe.gas.m, Kr.gas.m]
        ndot_tot_back = ndot_xe_back + ndot_kr_back
        @test ndot_xe_back / ndot_tot_back > 0.5
        @test ndot_kr_back / ndot_tot_back < 0.5
    end

    @testset "Simulation results" begin
        function eval_sims(prop1, prop2; rtol = 1.0e-6, check_kr_density = true)
            config_common = (;
                domain = (0.0, 0.08),
                discharge_voltage = 300.0,
                thruster = het.SPT_100,
                background_pressure_Torr = 1.0e-5,
                ion_wall_losses = false,
            )
            outputs = []
            for props in [prop1, prop2]
                config = het.Config(; propellants = props, config_common...)
                sim = het.run_simulation(config, simparams)
                avg = het.time_average(sim)
                push!(
                    outputs, (;
                        thrust = het.thrust(avg)[],
                        current = het.discharge_current(avg)[],
                        ion_current = het.ion_current(avg)[],
                        max_Te = maximum(avg.frames[].Tev),
                        max_E = maximum(avg.frames[].E),
                        max_ne = maximum(avg.frames[].ne),
                        max_ni_Xe = maximum(avg.frames[].ions[:Xe][1].n),
                        max_ni_Kr = if haskey(avg.frames[].ions, :Kr)
                            maximum(avg.frames[].ions[:Kr][1].n)
                        else
                            het.MIN_NUMBER_DENSITY
                        end,
                    )
                )
            end

            @test isapprox(outputs[1].thrust, outputs[2].thrust; rtol)
            @test isapprox(outputs[1].current, outputs[2].current; rtol)
            @test isapprox(outputs[1].ion_current, outputs[2].ion_current; rtol)
            @test isapprox(outputs[1].max_Te, outputs[2].max_Te; rtol)
            @test isapprox(outputs[1].max_E, outputs[2].max_E; rtol)
            @test isapprox(outputs[1].max_ne, outputs[2].max_ne; rtol)
            @test isapprox(outputs[1].max_ni_Xe, outputs[2].max_ni_Xe; rtol)
            if check_kr_density
                @test isapprox(outputs[1].max_ni_Kr, outputs[2].max_ni_Kr; rtol)
            end
        end

        # Reversing the order of propellants should not affect results
        outputs = eval_sims([Xe, Kr], [Kr, Xe])

        # Using a negligible fraction of Kr should not change the results compared to a pure xenon simulation
        Kr_small = het.Propellant(het.Krypton, flow_rate_kg_s = 1.0e-18, max_charge = 1)
        outputs = eval_sims([Xe], [Xe, Kr_small], check_kr_density = false)
    end

    return
end

function test_TOML_Read()
    @testset "TOML allowed_charges and max charge parsing" begin
        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                name = "Xenon"
                symbol = "Xe"
                flow_rate_kg_s = 5.0e-6
                allowed_charges = [-1, 1, 2]
                """
            )

            props = het.load_propellant_config(file)
            @test collect(props[1].allowed_charges) == [-1, 1, 2]
        end

        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                name = "Xenon"
                symbol = "Xe"
                flow_rate_kg_s = 5.0e-6
                max_charge = 3
                """
            )

            props = het.load_propellant_config(file)
            @test collect(props[1].allowed_charges) == [1, 2, 3]
        end

        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                name = "Xenon"
                symbol = "Xe"
                flow_rate_kg_s = 5.0e-6
                max_charge = 3
                allowed_charges = [-1, 1, 2]
                """
            )
            @test_throws Exception het.load_propellant_config(file)
        end

        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                name = "Xenon"
                symbol = "Xe"
                flow_rate_kg_s = 5.0e-6
                """
            )
            props = het.load_propellant_config(file)
            @test collect(props[1].allowed_charges) == [1]
        end

        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                symbol = "Xe"
                flow_rate_kg_s = 5.0e-6
                velocity_m_s = { xs = [0.0, 0.08], ys = [150.0, 300.0] }
                temperature_K = { xs = [0.0, 0.08], ys = [500.0, 800.0] }
                """
            )
            propellant = only(het.load_propellant_config(file))
            @test propellant.velocity_m_s.xs == [0.0, 0.08]
            @test propellant.velocity_m_s.ys == [150.0, 300.0]
            @test propellant.temperature_K.ys == [500.0, 800.0]
        end
    end

    return
end

test_TOML_Read()

function test_propellant_config_overrides()
    @testset "Propellant chemistry inheritance" begin
        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                symbol = "Xe"
                max_charge = 3
                excited_levels = [1, 2]
                excited_ion_levels = { 1 = [1], 2 = [1, 2], 3 = [1] }
                """
            )

            config_args = (;
                thruster = het.SPT_100,
                domain = (0.0, 0.08),
                discharge_voltage = 300.0,
                propellant_config = file,
            )

            # Physical overrides retain chemistry options omitted from Propellant.
            config = het.Config(;
                config_args...,
                propellants = [
                    het.Propellant("Xe"; flow_rate_kg_s = 5.0e-6, velocity_m_s = 150.0),
                ],
            )
            propellant = only(config.propellants)
            @test propellant.flow_rate_kg_s == 5.0e-6
            @test propellant.velocity_m_s(0.0) == 150.0
            @test propellant.allowed_charges == [1, 2, 3]
            @test propellant.excited_levels == [1, 2]
            @test propellant.excited_ion_levels == Dict(1 => [1], 2 => [1, 2], 3 => [1])

            # Specifying only max_charge still inherits both excited-state collections.
            config = het.Config(;
                config_args...,
                propellants = [
                    het.Propellant("Xe"; flow_rate_kg_s = 5.0e-6, max_charge = 3),
                ],
            )
            propellant = only(config.propellants)
            @test propellant.excited_levels == [1, 2]
            @test propellant.excited_ion_levels == Dict(1 => [1], 2 => [1, 2], 3 => [1])

            # Non-default chemistry values override the file; inherited ion levels
            # are restricted to the explicitly allowed charge states.
            config = het.Config(;
                config_args...,
                propellants = [
                    het.Propellant(
                        "Xe"; flow_rate_kg_s = 5.0e-6, max_charge = 2,
                        excited_levels = [4],
                    ),
                ],
            )
            propellant = only(config.propellants)
            @test propellant.allowed_charges == [1, 2]
            @test propellant.excited_levels == [4]
            @test propellant.excited_ion_levels == Dict(1 => [1], 2 => [1, 2])
        end


        mktempdir() do dir
            file = joinpath(dir, "propellant.toml")
            write(
                file, """
                [[species]]
                symbol = "Xe"
                max_charge = 1

                [[species]]
                symbol = "Kr"
                max_charge = 1
                """
            )
            xenon_velocity = het.LinearInterpolation(
                [0.0, 0.08], [150.0, 300.0],
            )
            config = het.Config(;
                thruster = het.SPT_100,
                domain = (0.0, 0.08),
                discharge_voltage = 300.0,
                propellant_config = file,
                propellants = [
                    het.Propellant(
                        het.Xenon, 5.0e-6; velocity_m_s = xenon_velocity,
                    ),
                ],
            )
            krypton = only(filter(p -> p.gas == het.Krypton, config.propellants))
            scale = sqrt(het.Xenon.M / het.Krypton.M)
            @test krypton.flow_rate_kg_s == 0.0
            @test krypton.velocity_m_s.xs == xenon_velocity.xs
            @test krypton.velocity_m_s.ys ≈ scale .* xenon_velocity.ys
        end
    end
    return
end

test_propellant_config_overrides()

function test_allowed_charges_initialization()
    @testset "Allowed charges initialization" begin
        Xe_default = het.Propellant(
            het.Xenon,
            flow_rate_kg_s = 5.0e-6,
        )
        @test collect(Xe_default.allowed_charges) == [1]

        Xe = het.Propellant(
            het.Xenon,
            flow_rate_kg_s = 5.0e-6,
            max_charge = 3,
        )
        @test collect(Xe.allowed_charges) == [1, 2, 3]

        Xe_neg = het.Propellant(
            het.Xenon,
            flow_rate_kg_s = 5.0e-6,
            allowed_charges = [-1, 1, 2, 3],
        )
        @test collect(Xe_neg.allowed_charges) == [-1, 1, 2, 3]
    end
    return
end

test_allowed_charges_initialization()
test_config_serialization()
test_configuration()
@testset "Multiple propellants" begin
    test_multiple_propellants()
end
