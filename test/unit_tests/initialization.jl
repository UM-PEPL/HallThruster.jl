using HallThruster: HallThruster as het

struct TestAnomModel <: het.AnomalousTransportModel end

het.num_anom_variables(::TestAnomModel) = 3

struct NewInitialization <: het.InitialCondition end

function test_sim_initialization()
    return @testset "Sim initialization" begin
        domain = (0.0, 0.08)
        thruster = het.SPT_100
        config = het.Config(;
            anom_model = TestAnomModel(),
            ncharge = 3,
            propellant = het.Xenon,
            thruster,
            domain,
            cathode_Tev = 5.0,
            anode_Tev = 3.0,
            neutral_velocity = 300.0,
            neutral_temperature_K = 100.0,
            ion_temperature_K = 300.0,
            initial_condition = het.DefaultInitialization(;
                max_electron_temperature = 10.0
            ),
            anode_mass_flow_rate = 3.0e-6,
            discharge_voltage = 500.0,
        )

        # TODO: mutliple props + fluid containers
        mi = config.propellants[1].gas.m

        ncells = 100
        fluids, fluid_ranges, _, _, _ = het.configure_fluids(config)
        grid = het.generate_grid(het.EvenGrid(ncells), config.thruster.geometry, domain)
        U, cache = het.allocate_arrays(grid, config)
        index = het.configure_index(fluids, fluid_ranges)

        params = (;
            het.params_from_config(config)...,
            index,
            cache,
            grid,
            min_Te = min(config.anode_Tev, config.cathode_Tev),
        )

        het.initialize!(U, params, config)

        @test abs((U[index.ρn, 1] / U[index.ρn, end]) - 100) < 1

        ne = [het.electron_density(U, params, i) for i in 1:(ncells + 2)]
        ϵ = params.cache.nϵ ./ ne

        max_Te = 2 / 3 * maximum(ϵ)
        @test 9 <= max_Te <= 11

        ui = [U[index.ρiui[Z], :] ./ U[index.ρi[Z], :] for Z in 1:(config.propellants[1].max_charge)]

        @test ui[1][1] ≈ -sqrt(het.e * config.anode_Tev / mi)
        @test ui[2][1] ≈ -sqrt(2 * het.e * config.anode_Tev / mi)
        @test ui[3][1] ≈ -sqrt(3 * het.e * config.anode_Tev / mi)

        @test abs(2 / 3 * ϵ[1] - config.anode_Tev) < 0.1
        @test abs(2 / 3 * ϵ[end] - config.cathode_Tev) < 0.1

        @test cache.anom_variables == [zeros(102) for _ in 1:3]

        @test_throws ArgumentError het.initialize!(U, params, config, NewInitialization())
    end
end

function test_anom_initialization()
    return @testset "Anom initialization" begin
        anom_model = het.NoAnom()

        config = het.Config(;
            thruster = het.SPT_100,
            domain = (0.0, 0.08),
            discharge_voltage = 300.0,
            anode_mass_flow_rate = 5.0e-6,
            wall_loss_model = het.WallSheath(het.BoronNitride, 0.15),
            anom_model,
        )

        ncells = 10
        sim_options = (; ncells, nsave = 2, verbose = false)
        initial_model = het.TwoZoneBohm(1 // 160, 1 / 16)

        # Check that anomalous transport is initialized to a two-zone Bohm approximation instead of the prescribed NoAnom.
        sol = het.run_simulation(
            config; ncells, dt = 0.0, duration = 0.0, sim_options...,
        )
        @test initial_model(zeros(ncells + 2), sol.params, config) == sol.params.cache.νan
        @test anom_model(zeros(ncells + 2), sol.params, config) != sol.params.cache.νan

        # Check that after one iteration, the anomalous transport is the correct value
        dt = 1.0e-8
        sol = het.run_simulation(config; ncells, dt, duration = dt, sim_options...)
        @test initial_model(zeros(ncells + 2), sol.params, config) != sol.params.cache.νan
        @test anom_model(zeros(ncells + 2), sol.params, config) == sol.params.cache.νan
    end
end

test_sim_initialization()
test_anom_initialization()
