using HallThruster: HallThruster as het
using Test

@testset "Smoothing stencil and Neumann boundaries" begin
    profile = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    # Explicit stencil values, with ghost cells equal to the filtered neighbors.
    expected = [1.0, 1.0, 2.48, 5.02, 10.04, 18.0, 18.0]
    actual = copy(profile)
    @test het.smooth!(actual, similar(actual); iters = 1) === actual
    @test actual ≈ expected

    # A second pass must filter the result of the first, including its new ghosts.
    expected_twice = [1.37, 1.37, 3.0992, 6.1528, 10.7696, 16.01, 16.01]
    filtered_twice = het.smooth!(copy(profile), similar(profile); iters = 2)
    @test filtered_twice ≈ expected_twice

    unchanged = copy(profile)
    scratch = fill(-1.0, length(profile))
    @test het.smooth!(unchanged, scratch; iters = 0) === unchanged
    @test unchanged == profile
    @test all(==(-1.0), scratch)

    for iters in (1, 3)
        filtered = het.smooth!(copy(profile), similar(profile); iters)
        @test filtered[1] == filtered[2]
        @test filtered[end] == filtered[end - 1]
        reflected = het.smooth!(reverse(profile), similar(profile); iters)
        @test filtered ≈ reverse(reflected)

        constant_profile = fill(7.0, length(profile))
        het.smooth!(constant_profile, similar(profile); iters)
        @test constant_profile ≈ fill(7.0, length(profile))
    end
end

@testset "Bohm smoothing does not compound across solver steps" begin
    for iters in (0, 1, 3), multiplier in (1.0, 2.0)
        @testset "iterations = $iters, multiplier = $multiplier" begin
            config = het.Config(
                domain = (0.0, 0.08),
                discharge_voltage = 300.0,
                thruster = het.SPT_100,
                anom_model = het.TwoZoneBohm(1 / 160, 1 / 16),
                propellants = [het.Propellant(het.Xenon, flow_rate_kg_s = 5e-6)],
                anom_smoothing_iters = iters,
            )
            sim = het.SimParams(
                duration = 1e-8,
                dt = 1e-10,
                adaptive = false,
                grid = het.EvenGrid(20),
                num_save = 2,
                verbose = false,
                current_control = het.NoController(),
            )
            params = het.setup_simulation(config, sim)
            base = config.anom_model(zeros(length(params.grid.cell_centers)), params)
            smoothed = het.smooth!(copy(base), similar(base); iters)
            params.cache.anom_multiplier[] = multiplier
            solver = het.SolverState(params; tspan = (0.0, sim.duration))

            # Cross iteration 5, when static Bohm models stop refreshing the base.
            for _ in 1:50
                @test het.step_solver!(params, config, solver)
                @test params.cache.anom_multiplier[] ≈ multiplier
                @test params.cache.νan_base ≈ base
                @test params.cache.νan ≈ multiplier .* smoothed
            end
            @test params.iteration[] == 51
            @test solver.retcode == :success
        end
    end
end
