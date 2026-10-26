using HallThruster: HallThruster as het
using Test

@testset "Every requested save time is respected" begin
    cases = (
        (; name = "first step spans multiple saves", duration = 1e-8, dt = 5e-9, num_save = 6, adaptive = false),
        (; name = "timestep exceeds duration", duration = 1e-9, dt = 2e-9, num_save = 6, adaptive = false),
        (; name = "initial and final frames only", duration = 1e-9, dt = 2e-9, num_save = 2, adaptive = false),
        (; name = "fixed timestep is restored", duration = 1e-9, dt = 3e-10, num_save = 3, adaptive = false),
        (; name = "nonuniform floating-point save intervals", duration = 1.1e-9, dt = 3e-10, num_save = 8, adaptive = false),
        (; name = "adaptive steps below min_dt", duration = 1e-12, dt = 1e-9, num_save = 6, adaptive = true),
        (; name = "adaptive initial step is clipped", duration = 1e-14, dt = 1e-9, num_save = 6, adaptive = true),
    )

    for case in cases
        @testset "$(case.name)" begin
            step_sizes = Float64[]
            last_iteration = Ref(0)
            source = function (_, params)
                # The heavy-species solver evaluates the source more than once
                # per step. Record the timestep actually used once per iteration.
                if params.iteration[] != last_iteration[]
                    push!(step_sizes, params.dt[])
                    last_iteration[] = params.iteration[]
                end
                return nothing
            end
            config = het.Config(
                domain = (0.0, 0.08),
                discharge_voltage = 300.0,
                thruster = het.SPT_100,
                anode_mass_flow_rate = 5e-6,
                source_heavy_species = source,
            )
            sim = het.SimParams(
                grid = het.EvenGrid(20),
                duration = case.duration,
                dt = case.dt,
                num_save = case.num_save,
                adaptive = case.adaptive,
                verbose = false,
            )
            sol = het.run_simulation(config, sim)
            expected_times = collect(range(0.0, case.duration; length = case.num_save))

            @test sol.retcode == :success
            @test length(sol.frames) == case.num_save
            @test length(sol.t) == length(sol.frames)
            @test sol.t == expected_times
            @test all(>(0.0), diff(sol.t))
            @test !isempty(step_sizes)
            @test all(dt -> isfinite(dt) && dt > 0, step_sizes)
            @test sum(step_sizes) ≈ case.duration
            @test maximum(step_sizes) <= maximum(diff(expected_times))
            @test sol.frames[1].dt[] == 0.0
            @test all(f -> 0 < f.dt[] <= maximum(diff(expected_times)), sol.frames[2:end])

            if case.name == "fixed timestep is restored"
                @test step_sizes ≈ [3e-10, 2e-10, 3e-10, 2e-10]
                @test [f.dt[] for f in sol.frames[2:end]] ≈ [2e-10, 2e-10]
            elseif case.adaptive
                @test all(<(sim.min_dt), step_sizes)
            end
        end
    end
end

@testset "Solver save-count override and completed steps" begin
    config = het.Config(
        domain = (0.0, 0.08),
        discharge_voltage = 300.0,
        thruster = het.SPT_100,
        anode_mass_flow_rate = 5e-6,
    )
    sim = het.SimParams(
        grid = het.EvenGrid(20), duration = 1e-9, dt = 5e-10,
        num_save = 2, adaptive = false, verbose = false,
    )
    params = het.setup_simulation(config, sim)
    tspan = (2e-9, 3e-9)
    solver = het.SolverState(params; tspan, num_save = 6)
    for i in 2:6
        previous_time = solver.t
        @test het.step_solver!(params, config, solver) == (i < 6)
        @test solver.t == solver.saveat[i]
        @test length(solver.frames) == i
        @test params.dt[] == solver.saveat[i] - previous_time
        @test solver.frames[end].dt[] == params.dt[]
    end
    @test solver.retcode == :success
    @test solver.times == collect(range(tspan...; length = 6))

    iteration = params.iteration[]
    @test !het.step_solver!(params, config, solver)
    @test params.iteration[] == iteration
    @test length(solver.frames) == 6

    for num_save in (0, 1, -2)
        @test_throws ArgumentError het.SolverState(params; tspan, num_save)
    end
    for invalid_tspan in ((0.0, 0.0), (1.0, 0.0), (0.0, Inf), (NaN, 1.0))
        @test_throws ArgumentError het.SolverState(params; tspan = invalid_tspan)
    end
    @test_throws ArgumentError het.SolverState(params; tspan = (1.0, nextfloat(1.0)), num_save = 3)
    params.simulation.num_save = 1
    @test_throws ArgumentError het.SolverState(params; tspan)
end

@testset "Timestep selection keeps the actual step separate from its limit" begin
    config = het.Config(
        domain = (0.0, 0.08), discharge_voltage = 300.0,
        thruster = het.SPT_100, anode_mass_flow_rate = 5e-6,
    )
    sim = het.SimParams(
        grid = het.EvenGrid(20), duration = 1e-8, dt = 3e-10,
        num_save = 2, adaptive = true, verbose = false,
        min_dt = 1e-10, max_dt = 1e-9,
    )
    params = het.setup_simulation(config, sim)
    solver = het.SolverState(params; tspan = (0.0, sim.duration))

    params.cache.dt_limit[] = 8e-10
    selection = het.select_timestep!(params, solver)
    @test params.dt[] == het.initial_timestep(sim)
    @test !selection.at_minimum
    @test params.cache.dt_limit[] == 8e-10

    # These selections do not advance the simulation or mutate the stability limit.
    params.iteration[] = 2
    for (limit, expected) in ((1e-12, sim.min_dt), (4e-10, 4e-10), (Inf, sim.max_dt))
        params.cache.dt_limit[] = limit
        selection = het.select_timestep!(params, solver)
        @test params.dt[] == expected
        @test selection.at_minimum == (expected == sim.min_dt)
        @test params.cache.dt_limit[] == limit
        @test solver.t == 0.0
        @test length(solver.times) == 1
    end

    # Fixed stepping ignores the previous actual step and the adaptive limit.
    params.simulation.adaptive = false
    params.dt[] = 1e-15
    het.select_timestep!(params, solver)
    @test params.dt[] == sim.dt

    # Frame output records a copy of the actual timestep, not the stability limit.
    frame = het.Frame(params, sim.dt)
    @test frame.dt[] == sim.dt
    params.dt[] = 9e-10
    @test frame.dt[] == sim.dt
end

@testset "Adaptive fallback expires even when dt equals min_dt" begin
    config = het.Config(
        domain = (0.0, 0.08), discharge_voltage = 300.0,
        thruster = het.SPT_100, anode_mass_flow_rate = 5e-6,
    )
    sim = het.SimParams(
        grid = het.EvenGrid(20), duration = 1e-8, dt = 1e-10,
        num_save = 2, adaptive = true, verbose = false,
        min_dt = 1e-10, max_dt = 1e-9, max_small_steps = 2,
    )
    params = het.setup_simulation(config, sim)
    solver = het.SolverState(params; tspan = (0.0, sim.duration))
    @test het.step_solver!(params, config, solver)
    @test solver.small_steps == 0
    @test solver.fallback_steps_remaining == 0

    # Two minimum steps trigger two fallback steps, then adaptation resumes.
    limits = [sim.min_dt, sim.min_dt, 3e-10, 3e-10, 3e-10]
    expected_steps = [sim.min_dt, sim.min_dt, sim.dt, sim.dt, 3e-10]
    expected_small_steps = [1, 0, 0, 0, 0]
    expected_remaining = [0, 2, 1, 0, 0]
    for i in eachindex(limits)
        params.cache.dt_limit[] = limits[i]
        @test het.step_solver!(params, config, solver)
        @test params.dt[] == expected_steps[i]
        @test solver.small_steps == expected_small_steps[i]
        @test solver.fallback_steps_remaining == expected_remaining[i]
    end
end

@testset "Checkpoint-shortened steps do not trigger adaptive fallback" begin
    config = het.Config(
        domain = (0.0, 0.08), discharge_voltage = 300.0,
        thruster = het.SPT_100, anode_mass_flow_rate = 5e-6,
    )
    sim = het.SimParams(
        grid = het.EvenGrid(20), duration = 1e-12, dt = 1e-9,
        num_save = 6, adaptive = true, verbose = false, max_small_steps = 2,
    )
    params = het.setup_simulation(config, sim)
    solver = het.SolverState(params; tspan = (0.0, sim.duration))
    for _ in 1:10
        params.cache.dt_limit[] = sim.min_dt / 2
        more_steps = het.step_solver!(params, config, solver)
        @test params.dt[] < sim.min_dt
        @test solver.small_steps == 0
        @test solver.fallback_steps_remaining == 0
        more_steps || break
    end
    @test solver.retcode == :success
    @test solver.times == solver.saveat
    @test length(solver.frames) == sim.num_save
end

@testset "Early termination returns only completed checkpoints" begin
    for retcode in (:failure, :error)
        @testset "$retcode" begin
            source = function (fluids, params)
                if params.iteration[] == 3
                    retcode == :error && error("Intentional source failure")
                    fluids.continuity[1].dens_ddt[3] = NaN
                end
                return nothing
            end
            config = het.Config(
                domain = (0.0, 0.08), discharge_voltage = 300.0,
                thruster = het.SPT_100, anode_mass_flow_rate = 5e-6,
                source_heavy_species = source,
            )
            sim = het.SimParams(
                grid = het.EvenGrid(20), duration = 1e-9, dt = 5e-10,
                num_save = 6, adaptive = false, verbose = false, print_errors = false,
            )
            sol = het.run_simulation(config, sim)
            expected_times = collect(range(0.0, sim.duration; length = sim.num_save))
            @test sol.retcode == retcode
            @test length(sol.t) == length(sol.frames) == 3
            @test sol.t == expected_times[1:3]
            @test length(het.time_average(sol).frames) == 1
        end
    end
end
