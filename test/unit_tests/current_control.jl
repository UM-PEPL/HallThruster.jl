using HallThruster: HallThruster as het
using Test

@testset "PID derivative with unequal timesteps" begin
    for error_function in (t -> 3t, t -> t^2)
        pid = het.PIDController(
            target_value = 0.0, proportional_constant = 0.0,
            derivative_constant = 2.0,
        )
        control = 0.0
        time = 0.0
        previous_error = 0.0
        # Establish the initial measurement before measuring a time derivative.
        control = het.apply_controller(pid, 0.0, control, 0.1)
        for dt in (0.1, 0.3, 0.005, 0.2)
            time += dt
            error = error_function(time)
            control = het.apply_controller(pid, -error, control, dt)
            # A derivative-only controller returns Kd times the latest slope,
            # independently of how the preceding interval was subdivided.
            @test control ≈ 2.0 * (error - previous_error) / dt rtol = 1.0e-11
            previous_error = error
        end
        encoded = het.serialize(pid)
        @test !haskey(encoded, "previous_dt")
        @test !haskey(encoded, "previous_error_slope")
        restored = het.deserialize(het.CurrentController, encoded)
        @test restored.previous_dt == 0.0
        @test restored.previous_error_slope == 0.0
        @test restored.errors == (0.0, 0.0, 0.0)
    end

    pid = het.PIDController(target_value = 4.0, derivative_constant = 1.0)
    control = het.apply_controller(pid, 0.0, 0.0, 100 * eps(Float64))
    @test control == 4.0
    @test het.apply_controller(pid, 0.0, control, 1.0e-9) == control
end

@testset "PID measurement smoothing" begin
    pid = het.PIDController(target_value = 0.0, smoothing_frequency = log(2.0))
    control = het.apply_controller(pid, 4.0, 0.0, 1.0)
    @test control ≈ -2.0
    @test het.apply_controller(pid, 4.0, control, 1.0) ≈ -3.0
end

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

function test_controller_serialization()
    return @testset "Serialization" begin
        none = het.NoController()
        test_roundtrip(het.CurrentController, none)

        pid = het.PIDController(target_value = 15.0)
        test_roundtrip(het.CurrentController, pid)

        pid_dict = het.Serialization.OrderedDict(
            "type" => "PIDController",
            "target_value" => 15.0,
            "integral_constant" => 1.0,
        )
        pid2 = het.deserialize(het.CurrentController, pid_dict)
        @test pid2 isa het.PIDController
        @test pid2.target_value == pid_dict["target_value"]
        @test pid2.integral_constant == pid_dict["integral_constant"]
        test_roundtrip(het.CurrentController, pid_dict)
    end
end

function test_pid_steady_state()
    return @testset "PID controller at steady state" begin
        # A PID controller that is already at the target value should effect no change
        # to the control variable
        target = 14.3
        controller = het.PIDController(
            target_value = target,
            proportional_constant = 1.5,
            integral_constant = 4.0,
            derivative_constant = 10.0,
        )

        starting_control = 1.0
        control = starting_control
        dt = 1.0
        val = target

        for _ in 1:10
            control = het.apply_controller(controller, val, control, dt)
        end

        @test control == starting_control
        @test controller.errors == (0.0, 0.0, 0.0)
    end
end

function update_springmass_system(state, controller, k, c, m, dt)
    # integrate spring-mass system
    # state: [x,v,f]
    # k: spring constant
    # c: damping constant
    # m: mass
    # f: external force
    # dt: timestep

    # differential equation:
    # m*x'' + c*x' + k*x = f(t)
    # => x'' = f/m + c/m*x' + k/m*x
    #
    # state space form:
    # x' = v
    # v' = (f - cv - kx)/m
    #
    # implicit discretization
    # x(t+dt) - x(t) = v(t) * dt
    # v(t+dt) - v(t) = (f - c*v(t+dt) - k*x(t+dt))/m * dt

    (x, v, f) = state

    # update position
    x = x + v * dt

    # apply force from PID controller
    f = het.apply_controller(controller, x, f, dt)

    # update velocity (implicitly)
    v = (v + (f - k * x) * (dt / m)) / (1 + dt * c / m)

    return (x, v, f)
end

function integrate_spring_mass_system(initial_condition, controller, k, c, m, t)
    N = length(t)
    x = zeros(N)
    v = zeros(N)
    (x[1], v[1], f) = initial_condition

    for i in 1:(N - 1)
        dt = t[i + 1] - t[i]
        x[i + 1], v[i + 1], f = update_springmass_system(
            (x[i], v[i], f), controller, k, c, m, dt,
        )
    end

    return x, v, f
end

function test_pid_springmass()
    # critical damping -> c = sqrt(4mk)
    k = 1
    m = 2
    c = sqrt(4 * m * k)
    dt = 0.2
    tmax = 240.0
    atol = 1.0e-6
    t = 0:dt:tmax

    target_value = 0.5

    # exact solutions (c.f https://content.oss.deltares.nl/delft3d/pid-controller.pdf)
    # for x(0) = 1, v(0) = 0, k = 1, m = 10, c = 2
    # K_p = 10.0, K_i = 0, K_d = 0 -> x = 0.45
    # K_p = 0, K_i = 0.05, K_d = 0 -> x = 0.5
    # K_p = 0, K_i = 0, K_d = 10.0 -> x = 0.0
    @testset "P only" begin
        # equilibrium: f = kx
        # Proportional only: f = K_p*(target-x)
        # K_p (target - x) = kx
        # => (K_p + k) x = target*K_p
        # => x =  0.5 K_p / (K_p + k)
        # Want x = 0.4 at equilibrium
        # => K_p = kx / (target - x) = x / 0.1 = 4

        equilibrium_x = 0.4
        K_p = k * equilibrium_x / (target_value - equilibrium_x)

        p_controller = het.PIDController(;
            target_value,
            proportional_constant = K_p,
            integral_constant = 0,
            derivative_constant = 0.0,
        )

        x, v, f = integrate_spring_mass_system((1.0, 0.0, 0.0), p_controller, k, c, m, t)

        @test isapprox(x[end], equilibrium_x; atol)
        @test isapprox(v[end], 0.0; atol)
        @test isapprox(f, k * x[end]; atol)
    end

    @testset "I only" begin
        # equilibrium: f = K_i ∫(target-x)dx = kx
        # => x = target

        equilibrium_x = target_value

        i_controller = het.PIDController(;
            target_value,
            proportional_constant = 0.0,
            integral_constant = 0.3,
            derivative_constant = 0.0,
        )

        x, v, f = integrate_spring_mass_system((1.0, 0.0, 0.0), i_controller, k, c, m, t)

        @test isapprox(x[end], equilibrium_x; atol)
        @test isapprox(v[end], 0.0; atol)
        @test isapprox(f, k * x[end]; atol)
    end

    @testset "D only" begin
        # equilibrium: f = K_p d(target-x)/dt = kx
        # => x = 0
        equilibrium_x = 0.0

        d_controller = het.PIDController(;
            target_value,
            proportional_constant = 0.0,
            integral_constant = 0.0,
            derivative_constant = 10.0,
        )

        x, v, f = integrate_spring_mass_system((1.0, 0.0, 0.0), d_controller, k, c, m, t)

        @test isapprox(x[end], equilibrium_x; atol)
        @test isapprox(v[end], 0.0; atol)
        @test isapprox(f, k * x[end]; atol)
    end

    return @testset "Tuned" begin
        equilibrium_x = target_value

        controller = het.PIDController(;
            target_value,
            proportional_constant = 4.0,
            integral_constant = 0.3,
            derivative_constant = 10.0,
        )

        x, v, f = integrate_spring_mass_system((1.0, 0.0, 0.0), controller, k, c, m, t)

        @test isapprox(x[end], equilibrium_x; atol)
        @test isapprox(v[end], 0.0; atol)
        @test isapprox(f, k * x[end]; atol)
    end
end

function test_current_control_bohmlike()
    return @testset "Bohm current-control multiplier does not compound" begin
        for multiplier in (0.5, 2.0)
            @testset "multiplier = $multiplier" begin
                config = het.Config(
                    domain = (0.0, 0.08),
                    discharge_voltage = 300.0,
                    thruster = het.SPT_100,
                    anom_model = het.TwoZoneBohm(1 / 160, 1 / 16),
                    propellants = [het.Propellant(het.Xenon, flow_rate_kg_s = 5.0e-6)],
                    anom_smoothing_iters = 0,
                )
                sim = het.SimParams(
                    duration = 1.0e-8,
                    dt = 1.0e-10,
                    adaptive = false,
                    grid = het.EvenGrid(20),
                    num_save = 2,
                    verbose = false,
                    current_control = het.NoController(),
                )
                params = het.setup_simulation(config, sim)
                # Evaluate the reference before the static model stops updating.
                base = config.anom_model(zeros(length(params.grid.cell_centers)), params)
                params.cache.anom_multiplier[] = multiplier
                solver = het.SolverState(params; tspan = (0.0, sim.duration))

                for _ in 1:50
                    @test het.step_solver!(params, config, solver)
                    @test params.cache.anom_multiplier[] ≈ multiplier
                    @test params.cache.νan_base ≈ base
                    @test params.cache.νan ≈ multiplier .* base
                end
                @test params.iteration[] == 51
                @test solver.retcode == :success
            end
        end
    end
end

function test_integral_current_control_bohmlike()
    return @testset "Integral current control updates the log multiplier" begin
        config = het.Config(
            domain = (0.0, 0.08),
            discharge_voltage = 300.0,
            thruster = het.SPT_100,
            propellants = [het.Propellant(het.Xenon, flow_rate_kg_s = 5.0e-6)],
        )
        controller = het.PIDController(
            target_value = 4.5,
            # Isolate the slow integral response from the proportional startup kick.
            proportional_constant = 0.0,
            integral_constant = 1.01 / 1.0e-4,
        )
        sim = het.SimParams(
            duration = 1.0e-8,
            dt = 1.0e-10,
            adaptive = false,
            grid = het.EvenGrid(20),
            num_save = 2,
            verbose = false,
            current_control = controller,
        )
        params = het.setup_simulation(config, sim)
        base = config.anom_model(zeros(length(params.grid.cell_centers)), params)
        # Setup invokes the controller once, starting with zero current and unit multiplier.
        setup_dt = min(het.initial_timestep(sim), 100 * eps(Float64))
        @test log(params.cache.anom_multiplier[]) ≈
            controller.integral_constant * controller.target_value * setup_dt atol = 1.0e-14
        solver = het.SolverState(params; tspan = (0.0, sim.duration))

        for _ in 1:50
            previous_log_multiplier = log(params.cache.anom_multiplier[])
            previous_current = params.cache.Id[]

            @test het.step_solver!(params, config, solver)
            expected_log_multiplier = previous_log_multiplier +
                controller.integral_constant *
                (controller.target_value - previous_current) * params.dt[]
            multiplier = params.cache.anom_multiplier[]
            @test isfinite(multiplier) && multiplier > 0
            @test log(multiplier) ≈ expected_log_multiplier atol = 1.0e-14
            @test params.cache.νan_base ≈ base
            @test params.cache.νan ≈ multiplier .* base
        end
        @test params.iteration[] == 51
        @test solver.retcode == :success
    end
end

test_current_control_bohmlike()
test_integral_current_control_bohmlike()
test_controller_serialization()
test_pid_steady_state()
test_pid_springmass()
