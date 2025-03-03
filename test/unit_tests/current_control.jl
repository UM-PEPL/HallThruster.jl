using HallThruster: HallThruster as ht
using Test

function test_controller_serialization()
    @testset "Serialization" begin
        none = ht.NoController()
        test_roundtrip(het.CurrentController, none)

        pid = ht.PIDController(target_value = 15.0)
        test_roundtrip(het.CurrentController, pid)

        pid_dict = ht.Serialization.OrderedDict(
            :type => "PIDController",
            :target_value => 15.0,
            :integral_constant => 1.0,
        )
        pid2 = ht.deserialize(ht.CurrentController, pid_dict)
        @test pid2 isa ht.PIDController
        @test pid2.target_value == pid_dict[:target_value]
        @test pid2.integral_constant == pid_dict[:integral_constant]
        test_roundtrip(ht.CurrentController, pid_dict)
    end
end

function test_pid_steady_state()
    @testset "PID controller at steady state" begin
        # A PID controller that is already at the target value should effect no change
        # to the control variable
        target = 14.3
        controller = ht.PIDController(
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
            control = ht.apply_controller(controller, val, control, dt)
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
    f = ht.apply_controller(controller, x, f, dt)

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
    atol = 1e-6
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

        p_controller = ht.PIDController(;
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

        i_controller = ht.PIDController(;
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

        d_controller = ht.PIDController(;
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

    @testset "Tuned" begin
        equilibrium_x = target_value

        controller = ht.PIDController(;
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

function test_current_control()
    @testset "Current control" begin
        test_controller_serialization()
        test_pid_steady_state()
        @testset "Spring-mass system" begin
            test_pid_springmass()
        end
    end
end
