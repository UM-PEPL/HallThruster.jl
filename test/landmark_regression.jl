using Printf, Test
using HallThruster: HallThruster as het

function run_landmark(
        duration = 1e-3; ncells = 200, nsave = 2, dt = 0.7e-8, CFL = 0.799, case = 1,)
    domain = (0.0, 0.05)

    #Landmark cases loss frequencies
    αϵ_in, αϵ_out = if case == 1
        (1.0, 1.0)
    elseif case == 2
        (0.5, 1.0)
    elseif case == 3
        (0.4, 1.0)
    end

    scheme = het.HyperbolicScheme(
        # We use global_lax_friedrichs here to better handle case 1, as it is very oscillatory and this
        # scheme is the most diffusive
        # In general, prefer rusanov or HLLE
        flux_function = het.global_lax_friedrichs,
        limiter = het.minmod,
        reconstruct = true,
    )

    ϵ_anode = 3.0
    ϵ_cathode = 3.0

    config = het.Config(;
        ncharge = 1,
        scheme,
        domain,
        anode_Tev = 2 / 3 * ϵ_anode,
        cathode_Tev = 2 / 3 * ϵ_cathode,
        discharge_voltage = 300.0,
        ionization_model = :Landmark,
        excitation_model = :Landmark,
        electron_neutral_model = :Landmark,
        electron_ion_collisions = false,
        wall_loss_model = het.ConstantSheathPotential(20, αϵ_in, αϵ_out),
        LANDMARK = true,
        neutral_velocity = 150.0,
        ion_temperature_K = 0.0,
        thruster = het.SPT_100,
        anode_mass_flow_rate = 5e-6,
        transition_length = 1e-3,
        ion_wall_losses = false,
        anom_model = het.TwoZoneBohm(1 / 160, 1 / 16),
        anode_boundary_condition = :dirichlet,
        conductivity_model = het.LANDMARK_conductivity(),
    )

    @time sol = het.run_simulation(
        config; duration, grid = het.EvenGrid(ncells), nsave,
        dt, dtmin = dt / 100, dtmax = dt * 10, adaptive = true, CFL, verbose = false,
    )
    return sol
end;

mean(x) = sum(x) / length(x)
function var(x)
    μ = mean(x)
    return mean((_x - μ)^2 for _x in x)
end
std(x) = sqrt(var(x))

@testset "LANDMARK regression tests" begin
    CFLs = [0.25, 0.799, 0.799]
    expected_thrusts = [90.135, 100.540, 98.406]
    expected_currents = [7.636, 8.027, 7.767]
    expected_ion_currents = [3.719, 3.642, 3.612]

    for (i, (CFL, thrust, current, ion_current)) in enumerate(zip(
        CFLs, expected_thrusts, expected_currents, expected_ion_currents,))
        nsave = 1000
        avg_start = 250
        n_avg = nsave - avg_start
        println("======================================")
        println("               Case $i                ")
        println("======================================")
        sol_info = @timed run_landmark(
            1e-3; ncells = 150, nsave = nsave, case = i, CFL = CFL,)
        sol = sol_info.value
        time = sol_info.time
        T = [het.thrust(sol, i) for i in avg_start:nsave] .* 1000
        T_mean = mean(T)
        T_err = std(T) / sqrt(n_avg)
        Id = [het.discharge_current(sol, i) for i in avg_start:nsave]
        Id_mean = mean(Id)
        Id_err = std(Id) / sqrt(n_avg)
        ji = [het.ion_current(sol, i) for i in avg_start:nsave]
        ji_mean = mean(ji)
        ji_err = std(ji) / sqrt(n_avg)
        @printf("Thrust: %.3f ± %.3f mN (expected %.3f mN)\n", T_mean, T_err, thrust)
        @printf("Discharge current: %.3f ± %.3f A (expected %.3f A)\n", Id_mean, Id_err,
            current)
        @printf("Ion current: %.3f ± %.3f A (expected %.3f A)\n", ji_mean, ji_err,
            ion_current)
        println()
        @test sol.retcode == :success
        @test isapprox(thrust, T_mean, atol = T_err)
        @test isapprox(current, Id_mean, atol = Id_err)
        @test isapprox(ion_current, ji_mean, atol = ji_err)
    end
end
