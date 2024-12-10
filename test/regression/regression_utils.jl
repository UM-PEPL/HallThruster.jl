using HallThruster: HallThruster as het
using Printf
using Test

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
end

function check_regression_case(case)
    (; file, thrust, current, ion_current) = case

    @testset "$(file)" begin
        println("======================================")
        println("        \"$(file)\"                   ")
        println("======================================")

        if haskey(case, :landmark_case)
            nsave = 1000
            sol_info = @timed run_landmark(
                1e-3; ncells = 150, nsave = nsave, case = case.landmark_case, CFL = case.CFL,
            )
        else
            file = "$(het.TEST_DIR)/regression/$(file)"
            sol_info = @timed het.run_simulation(file)
        end

        sol = sol_info.value

        nsave = length(sol.frames)
        avg_start = nsave ÷ 3
        n_avg = nsave - avg_start
        avg = het.time_average(sol, avg_start)

        T = het.thrust(sol) .* 1000
        T = [het.thrust(sol, i) for i in avg_start:nsave] .* 1000
        T_mean = het.mean(T)
        T_err = het.std(T) / sqrt(n_avg)
        Id = [het.discharge_current(sol, i) for i in avg_start:nsave]
        Id_mean = het.mean(Id)
        Id_err = het.std(Id) / sqrt(n_avg)
        ji = [het.ion_current(sol, i) for i in avg_start:nsave]
        ji_mean = het.mean(ji)
        ji_err = het.std(ji) / sqrt(n_avg)

		println("Performance:")
        @printf("Thrust: %.3f ± %.3f mN (expected %.3f mN)\n",
            T_mean, T_err, thrust)
        @printf("Discharge current: %.3f ± %.3f A (expected %.3f A)\n",
            Id_mean, Id_err, current)
        @printf("Ion current: %.3f ± %.3f A (expected %.3f A)\n",
            ji_mean, ji_err, ion_current)
        @test isapprox(thrust, T_mean, atol = T_err)
        @test isapprox(current, Id_mean, atol = Id_err)
        @test isapprox(ion_current, ji_mean, atol = ji_err)

		efficiencies = Dict(
			"Mass" => het.mass_eff,
			"Current" => het.current_eff,
			"Voltage" => het.voltage_eff,
			"Divergence" => het.divergence_eff,
			"Anode" => het.anode_eff,
		)

		println("\nEfficiencies:")

		for (eff_name, eff_func) in efficiencies
			eff = eff_func(sol)
			eff_mean = het.mean(eff)
			eff_err = het.std(eff) / sqrt(n_avg)
			eff_expected = case.efficiencies[eff_name]
			@printf("%s: %.1f ±  %.1f%% (expected %.1f%%)\n", 
		   		eff_name, eff_mean * 100, eff_err * 100, eff_expected * 100)
			@test isapprox(eff_mean, eff_expected, rtol = 1e-2)
		end

        max_Te = maximum(avg[:Tev][])
        max_E = maximum(avg[:E][])
        max_nn = maximum(avg[:nn][])
        max_ni = maximum(avg[:ni][])

		println("\nPlasma properties:")
        @printf("Peak electron temp: %.3f eV (expected %.3f eV)\n",
            max_Te, case.max_Te)
        @printf("Peak electric field: %.3e V/m (expected %.3e V/m)\n",
            max_E, case.max_E)
        @printf("Peak neutral density: %.3e m^-3 (expected %.3e m^-3)\n",
            max_nn, case.max_nn)
        @printf("Peak ion density: %.3e m^-3 (expected %.3e m^-3)\n",
            max_ni, case.max_ni)
        println()
        @test sol.retcode == :success

        @test isapprox(max_Te, case.max_Te, rtol = 1e-2)
        @test isapprox(max_E, case.max_E, rtol = 1e-2)
        @test isapprox(max_nn, case.max_nn, rtol = 1e-2)
        @test isapprox(max_ni, case.max_ni, rtol = 1e-2)
    end
end
