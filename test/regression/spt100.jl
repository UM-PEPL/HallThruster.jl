using HallThruster: HallThruster as het
using Test
using Printf

function test_spt100_case(case)
    (; file, thrust, current, ion_current) = case

    @testset "$(file)" begin
        file = "$(het.TEST_DIR)/regression/$(file)"

        sol_info = @timed het.run_simulation(file)
        sol = sol_info.value

        nsave = length(sol.savevals)
        avg_start = nsave ÷ 3
        n_avg = nsave - avg_start
        avg = het.time_average(sol, avg_start)

        T = HallThruster.thrust(sol) .* 1000
        T = [HallThruster.thrust(sol, i) for i in avg_start:nsave] .* 1000
        T_mean = het.mean(T)
        T_err = het.std(T) / sqrt(n_avg)
        Id = [HallThruster.discharge_current(sol, i) for i in avg_start:nsave]
        Id_mean = het.mean(Id)
        Id_err = het.std(Id) / sqrt(n_avg)
        ji = [HallThruster.ion_current(sol, i) for i in avg_start:nsave]
        ji_mean = het.mean(ji)
        ji_err = het.std(ji) / sqrt(n_avg)

        max_Te = maximum(avg[:Tev][])
        max_E = maximum(avg[:E][])
        max_nn = maximum(avg[:nn][])
        max_ni = maximum(avg[:ni][])

        println("======================================")
        println("           Case $(case.file)          ")
        println("======================================")
        @printf("Thrust: %.3f ± %.3f mN (expected %.3f mN)\n",
            T_mean, T_err, thrust)
        @printf("Discharge current: %.3f ± %.3f A (expected %.3f A)\n",
            Id_mean, Id_err, current)
        @printf("Ion current: %.3f ± %.3f A (expected %.3f A)\n",
            ji_mean, ji_err, ion_current)
        @printf("Peak electron temp: %.3f eV (expected %.3f eV)\n",
            max_Te, case.max_Te)
        @printf("Peak electric field: %.3g V/m (expected %.3g V/m)\n",
            max_E, case.max_E)
        @printf("Peak neutral density: %.3e m^-3 (expected %.3e m^-3)\n",
            max_nn, case.max_nn)
        @printf("Peak ion density: %.3e m^-3 (expected %.3e m^-3)\n",
            max_ni, case.max_ni)
        println()
        @test sol.retcode == :success

        @test isapprox(thrust, T_mean, atol = T_err)
        @test isapprox(current, Id_mean, atol = Id_err)
        @test isapprox(ion_current, ji_mean, atol = ji_err)
        @test isapprox(max_Te, case.max_Te, rtol = 1e-2)
        @test isapprox(max_E, case.max_E, rtol = 1e-2)
        @test isapprox(max_nn, case.max_nn, rtol = 1e-2)
        @test isapprox(max_ni, case.max_ni, rtol = 1e-2)
    end
end

function test_spt100_regression()
    @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 85.412,
            current = 4.612,
            ion_current = 3.925,
            max_Te = 24.078,
            max_E = 63010.527,
            max_nn = 2.096424e19,
            max_ni = 8.74894e17,
        )
        test_spt100_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 102.769,
            current = 4.989,
            ion_current = 4.045,
            max_Te = 27.789,
            max_E = 96495.372,
            max_nn = 2.1017e19,
            max_ni = 1.013126e18,
        )
        test_spt100_case(with_plume)
    end
end
