using HallThruster: HallThruster as het
using Printf
using Test
using CairoMakie: Makie as mk
using DelimitedFiles

function load_landmark_data(case)
    codes = ["fluid_1", "fluid_2", "hybrid"]
    fields = ["electric_field", "energy", "ionization", "neutral_density", "plasma_density", "potential"]

    field_type = Dict{String, @NamedTuple{x::Vector{Float64}, y::Vector{Float64}}}
    output = Dict{String, field_type}()

    for field in fields
        f = field_type()
        for code in codes
            x = readdlm(joinpath(het.LANDMARK_FOLDER, "case_$(case)", "$(field)_$(code).csv"), ',')
            f[code] = (; x = x[:, 1], y = x[:, 2])
        end
        output[field] = f
    end
    return output
end

function plot_landmark(sol, case)

    avg = het.time_average(sol, 5.0e-4)
    data = load_landmark_data(case)

    fig = mk.Figure()
    xlabel = "z [cm]"
    ax_nn = mk.Axis(fig[1, 1]; xlabel, ylabel = "Density [10¹⁹ m³/s]")
    ax_ne = mk.Axis(fig[1, 2]; xlabel, ylabel = "Density [10¹⁸ m³/s]")
    ax_E = mk.Axis(fig[2, 1]; xlabel, ylabel = "Electric field [kV/m]")
    ax_Te = mk.Axis(fig[2, 2]; xlabel, ylabel = "Electron temperature [eV]")

    function to_title(label)
        return titlecase(replace(label, "_" => " "))
    end

    lines = []
    labels = String[]

    colors = [:red, :green, :blue]
    style = (; linewidth = 2, color = :black)

    for (i, (k, v)) in enumerate(pairs(data["neutral_density"]))
        l = mk.lines!(ax_nn, v.x, v.y ./ 1.0e19; color = colors[i], linestyle = :dash)
        push!(lines, l)
        push!(labels, to_title(k))
    end
    l = mk.lines!(ax_nn, avg[:z], avg[:nn][] ./ 1.0e19; style...)
    push!(lines, l)
    push!(labels, "HallThruster.jl")

    for (i, (k, v)) in enumerate(pairs(data["plasma_density"]))
        mk.lines!(ax_ne, v.x, v.y ./ 1.0e18; color = colors[i], linestyle = :dash)
    end
    mk.lines!(ax_ne, avg[:z], avg[:ne][] ./ 1.0e18; style...)

    for (i, (k, v)) in enumerate(pairs(data["electric_field"]))
        mk.lines!(ax_E, v.x, v.y ./ 1000; color = colors[i], linestyle = :dash)
    end
    mk.lines!(ax_E, avg[:z], avg[:E][] ./ 1000; style...)

    for (i, (k, v)) in enumerate(pairs(data["energy"]))
        mk.lines!(ax_Te, v.x, v.y / 1.5; color = colors[i], linestyle = :dash)
    end
    mk.lines!(ax_Te, avg[:z], avg[:Tev][]; style...)

    mk.Legend(fig[0, :], lines, labels, orientation = :horizontal)

    return mk.save("landmark_$(case).png", fig)
end

function run_landmark(
        duration = 1.0e-3; ncells = 200, nsave = 2, dt = 0.7e-8, CFL = 0.799, case = 1,
    )
    domain = (0.0, 0.05)

    #Landmark cases loss frequencies
    αϵ_in, αϵ_out = if case == 1
        (1.0, 1.0)
    elseif case == 2
        (0.5, 1.0)
    elseif case == 3
        (0.4, 1.0)
    end

    ϵ_anode = 3.0
    ϵ_cathode = 3.0

    config = het.Config(;
        ncharge = 1,
        domain,
        reconstruct = case > 1,
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
        anode_mass_flow_rate = 5.0e-6,
        transition_length = 1.0e-3,
        ion_wall_losses = false,
        anom_model = het.TwoZoneBohm(1 / 160, 1 / 16),
        anode_boundary_condition = :dirichlet,
        conductivity_model = het.LANDMARK_conductivity(),
    )

    @time sol = het.run_simulation(
        config; duration, grid = het.UnevenGrid(ncells), nsave,
        dt, dtmin = dt / 100, dtmax = dt * 10, adaptive = true, CFL, verbose = false,
    )

    plot_landmark(sol, case)

    return sol
end

function check_regression_case(case)
    (; file, thrust, current, ion_current) = case

    return @testset "$(file)" begin
        println("======================================")
        println("        \"$(file)\"                   ")
        println("======================================")

        if haskey(case, :landmark_case)
            nsave = 1000
            ncells = case.landmark_case == 1 ? 250 : 150
            sol_info = @timed run_landmark(
                1.0e-3; ncells, nsave = nsave, case = case.landmark_case, CFL = case.CFL,
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
        @printf(
            "Thrust: %.3f ± %.3f mN (expected %.3f mN)\n",
            T_mean, T_err, thrust
        )
        @printf(
            "Discharge current: %.3f ± %.3f A (expected %.3f A)\n",
            Id_mean, Id_err, current
        )
        @printf(
            "Ion current: %.3f ± %.3f A (expected %.3f A)\n",
            ji_mean, ji_err, ion_current
        )
        @test isapprox(thrust, T_mean, atol = T_err)
        @test isapprox(current, Id_mean, atol = Id_err)
        @test isapprox(ion_current, ji_mean, atol = ji_err)

        efficiency_funcs = Dict(
            "Mass" => het.mass_eff,
            "Current" => het.current_eff,
            "Voltage" => het.voltage_eff,
            "Divergence" => het.divergence_eff,
            "Anode" => het.anode_eff,
        )

        println("\nEfficiencies:")

        efficiencies = Dict{String, Float64}()

        for (eff_name, eff_func) in efficiency_funcs
            eff = eff_func(sol)
            eff_mean = het.mean(eff)
            eff_err = het.std(eff) / sqrt(n_avg)
            eff_expected = case.efficiencies[eff_name]
            @printf(
                "%s: %.1f ±  %.1f%% (expected %.1f%%)\n",
                eff_name, eff_mean * 100, eff_err * 100, eff_expected * 100
            )
            @test isapprox(eff_mean, eff_expected, rtol = 1.0e-2)
            efficiencies[eff_name] = eff_mean
        end

        max_Te = maximum(avg[:Tev][])
        max_E = maximum(avg[:E][])
        max_nn = maximum(avg[:nn][])
        max_ni = maximum(avg[:ni][])

        println("\nPlasma properties:")
        @printf(
            "Peak electron temp: %.3f eV (expected %.3f eV)\n",
            max_Te, case.max_Te
        )
        @printf(
            "Peak electric field: %.3e V/m (expected %.3e V/m)\n",
            max_E, case.max_E
        )
        @printf(
            "Peak neutral density: %.3e m^-3 (expected %.3e m^-3)\n",
            max_nn, case.max_nn
        )
        @printf(
            "Peak ion density: %.3e m^-3 (expected %.3e m^-3)\n",
            max_ni, case.max_ni
        )
        println()
        @test sol.retcode == :success

        @test isapprox(max_Te, case.max_Te, rtol = 1.0e-2)
        @test isapprox(max_E, case.max_E, rtol = 1.0e-2)
        @test isapprox(max_nn, case.max_nn, rtol = 1.0e-2)
        @test isapprox(max_ni, case.max_ni, rtol = 1.0e-2)

        # print output to replace what we had
        println("-----")
        println("Replace contents of benchmark with the following if necessary")
        println("-----")
        @printf("thrust = %.3f,\n", T_mean)
        @printf("current = %.3f,\n", Id_mean)
        @printf("ion_current = %.3f,\n", ji_mean)
        @printf("max_Te = %.3f,\n", max_Te)
        @printf("max_E = %.3e,\n", max_E)
        @printf("max_nn = %.3e,\n", max_nn)
        @printf("max_ni = %.3e,\n", max_ni)
        @printf("efficiencies = Dict(\n")
        @printf("   \"Mass\" => %.3f,\n", efficiencies["Mass"])
        @printf("   \"Current\" => %.3f,\n", efficiencies["Current"])
        @printf("   \"Divergence\" => %.3f,\n", efficiencies["Divergence"])
        @printf("   \"Voltage\" => %.3f,\n", efficiencies["Voltage"])
        @printf("   \"Anode\" => %.3f,\n", efficiencies["Anode"])
        @printf("),\n")
        println("-----")
    end
end
