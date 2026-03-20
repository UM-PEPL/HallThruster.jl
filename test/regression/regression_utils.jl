using HallThruster: HallThruster as het, JSON
using Printf
using Test
using CairoMakie: Makie as mk
using DelimitedFiles

struct Oscillations
    time::Vector{Float64}
    thrust::Vector{Float64}
    discharge_current::Vector{Float64}
    ion_current::Vector{Float64}
end

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

to_title(label) = uppercasefirst(replace(label, "_" => " "))

function plot_sim(avg, ref, filename, landmark_case = nothing)
    fig = mk.Figure(size = (1200, 600))
    xlabel = "z [cm]"
    ax_nn = mk.Axis(fig[1, 1]; xlabel, ylabel = "Neutral density (10¹⁹ m³/s)")
    ax_ne = mk.Axis(fig[2, 1]; xlabel, ylabel = "Plasma density [10¹⁸ m³/s)")
    ax_E = mk.Axis(fig[1, 2]; xlabel, ylabel = "Electric field (kV/m)")
    ax_Te = mk.Axis(fig[2, 2]; xlabel, ylabel = "Electron temperature (eV)")
    ax_ui = mk.Axis(fig[1, 3]; xlabel, ylabel = "Ion velocity (km/s)")
    ax_ue = mk.Axis(fig[2, 3]; xlabel, ylabel = "Electron velocity (km/s)")
    ax_phi = mk.Axis(fig[1, 4]; xlabel, ylabel = "Electric potential (V)")

    function plot_sim_fields!(lines, labels, sim; style, label = "")
        frame = sim.frames[]
        z = sim.grid

        l = mk.lines!(ax_nn, z, frame.neutrals[:Xe].n ./ 1.0e19; style...)
        push!(lines, l)
        push!(labels, label)

        mk.lines!(ax_ne, z, frame.ne ./ 1.0e18; style...)
        mk.lines!(ax_E, z, frame.E ./ 1000; style...)
        mk.lines!(ax_Te, z, frame.Tev; style...)
        mk.lines!(ax_ui, z, frame.ions[:Xe][1].u ./ 1000; style...)
        mk.lines!(ax_ue, z, frame.ue ./ 1000; style...)
        mk.lines!(ax_phi, z, frame.potential; style...)
        return lines, labels
    end

    lines = []
    labels = String[]

    if ref !== nothing
        plot_sim_fields!(lines, labels, ref; style = (; linewidth = 2, color = :red), label = "Reference")
    end

    plot_sim_fields!(lines, labels, avg; style = (; linewidth = 2, color = :black), label = "HallThruster.jl")

    if landmark_case !== nothing
        landmark_colors = [:red, :green, :blue]
        data = load_landmark_data(landmark_case)
        for (i, (k, v)) in enumerate(pairs(data["neutral_density"]))
            l = mk.lines!(ax_nn, v.x, v.y ./ 1.0e19; color = landmark_colors[i], linestyle = :dash)
            push!(lines, l)
            push!(labels, to_title(k))
        end
        for (i, (_, v)) in enumerate(pairs(data["plasma_density"]))
            mk.lines!(ax_ne, v.x, v.y ./ 1.0e18; color = landmark_colors[i], linestyle = :dash)
        end
        for (i, (_, v)) in enumerate(pairs(data["electric_field"]))
            mk.lines!(ax_E, v.x, v.y ./ 1000; color = landmark_colors[i], linestyle = :dash)
        end
        for (i, (_, v)) in enumerate(pairs(data["energy"]))
            mk.lines!(ax_Te, v.x, v.y / 1.5; color = landmark_colors[i], linestyle = :dash)
        end
        for (i, (_, v)) in enumerate(pairs(data["potential"]))
            mk.lines!(ax_phi, v.x, v.y; color = landmark_colors[i], linestyle = :dash)
        end
    end

    mk.Legend(fig[2, 4], lines, labels, tellwidth = false)
    mk.save(filename, fig)
    return
end

function plot_oscillations(osc::Oscillations, ref::Oscillations, filename::String)
    fig = mk.Figure(size = (1000, 800))
    ax_T = mk.Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Thrust (N)")
    ax_T_zoom = mk.Axis(fig[1, 2], xlabel = "Time (s)", ylabel = "Thrust (N)")
    ax_I = mk.Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Discharge current (A)")
    ax_I_zoom = mk.Axis(fig[2, 2], xlabel = "Time (s)", ylabel = "Discharge current (A)")
    ax_J = mk.Axis(fig[3, 1], xlabel = "Time (s)", ylabel = "Ion current (A)")
    ax_J_zoom = mk.Axis(fig[3, 2], xlabel = "Time (s)", ylabel = "Ion current (A)")

    lines = []
    labels = []

    start_ind = length(osc.time) * 3 ÷ 4

    for (sim, name, linestyle, color) in zip([osc, ref], ["Current", "Reference"], [:solid, :dash], [:black, :red])
        kwargs = (; linestyle, color)
        t = sim.time
        t_zoom = t[start_ind:end]

        l = mk.lines!(ax_T, t, sim.thrust; kwargs...)
        push!(lines, l)
        push!(labels, name)
        mk.lines!(ax_T_zoom, t_zoom, sim.thrust[start_ind:end]; kwargs...)

        mk.lines!(ax_I, t, sim.discharge_current; kwargs...)
        mk.lines!(ax_I_zoom, t_zoom, sim.discharge_current[start_ind:end]; kwargs...)

        mk.lines!(ax_J, t, sim.ion_current; kwargs...)
        mk.lines!(ax_J_zoom, t_zoom, sim.ion_current[start_ind:end]; kwargs...)
    end

    for (full, zoom) in zip([ax_T, ax_I, ax_J], [ax_T_zoom, ax_I_zoom, ax_J_zoom])
        mk.reset_limits!(zoom)
        zoom_lims = zoom.finallimits[]
        mk.poly!(full, zoom_lims, color = (:white, 0), strokecolor = :orange, strokewidth = 2)
    end

    mk.Legend(fig[0, :], lines, labels, orientation = :horizontal)
    mk.save(filename, fig)
    return
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
        config; duration, grid = het.EvenGrid(ncells), nsave,
        dt, dtmin = dt / 100, dtmax = dt * 10, adaptive = true, CFL, verbose = false,
    )

    return sol
end

function p2p(x)
    _min, _max = extrema(x)
    return _max - _min
end

function check_and_print(description, reduction, x::Vector{T}, y::Vector{T}; atol = zero(T), rtol = (atol > 0 ? zero(T) : sqrt(eps(T))), exponential = false) where {T <: Real}
    rx, ry = reduction(x), reduction(y)
    @test isapprox(rx, ry; rtol, atol)
    return if exponential
        @printf("%s: %.4g (expected %.4g)\n", description, rx, ry)
    else
        @printf("%s: %.3f (expected %.3f)\n", description, rx, ry)
    end
end

HEADER_WIDTH = 40

print_divider(c = '=', w = HEADER_WIDTH) = println(c^w)

function print_centered(s, w = HEADER_WIDTH)
    num_left = max(0, floor(Int, (w - length(s)) / 2))
    num_right = max(0, w - length(s) - num_left)
    return println(" "^num_left * s * " "^num_right)
end

function print_header(s, c = '='; w = HEADER_WIDTH, top = true, bottom = true)
    top && print_divider(c, w)
    print_centered(s, w)
    return bottom && print_divider(c, w)
end

function check_regression_case(case; fix = false)
    (; file) = case

    @testset "$(file)" begin
        REGRESSION_DIR = joinpath(het.TEST_DIR, "regression")
        OUTPUT_DIR = joinpath(REGRESSION_DIR, "output")
        mkpath(OUTPUT_DIR)

        print_header(file; w = 60)

        casename, landmark_case, sol_info = if haskey(case, :landmark_case)
            nsave = 1000
            ncells = 200
            sol_info = @timed run_landmark(2.0e-3; ncells, nsave = nsave, case = case.landmark_case, CFL = 0.5)
            "landmark_$(case.landmark_case)", case.landmark_case, sol_info
        else
            casename = splitext(basename(file))[1]
            sol_info = @timed het.run_simulation(joinpath(REGRESSION_DIR, file))
            splitext(basename(file))[1], nothing, sol_info
        end

        sol = sol_info.value
        avg_start_time = 1.5e-3
        nsave = sol.simulation.num_save
        avg_start_ind = floor(Int, nsave * avg_start_time / sol.t[end])
        n_avg = nsave - avg_start_ind
        avg = het.time_average(sol, avg_start_ind)

        ϕ = avg[:ϕ][]
        Tev = avg[:Tev][]

        # Check for sim success
        @test sol.retcode == :success

        # Check potential boundary condition
        @test ϕ[end] ≈ avg.config.cathode_coupling_voltage atol = 0.01

        # Check temperature boundary condition
        diff = abs(Tev[end] - Tev[end - 1])
        @test Tev[end] ≈ avg.config.cathode_Tev atol = 0.2 * diff

        thrust = het.thrust(sol)
        discharge_current = het.discharge_current(sol)
        ion_current = het.ion_current(sol)

        oscillations = Oscillations(sol.t, thrust, discharge_current, ion_current)

        # Overwrite existing benchmarks if we so choose
        COMPARISON_FILE = joinpath(OUTPUT_DIR, "ref_$(casename).json")
        OSCILLATIONS_FILE = joinpath(OUTPUT_DIR, "oscillations_$(casename).json")
        if fix
            open(joinpath(OUTPUT_DIR, COMPARISON_FILE), "w") do f
                JSON.write_json(f, het.serialize(avg))
            end
            open(OSCILLATIONS_FILE, "w") do f
                JSON.write_json(f, oscillations)
            end
        end

        # Load comparison files
        ref_sim_dict = JSON.parsefile(COMPARISON_FILE)
        ref_sim = het.deserialize(het.Solution, ref_sim_dict)
        ref_oscillations_dict = JSON.parsefile(OSCILLATIONS_FILE)

        ref_oscillations = Oscillations(
            ref_oscillations_dict["time"],
            ref_oscillations_dict["thrust"],
            ref_oscillations_dict["discharge_current"],
            ref_oscillations_dict["ion_current"]
        )

        # Plot comparison of time-averaged properties
        plot_sim(avg, ref_sim, joinpath(OUTPUT_DIR, "ref_$(casename).png"), landmark_case)

        # Replace landmark case plots if requested
        if fix && landmark_case !== nothing
            plot_sim(avg, nothing, joinpath(het.PACKAGE_ROOT, "docs/src/assets/$(casename).svg"), landmark_case)
        end

        # Plot comparison of oscillatory properties
        plot_oscillations(oscillations, ref_oscillations, joinpath(OUTPUT_DIR, "oscillations_$(casename).png"))

        @show avg_start_ind

        print_header("Oscillations", '-')
        for key in [:thrust, :discharge_current, :ion_current]
            sim_osc = getfield(oscillations, key)
            ref_osc = getfield(ref_oscillations, key)

            sim_osc_steady = sim_osc[avg_start_ind:end]
            ref_osc_steady = ref_osc[avg_start_ind:end]
            osc_err = het.std(ref_osc_steady) / sqrt(n_avg)

            # Check for 3 things - height of transient, p2p amplitude, and mean value
            name = to_title(string(key))
            check_and_print(name * " max", maximum, sim_osc, ref_osc; rtol = 0.05)
            check_and_print(name * " steady p2p", p2p, sim_osc_steady, ref_osc_steady; rtol = 0.05)
            check_and_print(name * " steady mean", het.mean, sim_osc_steady, ref_osc_steady; atol = osc_err)
        end
        println()

        print_header("Efficiencies", '-')
        efficiency_funcs = Dict(
            "Mass" => het.mass_eff,
            "Current" => het.current_eff,
            "Voltage" => het.voltage_eff,
            "Divergence" => het.divergence_eff,
            "Anode" => het.anode_eff,
        )
        for (eff_name, eff_func) in efficiency_funcs
            sim_eff = eff_func(avg)
            ref_eff = eff_func(ref_sim)
            check_and_print(eff_name, x -> x[], sim_eff, ref_eff, rtol = 0.01)
        end
        println()

        print_header("Plasma properties", '-')
        sim_frame = avg.frames[]
        ref_frame = ref_sim.frames[]
        field_names = ["Tev", "E", "ne", "nn"]
        field_fns = [:Tev, :E, :ne, x -> x.neutrals[:Xe].n]
        for Z in eachindex(ref_frame.ions[:Xe])
            push!(field_names, "ni_$(Z)")
            push!(field_fns, x -> x.ions[:Xe][Z].n)
        end
        reductions = [maximum, minimum, het.mean, only]
        reduction_names = ["max", "min", "mean", "L2"]
        for (field_name, field_fn) in zip(field_names, field_fns)
            if field_fn isa Symbol
                sim_qty = getfield(sim_frame, field_fn)
                ref_qty = getfield(ref_frame, field_fn)
            else
                sim_qty = field_fn(sim_frame)
                ref_qty = field_fn(ref_frame)
            end
            for (reduction_name, reduction) in zip(reduction_names, reductions)
                exponential = field_name != "Tev" && field_name != "E" && reduction_name != "L2"
                name = "$(field_name) $(reduction_name)"
                if reduction_name == "L2"
                    L2 = sqrt(sum((sim_qty .- ref_qty) .^ 2) / sum(ref_qty .^ 2))
                    check_and_print(name, reduction, [L2], [0.0]; atol = 0.01, exponential)
                else
                    check_and_print(name, reduction, sim_qty, ref_qty; rtol = 1.0e-2, exponential)
                end
            end
            println()
        end
    end

    return
end
""
