using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase, DataFrames, CSV, JLD2

using DelimitedFiles

function landmark_references(case, variable)
    if case > 3
        return String[]
    else
        suffixes = "fluid_1", "fluid_2", "hybrid"
        return [
            joinpath("landmark", "case_$case", "$(variable)_$(suffix).csv")
            for suffix in suffixes
        ]
    end
end

function landmark_styles()
    common_options = (
        linewidth = 1.5,
        linestyle = :dash,
    )
    return [
        (color = :red, label = "Fluid (δ = 1 mm)", common_options...),
        (color = :green, label = "Fluid (δ = 0.5 mm)", common_options...),
        (color = :blue, label = "Hybrid", common_options...)
    ]
end


function plot_quantity(u, z = nothing, zmin = 0.0, zmax = 0.05; ref_paths = String[], ref_styles = nothing, hallis = nothing, hallisvar = nothing, kwargs...)
    if isnothing(z)
        z = LinRange(zmin, zmax, length(u))
    end
    p = plot()

    for (i, ref_path) in enumerate(ref_paths)
        if ispath(ref_path)
            ref_data = readdlm(ref_path, ',')
            z_ref, q_ref = ref_data[:, 1], ref_data[:, 2]
            if isnothing(ref_styles)
                plot!(p, z_ref, q_ref, label = "Reference $i", lw = 1, lc = :red, ls = :dash)
            else
                plot!(p, z_ref, q_ref; ref_styles[i]...)
            end
        else
            @warn("File $ref_path not found, skipping plot.")
        end
    end
    if !isnothing(hallis)
        z_ref, q_ref = hallis.z, hallisvar
        plot!(p, z_ref, q_ref, label = "Reference", lw = 2, lc = :red, ls = :dash)
    end

    plot!(
        p, z, u; label = "HallThruster.jl", xlabel = "z (m)", legend = :outertop, margin = 8Plots.mm, lw = 2,
        color = :black, linestyle = :solid,
        kwargs...
    )
    return p
end

function load_hallis_output(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei,
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    return output[1:end-1, :]
end

function plot_solution(u, saved_values, z, case = 1)
    mi = HallThruster.Xenon.m
    coeff = HallThruster.load_landmark()
    (;Tev, ue, ϕ_cell, ∇ϕ, ne, pe, ∇pe) = saved_values
    #z_edge = [z[1]; [0.5 * (z[i] + z[i+1]) for i in 2:length(z)-2]; z[end]]
    ionization_rate = [coeff.rate_coeff(3/2 * Tev[i])*u[1, i]*ne[i]/mi for i in 1:size(u, 2)]

    ref_styles = landmark_styles()

    p_nn = plot_quantity(
        u[1, :] / mi, z; title = "Neutral density", ylabel = "nn (m⁻³)",
        ref_paths = landmark_references(case, "neutral_density"), ref_styles
    )

    p_ne = plot_quantity(
        ne, z; title = "Plasma density", ylabel = "ne (m⁻³)",
        ref_paths = landmark_references(case, "plasma_density"), ref_styles
    )

    p_ui = plot_quantity(u[3, :] ./ u[2, :] ./ 1000, z; title = "Ion velocity", ylabel = "ui (km/s)")

    p_iz = plot_quantity(
        ionization_rate, z; title = "Ionization rate", ylabel = "nϵ (eV m⁻³)",
        ref_paths = landmark_references(case, "ionization"), ref_styles
    )

    p_ϵ  = plot_quantity(
        u[4, :] ./ ne, z; title = "Electron energy (3/2 Te) (eV)", ylabel = "ϵ (eV)",
        ref_paths = landmark_references(case, "energy"), ref_styles
    )

    p_ue = plot_quantity(ue ./ 1000, z; title = "Electron velocity", ylabel = "ue (km/s)")
    p_ϕ  = plot_quantity(
        ϕ_cell, z; title = "Potential", ylabel = "ϕ (V)",
        ref_paths = landmark_references(case, "potential"), ref_styles
    )

    p_E  = plot_quantity(
        -∇ϕ, z; title = "Electric field", ylabel = "E (V/m)",
        ref_paths = landmark_references(case, "electric_field"), ref_styles
    )

    #p_pe  = plot_quantity(HallThruster.e * pe, z; title = "Electron pressure", ylabel = "∇pe (Pa)")
    #p_∇pe  = plot_quantity(HallThruster.e * ∇pe, z; title = "Pressure gradient", ylabel = "∇pe (Pa/m)")
    plot(p_nn, p_ne, p_ui, p_ϕ, #=p_pe,=# p_iz, p_ϵ, p_ue, p_E, #=p_∇pe,=# layout = (2, 4), size = (2500, 1000))
end

function plot_solution_OVS(u, z = nothing, case = 1)
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_0$(case).out")
    coeff = HallThruster.load_landmark()
    mi = HallThruster.Xenon.m
    ionization_rate = [coeff.rate_coeff(u[5, i])*u[1, i]*u[2, i]/mi/mi for i in 1:size(u, 2)]
    p_nn = plot_quantity(u[1, :] / mi, z; title = "Neutral density", ylabel = "nn (m⁻³)", hallis = hallis, hallisvar = hallis.nn)
    p_ne = plot_quantity(u[2, :] / mi, z; title = "Plasma density", ylabel = "ne (m⁻³)", hallis = hallis, hallisvar = hallis.ne)
    p_ui = plot_quantity(u[3, :] ./ u[2, :] ./ 1000, z; title = "Ion velocity", ylabel = "ui (km/s)")
    p_nϵ = plot_quantity(u[4, :] , z; title = "nepsilon", ylabel = "nepsilon", hallis = hallis, hallisvar = hallis.ne)
    p_ϵ  = plot_quantity(u[5, :], z; title = "Electron temperature (eV)", ylabel = "ϵ (eV)", hallis = hallis, hallisvar = hallis.Te)
    p_ue = plot_quantity(u[10, :] ./ 1000, z; title = "Electron velocity", ylabel = "ue (km/s)")
    p_ϕ  = plot_quantity(u[8, :], z; title = "Potential", ylabel = "ϕ (V)", hallis = hallis, hallisvar = hallis.ϕ)
    p_E  = plot_quantity(-u[9, :], z; title = "Electric field", ylabel = "E (V/m)", hallis = hallis, hallisvar = hallis.Ez)
    p = plot(p_nn, p_ne, p_ui, p_nϵ, p_ϵ, p_ue, p_ϕ, p_E, layout = (2, 4), size = (2000, 1000))
    png(p, "last")
    return p
end

function animate_solution_OVS(sol, z = nothing)
    @gif for (u, t) in zip(sol.u, sol.t)
        plot_solution_OVS(u, z)
    end
end


function animate_solution_all(sol, case, z = nothing)
    @gif for i in 1:length(sol.u)
        plot_solution(sol.u[i], sol.savevals[i], sol.params.z_cell, case)
    end
end

function write_sol_csv(filename, sol)
    CSV.write(filename*".csv", DataFrame(sol), header = false)
end

function write_sol_jld2(filename, sol)
    jldsave(filename*".jld2"; sol)
end

function read_csv(filename)
    sol = CSV.read(filename, DataFrame, header = false)
    return sol
end

function read_jld2(filename)
    f = jldopen(filename*".jld2", "r")
    sol = read(f, "sol")
    return sol
end

#=analyse simulation
make time averaged quantities and plot them, to compare with Landmark
make time resolved discharge current plots
make x, t space plots
include ionization rate, can be inferred after from nn, ne, and Tev
get boundaries right with not fixing value ion density
run simulation under three test cases
return to implicit and other general framework, see what happens
=#

function plot_current(current, sol)

    min_I, max_I = extrema(current)
    mid_I = (min_I + max_I)/2
    range = min(30, 5 + max_I - min_I)
    ylims = (mid_I - range/2, mid_I + range/2)
    p1 = plot(;ylims)
    plot!(p1, sol.t, current[1, :], title = "Currents at right boundary", label = ["Iᵢ" ""])
    plot!(p1, sol.t, current[2, :], label = ["Iₑ" ""])
    plot!(p1, sol.t, current[3, :], label = ["I total" ""])
    png(p1, "currents")
    return p1
end

function load_hallis_for_input()
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_01.out")
    ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, hallis.ϕ)
    grad_ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, -hallis.Ez)
    return ϕ_hallis, grad_ϕ_hallis
end

Plots.plot(sol::HallThruster.HallThrusterSolution, frame = length(sol.savevals); case) = plot_solution(sol.u[frame], sol.savevals[frame], sol.params.z_cell, case)

function plot_timeaveraged(sol, case, start_ind)
    avg, avg_savevals = HallThruster.timeaveraged(sol, start_ind)
    plot_solution(avg, avg_savevals, sol.params.z_cell, case)
end