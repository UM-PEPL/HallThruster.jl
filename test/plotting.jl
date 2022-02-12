using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase, DataFrames, CSV, JLD2


function animate_solution_individual(sol)
    mi = HallThruster.Xenon.m
    @gif for (u, t) in zip(sol.u, sol.t)
        p = plot(ylims = (1e13, 1e20))
        plot!(p, u[1, :] / mi, yaxis = :log, title = "Neutral and ion densities [n/m^3]", label = ["nₙ" ""])
        plot!(p, u[2, :] / mi, label = ["nᵢ" ""])
    end
    @gif for (u, t) in zip(sol.u, sol.t)
        p = plot(ylims = (0, 3e4))
        plot!(p, u[3, :] ./ u[2, :], title = "Ion velocity [m/s]", label = ["vᵢ" ""])
    end
    @gif for (u, t) in zip(sol.u, sol.t) #nϵ
        p = plot(ylims = (0, 20))
        plot!(p, u[4, :], title = "Internal electron energy [eV*n/m^3]", label = ["nϵ" ""])
    end
    @gif for (u, t) in zip(sol.u, sol.t) #Tev
        p = plot(ylims = (0, 120000))
        plot!(p, u[5, :], title = "Electron temperature [eV]", label = ["Tev" ""])
    end
    @gif for (u, t) in zip(sol.u, sol.t) #ne
        p = plot(ylims = (1e16, 1e20))
        plot!(p, u[6, :], yaxis = :log, title = "Electron density and pressure", label = ["nₑ [n/m^3]" ""])
        plot!(p, u[7, :] ./ HallThruster.e, label = ["pₑ [n*eV/m^3]" ""])
    end
    @gif for (u, t) in zip(sol.u, sol.t) #pe
        p = plot(ylims = (1e13, 1e22))
        plot!(p, u[7, :] ./ HallThruster.e, yaxis = :log, title = "Electron pressure", label = ["pₑ [n*eV/m^3]" ""])
    end
    @gif for (u, t) in zip(sol.u, sol.t) #ϕ
        p = plot(ylims = (-100, 400))
        plot!(p, u[8, :], title = "Potential", label = ["ϕ [V]" ""])
    end
end

using DelimitedFiles

function plot_quantity(u, zmin = 0.0, zmax = 0.05; ref_path = nothing, hallis = nothing, hallisvar = nothing, kwargs...)
    z = LinRange(zmin, zmax, length(u))
    p = plot()
    plot!(
        p, z, u; label = "HallThruster.jl", xlabel = "z (m)", legend = :outertop, bottommargin = 7Plots.mm, topmargin = 7Plots.mm, lw = 2,
        kwargs...
    )
    if !isnothing(ref_path)
        ref_data = readdlm(ref_path, ',')
        z_ref, q_ref = ref_data[:, 1], ref_data[:, 2]
        plot!(p, z_ref, q_ref, label = "Reference", lw = 2, lc = :red, ls = :dash)
    elseif !isnothing(hallis)
        z_ref, q_ref = hallis.z_cm, hallisvar
        plot!(p, z_ref, q_ref, label = "Reference", lw = 2, lc = :red, ls = :dash)
    end
    return p
end

function load_hallis_output(output_path)
    output_headers = [
        :z_cm, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei,
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    return output[1:end-1, :]
end

function plot_solution(u, case = 1)
    mi = HallThruster.Xenon.m
    p_nn = plot_quantity(u[1, :] / mi, title = "Neutral density", ylabel = "nn (m⁻³)", ref_path = "landmark/landmark_neutral_density_$(case).csv")
    p_ne = plot_quantity(u[2, :] / mi, title = "Plasma density", ylabel = "ne (m⁻³)", ref_path = "landmark/landmark_plasma_density_$(case).csv")
    p_ui = plot_quantity(u[3, :] ./ u[2, :] ./ 1000, title = "Ion velocity", ylabel = "ui (km/s)")
    p_nϵ = plot_quantity(u[4, :], yaxis = :log, title = "Energy density", ylabel = "nϵ (eV m⁻³)")
    p_ϵ  = plot_quantity(u[5, :], title = "Electron temperature (eV)", ylabel = "ϵ (eV)", ref_path = "landmark/landmark_electron_temperature_$(case).csv")
    p_ue = plot_quantity(u[10, :] ./ 1000, title = "Electron velocity", ylabel = "ue (km/s)")
    p_ϕ  = plot_quantity(u[8, :], title = "Potential", ylabel = "ϕ (V)", ref_path = "landmark/landmark_potential_$(case).csv")
    p_E  = plot_quantity(-u[9, :], title = "Electric field", ylabel = "E (V/m)", ref_path = "landmark/landmark_electric_field_$(case).csv")
    plot(p_nn, p_ne, p_ui, p_nϵ, p_ϵ, p_ue, p_ϕ, p_E, layout = (2, 4), size = (2000, 1000))
end

function plot_solution_real(u, case = 1)
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_0$(case).out")
    coeff = HallThruster.load_landmark()
    mi = HallThruster.Xenon.m
    for i in 1:length(u[4, :])
        u[4, i] = coeff.rate_coeff(u[5, i])*u[1, i]*u[2, i]/mi/mi
    end
    p_nn = plot_quantity(u[1, :] / mi, title = "Neutral density", ylabel = "nn (m⁻³)", hallis = hallis, hallisvar = hallis.nn)
    p_ne = plot_quantity(u[2, :] / mi, title = "Plasma density", ylabel = "ne (m⁻³)", hallis = hallis, hallisvar = hallis.ne)
    p_ui = plot_quantity(u[3, :] ./ u[2, :] ./ 1000, title = "Ion velocity", ylabel = "ui (km/s)")
    p_nϵ = plot_quantity(u[4, :], title = "Ionization_rate", ylabel = "nϵ (eV m⁻³)", hallis = hallis, hallisvar = hallis.ndot)
    p_ϵ  = plot_quantity(u[5, :], title = "Electron temperature (eV)", ylabel = "ϵ (eV)", hallis = hallis, hallisvar = hallis.Te)
    p_ue = plot_quantity(u[10, :] ./ 1000, title = "Electron velocity", ylabel = "ue (km/s)")
    p_ϕ  = plot_quantity(u[8, :], title = "Potential", ylabel = "ϕ (V)", hallis = hallis, hallisvar = hallis.ϕ)
    p_E  = plot_quantity(-u[9, :], title = "Electric field", ylabel = "E (V/m)", hallis = hallis, hallisvar = hallis.Ez)
    plot(p_nn, p_ne, p_ui, p_nϵ, p_ϵ, p_ue, p_ϕ, p_E, layout = (2, 4), size = (2000, 1000))
end

function animate_solution_all(sol)
    @gif for (u, t) in zip(sol.u, sol.t)
        plot_solution(u)
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

function timeaveraged(sol, tstampstart)
    avg = zeros(size(sol.u[1]))
    tstamps = length(sol.t)
    for i in tstampstart:length(sol.t)
        avg .+= sol.u[i]
    end
    avg /= (tstamps - tstampstart + 1)
    return avg
end

function calc_current(sol) #need to divide by area?
    current = zeros(2, length(sol.t))
    area = pi*(0.05^2 - 0.025^2)
    distance = 0.050 - 0.0350
    for i in 1:length(sol.t)
        current[1, i] = sol.u[i][3, end]*HallThruster.e/HallThruster.Xenon.m*area
        current[2, i] = -sol.u[i][6, end]*HallThruster.e*sol.u[i][10, end]*area
    end
    return current
end

function plot_current(current, sol)
    p1 = plot(ylims = (0, 30))
    plot!(p1, sol.t, current[1, :], title = "Currents at right boundary", label = ["Iᵢ" ""])
    plot!(p1, sol.t, current[2, :], label = ["Iₑ" ""])
    plot!(p1, sol.t, current[2, :] + current[1, :], label = ["I total" ""])
    png(p1, "currents")
    return p1
end