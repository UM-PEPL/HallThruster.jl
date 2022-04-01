using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase, DataFrames, CSV, JLD2, FFTW

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


function plot_quantity(u, z = nothing, zmin = 0.0, zmax = 0.05; normalize_z_factor = 1.0, ref_paths = String[], ref_styles = nothing, hallis = nothing, hallisvar = nothing, slicezr = nothing, slicezrvar = nothing, slice = nothing, slicevar = nothing, label = "HallThruster.jl", kwargs...)
    if isnothing(z)
        z = LinRange(zmin, zmax, length(u))
    end

    z ./= normalize_z_factor

    p = plot(; xlims = (z[1], z[end]), kwargs...)

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
        z_ref, q_ref = hallis.z / normalize_z_factor, hallisvar
        plot!(p, z_ref, q_ref, label = "Hallis 2D", lw = 2, lc = :red, ls = :dash)
    end
    if !isnothing(slicezr)
        z_ref, q_ref = slicezr.z / normalize_z_factor, slicezrvar
        plot!(p, z_ref, q_ref, label = "Slicezr", lw = 2, lc = :green, ls = :dash)
    end
    if !isnothing(slice)
        z_ref, q_ref = slice.z / normalize_z_factor, slicevar
        plot!(p, z_ref, q_ref, label = "Hall2De", lw = 2, lc = :blue, ls = :dashdot)
    end

    plot!(
        p, z, u; label = label, legend = :outertop, margin = 8Plots.mm, lw = 2,
        color = :black, linestyle = :solid,
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

function load_hallis_output_2D(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei, :f_e, :f_en, :f_bohm, :f_wall, :f_ei
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    #output.Te .*= 2/3
    return output[1:end-1, :]
end

function plot_solution(u, saved_values, z, case = 1)
    mi = HallThruster.Xenon.m
    Xe_0 = HallThruster.Xenon(0)
    Xe_I = HallThruster.Xenon(1)
    rxn = HallThruster.load_ionization_reactions(HallThruster.LandmarkIonizationLUT(), [Xe_0, Xe_I])[1]
    (;Tev, ue, ϕ_cell, ∇ϕ, ne, pe, ∇pe) = saved_values
    ionization_rate = [rxn.rate_coeff(3/2 * Tev[i])*u[1, i]*ne[i]/mi for i in 1:size(u, 2)]

    ref_styles = landmark_styles()

    p_nn = plot_quantity(
        u[1, :] / mi, z; title = "Neutral density",  ylabel = "nn (m⁻³)",
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
        u[end, :] ./ ne, z; title = "Electron energy (3/2 Te) (eV)", ylabel = "ϵ (eV)",
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

function plot_solution_hallis_and_2D(u, saved_values, z, case = 1, hallis = nothing, slicezr = nothing, slice = nothing)
    mi = HallThruster.Xenon.m
    Xe_0 = HallThruster.Xenon(0)
    Xe_I = HallThruster.Xenon(1)
    rxn = HallThruster.load_ionization_reactions(HallThruster.LandmarkIonizationLUT(), [Xe_0, Xe_I])[1]
    (;Tev, ue, ϕ_cell, ∇ϕ, ne, pe, ∇pe) = saved_values
    ionization_rate = [rxn.rate_coeff(3/2 * Tev[i])*u[1, i]*ne[i]/mi for i in 1:size(u, 2)]

    ref_styles = landmark_styles()

    p_nn = plot_quantity(
        u[1, :] / mi, z; title = "Neutral density",  ylabel = "nn (m⁻³)",
        ref_paths = landmark_references(case, "neutral_density"), ref_styles, hallis = hallis, hallisvar = hallis.nn,
        slice = slice, slicevar = slice.nn
    )

    p_ne = plot_quantity(
        ne, z; title = "Plasma density", ylabel = "ne (m⁻³)",
        ref_paths = landmark_references(case, "plasma_density"), ref_styles, hallis = hallis, hallisvar = hallis.ne,
        slice = slice, slicevar = slice.ne
    )

    p_ui = plot_quantity(u[3, :] ./ u[2, :] ./ 1000, z; title = "Ion velocity", ylabel = "ui (km/s)", slice = slice, slicevar = slice.uiz_1_1)

    p_iz = plot_quantity(
        ionization_rate, z; title = "Ionization rate", ylabel = "nϵ (eV m⁻³)",
        ref_paths = landmark_references(case, "ionization"), ref_styles, hallis = hallis, hallisvar = hallis.ndot
    )

    p_ϵ  = plot_quantity(
        u[4, :] ./ ne, z; title = "Electron energy (3/2 Te) (eV)", ylabel = "ϵ (eV)",
        ref_paths = landmark_references(case, "energy"), ref_styles, hallis = hallis, hallisvar = hallis.Te,
        slice = slice, slicevar = slice.Te
    )

    p_ue = plot_quantity(ue ./ 1000, z; title = "Electron velocity", ylabel = "ue (km/s)")
    p_ϕ  = plot_quantity(
        ϕ_cell, z; title = "Potential", ylabel = "ϕ (V)",
        ref_paths = landmark_references(case, "potential"), ref_styles, hallis = hallis, hallisvar = hallis.ϕ,
        slice = slice, slicevar = slice.ϕ
    )

    p_E  = plot_quantity(
        -∇ϕ, z; title = "Electric field", ylabel = "E (V/m)",
        ref_paths = landmark_references(case, "electric_field"), ref_styles
    )

    #p_pe  = plot_quantity(HallThruster.e * pe, z; title = "Electron pressure", ylabel = "∇pe (Pa)")
    #p_∇pe  = plot_quantity(HallThruster.e * ∇pe, z; title = "Pressure gradient", ylabel = "∇pe (Pa/m)")
    plot(p_nn, p_ne, p_ui, p_ϕ, #=p_pe,=# p_iz, p_ϵ, p_ue, p_E, #=p_∇pe,=# layout = (2, 4), size = (2500, 1000))
end

function plot_multiple_solution(sols, labels, case = 1, timeaveraged_start = nothing)
    
    if timeaveraged_start !== nothing
        u, saved_values = HallThruster.timeaveraged(sols[1], timeaveraged_start)
        z = sols[1].params.z_cell
    else
        u, saved_values, z = sols[1].u[end], sols[1].savevals[end], sols[1].params.z_cell
    end
    mi = HallThruster.Xenon.m
    Xe_0 = HallThruster.Xenon(0)
    Xe_I = HallThruster.Xenon(1)
    rxn = HallThruster.load_ionization_reactions(HallThruster.LandmarkIonizationLUT(), [Xe_0, Xe_I])[1]
    (;Tev, ue, ϕ_cell, ∇ϕ, ne, pe, ∇pe) = saved_values
    ionization_rate = [rxn.rate_coeff(3/2 * Tev[i])*u[1, i]*ne[i]/mi for i in 1:size(u, 2)]

    ref_styles = landmark_styles()

    p_nn = plot_quantity(
        u[1, :] / mi, z; title = "Neutral density", #=yaxis = :log, ylims = (1e16, 1e21),=#  ylabel = "nn (m⁻³)",
        ref_paths = landmark_references(case, "neutral_density"), ref_styles, label = labels[1]
    )

    p_ne = plot_quantity(
        ne, z; title = "Plasma density", ylabel = "ne (m⁻³)", #=yaxis = :log, ylims = (1e16, 1e21),=#
        ref_paths = landmark_references(case, "plasma_density"), ref_styles, label = labels[1]
    )

    p_ui = plot_quantity(u[3, :] ./ u[2, :] ./ 1000, z; title = "Ion velocity", ylabel = "ui (km/s)", label = labels[1])

    p_iz = plot_quantity(
        ionization_rate, z; title = "Ionization rate", ylabel = "nϵ (eV m⁻³)",
        ref_paths = landmark_references(case, "ionization"), ref_styles, label = labels[1]
    )

    p_ϵ  = plot_quantity(
        u[4, :] ./ ne, z; title = "Electron energy (3/2 Te) (eV)", ylabel = "ϵ (eV)",
        ref_paths = landmark_references(case, "energy"), ref_styles, label = labels[1]
    )

    p_ue = plot_quantity(ue ./ 1000, z; title = "Electron velocity", ylabel = "ue (km/s)", label = labels[1])
    p_ϕ  = plot_quantity(
        ϕ_cell, z; title = "Potential", ylabel = "ϕ (V)",
        ref_paths = landmark_references(case, "potential"), ref_styles, label = labels[1]
    )

    p_E  = plot_quantity(
        -∇ϕ, z; title = "Electric field", ylabel = "E (V/m)",
        ref_paths = landmark_references(case, "electric_field"), ref_styles, label = labels[1]
    )

    #
    for i in 2:length(sols)
        if timeaveraged_start !== nothing
            u, saved_values = HallThruster.timeaveraged(sols[i], timeaveraged_start)
            z = sols[i].params.z_cell
        else
            u, saved_values, z = sols[i].u[end], sols[i].savevals[end], sols[i].params.z_cell
        end
        mi = HallThruster.Xenon.m
        (;Tev, ue, ϕ_cell, ∇ϕ, ne, pe, ∇pe) = saved_values
        ionization_rate = [rxn.rate_coeff(3/2 * Tev[i])*u[1, i]*ne[i]/mi for i in 1:size(u, 2)]
    
        plot!(p_nn, z, u[1, :] / mi, label = labels[i])
        plot!(p_ne, z, ne, label = labels[i])
        plot!(p_ui, z, u[3, :] ./ u[2, :] ./ 1000, label = labels[i])
        plot!(p_iz, z, ionization_rate, label = labels[i])
        plot!(p_ϵ, z, u[4, :] ./ ne, label = labels[i])
        plot!(p_ue, z, ue ./ 1000, label = labels[i])
        plot!(p_ϕ, z, ϕ_cell, label = labels[i])
        plot!(p_E, z, -∇ϕ, label = labels[i])
    end

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

    p1 = plot()
    plot!(p1, sol.t, current[1, :], title = "Currents at right boundary", label = ["Iᵢ" ""])
    plot!(p1, sol.t, current[2, :], label = ["Iₑ" ""])
    plot!(p1, sol.t, current[3, :], label = ["I total" ""])
    return p1
end

function plot_current_compare(currents, sols, labels)
    min_I, max_I = extrema(currents[1])
    mid_I = (min_I + max_I)/2
    range = min(30, 5 + max_I - min_I)
    ylims = (mid_I - range/2, mid_I + range/2)
    p1 = plot(;ylims, size = (1000, 500), xlabel = "t [s]", ylabel = "I [A]", margin = 5Plots.mm)
    #plot!(p1, sols[1].t, currents[1][1, :], linestyle = :solid, color = :green, title = "Total current vs grid size", label = labels[1] * "  Iᵢ")
    #plot!(p1, sols[1].t, currents[1][2, :], linestyle = :solid, color = :blue, label =  labels[1] * " Iₑ")
    plot!(p1, sols[1].t, currents[1][3, :], linestyle = :solid, color = :black, label =  labels[1] * " I total", title = "Total current vs grid size")
    linestyles = (:placeholder, :dashdot, :dot, :dash)
    colors = (:placeholder, :blue, :green, :red)
    for i in 2:length(currents)
        #plot!(p1, sols[i].t, currents[i][1, :], linestyle = linestyles[i], color = :green, label = labels[i] * "  Iᵢ")
        #plot!(p1, sols[i].t, currents[i][2, :], linestyle = linestyles[i], color = :blue, label =  labels[i] * " Iₑ")
        plot!(p1, sols[i].t, currents[i][3, :], linestyle = linestyles[i], color = colors[i], label =  labels[i] * " I total")
    end
    return p1
end
    

function load_hallis_for_input()
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_01.out")
    ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, hallis.ϕ)
    grad_ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, -hallis.Ez)
    return ϕ_hallis, grad_ϕ_hallis
end

Plots.plot(sol::HallThruster.HallThrusterSolution, frame = length(sol.savevals); case = 43) = plot_solution(sol.u[frame], sol.savevals[frame], sol.params.z_cell, case)

function plot_timeaveraged(sol, case, start_ind)
    avg, avg_savevals = HallThruster.timeaveraged(sol, start_ind)
    plot_solution(avg, avg_savevals, sol.params.z_cell, case)
end

function plot_thrust(thrust)
    plot(sol.t, thrust, title = "Thrust", label = "Thrust", xlabel = "t [s]", ylabel = "F [N]")
end

function plot_thrust_compare(thrusts, labels)
    p1 = plot(sol.t, thrust[1], title = "Thrust", label = thrusts[1], xlabel = "t [s]", ylabel = "F [N]")
    for i in 2:length(thrusts)
        plot!(p1, sol.t, thrusts[i], label = labels[i])
    end
end

function load_slicezr(slicezr_path)
    slicezr_headers = [
        :z, :ne, :nn, :Te, :ϕ, :Ωe, :ni_1_1, :ni_1_2, :ni_2_1, :uiz_1_1, :uPIC, :uiz_2_1
    ]

    slicezr = DataFrame(readdlm(slicezr_path, Float64), slicezr_headers)
    slicezr.z = slicezr.z./100
    slicezr.Te = slicezr.Te .* 3/2
    return slicezr
end

function load_slice(slice_path)
    slice_headers = [
        :z, :ne, :nn, :Te, :ϕ, :f_en, :f_ei, :f_wall, :f_iz, :Ωe, :ni_1_1, :ni_2_1, :ni_1_3, :uiz_1_1, :uir_1_1
    ]

    slice = DataFrame(readdlm(slice_path, Float64), slice_headers)
    slice.f_e = slice.f_en + slice.f_ei + slice.f_wall
    slice.ωce = slice.Ωe .* slice.f_e
    slice.z = slice.z/100
    slice.Te = slice.Te .* 3/2
    slice.uiz_1_1 = slice.uiz_1_1./1000
    return slice
end

function plot_PSD_current(current, sol)
    # Number of points 
    N = size(current)[2]
    # Sample period
    Ts = sol.t[2] - sol.t[1]
    #time
    t = sol.t

    # signal 
    signal = current[3, :]

    # Fourier Transform of it 
    F = fft(signal) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

    #return resolution limits
    #max determined by period of samples
    max_frequ = 1/2/Ts
    #min by length of aquisition
    min_frequ = 2/t[end]
    # plots
    time_domain = plot(t, signal, title = "Signal", label = "Total current", ylabel = "[A]", xlabel = "time [s]", margin = 5Plots.mm)
    freq_domain = plot(freqs, abs.(F), title = "Spectrum", ylim =(0, 2e4), xlim=(min_frequ, 20000), label = "Power spectral density", xlabel = "frequency [Hz]", ylabel = "[dB/Hz]") 
    plot(time_domain, freq_domain, layout = 2, size = (1000, 500))
end