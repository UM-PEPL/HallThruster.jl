using Statistics, HallThruster

function Lp_norm(v, p)
    N = length(v)
    return N^(-1/p) * norm(v, p)
end

function sin_wave(var; amplitude, phase, nwaves, offset = 0.0)
    return amplitude * sin(nwaves * (2π * var) + phase) + offset
end

function test_refinements(verification_func, refinements, norm_orders)
    norms = [
        let results = verification_func(ncells)
            [Lp_norm(res.sim .- res.exact, p) for res in results, p in norm_orders]
        end
        for ncells in refinements
    ] |> unzip

    slopes = [expsmooth(compute_slope(refinements, norm), 0.75)[end] for norm in norms]
    return slopes, norms
end

function unzip(v)
    ncolumns = length(v[1])
    nrows = length(v)
    return [
        [v[j][i] for j in 1:nrows] for i in 1:ncolumns
    ]
end

function compute_slope(refinements, errors)
    q = [
        #log(abs(errors[i+2]-errors[i+1])/abs(errors[i+1]-errors[i]))/log(0.5) for i in 1:length(errors)-2
        log(errors[i+1] / errors[i]) /
        log(refinements[i] / refinements[i+1])
        for i in 1:length(errors)-1
    ]
    return q
end

function expsmooth(xs, α)
    smoothed = copy(xs)
    for i in 2:length(xs)
        smoothed[i] = α * xs[i] + (1 - α) * smoothed[i-1]
    end
    return smoothed
end

function plot_convergence(refinements, errors)
    normalized_errors = errors ./ errors[1]
    first_order_errors = [(1 / ncells) for ncells in refinements]
    first_order_errors ./= first_order_errors[1]
    second_order_errors = [(1 / ncells)^2 for ncells in refinements]
    second_order_errors ./= second_order_errors[1]

    p = plot(
        yaxis = (:log10, "Relative error"), xaxis = (:log10, "Number of cells"), legend = :bottomleft
    )
    plot!(p, refinements, first_order_errors, linestyle = :dash, label = "1st order")
    plot!(p, refinements, second_order_errors, linestyle = :dash, label = "2nd order")
    plot!(p, refinements, normalized_errors, label = "Observed error")
    return p
end

titles_ions = ("Neutral continuity", "Ion continuity", "Ion momentum", "Neutral continuity", "Ion continuity", "Ion momentum", "Neutral continuity", "Ion continuity", "Ion momentum")
titles_ϕ = ("", "", "")
titles_pot = ["Potential"]
titles_grad = ("∇ϕ", "∇pe", "ue", "∇ϕ", "∇pe", "ue", "∇ϕ", "∇pe", "ue")
titles_norms = ("L1", "L2", "LInf")
titles_energy = ["Energy Implicit"]

function plot_entire_range_convergence(refinements, norms, titles)
    N = length(norms)
    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, N)
    for i in 1:N
        p[i] = plot_convergence(refinements, norms[i])
        p[i] = plot!(p[i], title = titles[i])
    end
    return p
end

function plot_and_save_OVS_pot(prange, titles, saveloc, addition)
    p1 = plot(p[1], p[2], p[3], layout = (1, 3), size = (1200, 400), margin = 6Plots.mm, plot_title = titles[1])
    png(p1, saveloc*"ovs"*addition*".png")
end

function plot_and_save_OVSionneutral_grad(prange, titles_norms, saveloc, addition) #assuming L1, L2, LInf
    p1 = plot(p[1], p[2], p[3], layout = (1, 3), size = (1200, 400), margin = 6Plots.mm, plot_title = titles_norms[1])
    p2 = plot(p[4], p[5], p[6], layout = (1, 3), size = (1200, 400), margin = 6Plots.mm, plot_title = titles_norms[2])
    p3 = plot(p[7], p[8], p[9], layout = (1, 3), size = (1200, 400), margin = 6Plots.mm, plot_title = titles_norms[3])
    png(p1, saveloc*"ovs_ions_L1"*addition*".png")
    png(p2, saveloc*"ovs_ions_L2"*addition*".png")
    png(p3, saveloc*"ovs_ions_Linf"*addition*".png")
end


function refines(num_refinements, initial, factor)
    return [
        initial * factor^(p-1)
        for p in 1:num_refinements
    ]
end

struct OVS_Ionization <: HallThruster.IonizationModel end
struct OVS_Excitation <: HallThruster.ExcitationModel end

import HallThruster.load_reactions

function OVS_rate_coeff_iz(ϵ)
    return 1e-12 * exp(-12.12/ϵ)
end

function OVS_rate_coeff_ex(ϵ)
    return 1e-12 * exp(-8.32/ϵ)
end

function HallThruster.load_reactions(::OVS_Ionization, species)
    return [HallThruster.IonizationReaction(12.12, HallThruster.Xenon(0), HallThruster.Xenon(1), OVS_rate_coeff_iz)]
end

function HallThruster.load_reactions(::OVS_Excitation, species)
    return [HallThruster.ExcitationReaction(8.32, HallThruster.Xenon(0), OVS_rate_coeff_ex)]
end