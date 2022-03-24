using Statistics

function Lp_norm(v, p)
    N = length(v)
    return N^(-1/p) * norm(v, p)
end

function sin_wave(var; amplitude, phase, nwaves, offset = 0.0)
    return amplitude * sin(2π * nwaves * var + phase) + offset
end

function test_refinements(verification_func, refinements, norm_orders)
    norms = [
        let results = verification_func(ncells)
            [Lp_norm(res.sim .- res.exact, p) for res in results, p in norm_orders]
        end
        for ncells in refinements
    ] |> unzip

    slopes = compute_slope.(norms)
    return slopes, norms
end

function unzip(v)
    ncolumns = length(v[1])
    nrows = length(v)
    return [
        [v[j][i] for j in 1:nrows] for i in 1:ncolumns
    ]
end

function compute_slope(errors)
    p = [
        log(abs(errors[i+2]-errors[i+1])/abs(errors[i+1]-errors[i]))/log(0.5) for i in 1:length(errors)-2
    ]
    return mean(p)
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

function refines(num_refinements, initial, factor)
    return [
        initial * factor^p
        for p in 1:num_refinements
    ]
end