using OrdinaryDiffEq, UnicodePlots

function update_energy!(du, u, params, t)
    ue, ne, nn, E, μ, ϕ = params.ue, params.ne, params.nn, params.E, params.μ, params.ϕ
    ε = params.ε
    z = params.zs
    dz = z[2] - z[1]
    ncells = length(z)
    K = params.K
    νε = 1e7
    U = 20.0

    for i in 2:ncells-1
        E[i] = -(ϕ[i+1] - ϕ[i-1])/2/dz
        ue[i] = μ[i] * (-E[i] - (u[1, i+1] - u[1, i-1])/(2dz * ne[i]))
        ε[i] = u[1, i] / ne[i]
        u[2, i] = ue[i]
    end

    for i in 2:ncells-1

        advection_term = -(ue[i] * u[1, i] - ue[i+1] * u[1, i+1]) / dz

        diffusion_term = (μ[i+1] * u[1, i+1] - μ[i] * u[1, i]) * (ε[i+1] - ε[i])
        diffusion_term += μ[i] * u[1, i] * (ε[i-1] - 2ε[i] + ε[i+1])
        diffusion_term /= dz^2

        W = νε * ε[i] * exp(-U / ε[i])
        source_term = ne[i] * (-ue[i] * E[i] - nn[i] * K(ε[i]) - W)
        du[1, i] = -5/3 * advection_term + 10/9 * diffusion_term + source_term
    end
end

function solve_energy(;ncells = 50, tmax = 1e-3, dt = 1e-9)
    d = 0.05
    V = 300.0

    e = 1.6e-19
    me = 9.1e-31

    zs = LinRange(0.0, d, ncells)
    Bmax = 0.015
    l = 0.025
    δB(z) = z ≤ l ? 0.011 : 0.018
    B = [Bmax * exp(-(z - l)^2 / 2δB(z)) for z in zs]

    ne_f(z) = 1.2 * exp(-(z - 0.014)^2 / (0.01)^2)
    ne = [1e18 * (z > l ? (l - z) + ne_f(l) : ne_f(z)) for z in zs]
    nn = [4e19 * (1 - 1 / (1 + exp(-(z - 0.01) / 0.002))) for z in zs]
    ϕ = [V * (1 - 1 / (1 + exp(-(z - l) / 0.005))) for z in zs]
    ϕ[1] = V
    ϕ[end] = 0.0

    ν = [
        2.5e-13 * nn[i] +
        (
            z < l ?
                1e7 + 0.1/16 * e * B[i] / me :
                1/16 * e * B[i]/ me
        )
        for (i, z) in enumerate(zs)
    ]
    μ = e / me * (ν ./ (ν.^2 + (e * B / me).^2))

    # this is a fit function to the LANDMARK data
    K(ε) = 6e-12 * exp(-39.0 / (ε + 3.0))

    dz = zs[2] - zs[1]
    E = zeros(ncells)
    for i in 2:ncells-1
        E[i] = -(ϕ[i+1] - ϕ[i-1])/(2dz)
    end
    E[1] = (ϕ[2] - ϕ[1])/dz
    E[2] = (ϕ[end] - ϕ[end-1])/dz

    ue = zeros(ncells)
    ε = 3.0 .* ones(ncells)

    for i in 2:ncells-1
        ue[i] = μ[i] * (-E[i] - (ne[i+1] * ε[i+1] - ne[i-1] * ε[i-1])/2dz/ne[i])
    end

    u0 = vcat((ne .* ε)', ue')

    params = (;
        zs, μ, ϕ, nn, ne, B, K, E, ue, ε, ν
    )

    tmax = 1e-4
    tspan = (0.0, tmax)
    prob = ODEProblem{true}(update_energy!, u0, tspan, params)
    @time sol = solve(prob, SSPRK22(), dt = dt, saveat = 0.0:1e-6:tmax)
    lineplot(zs ./ 100, sol.u[end][1,:] ./ params.ne, xlabel = "z (cm)", ylabel = "ε (eV)") |> display
    return sol, params
end

