using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra

function source!(Q, U, params, ϕ, Tev, i)
    HallThruster.apply_reactions!(Q, U, params, Tev, i)
    HallThruster.apply_ion_acceleration!(Q, U, params, ϕ, i)
    return Q
end

function source_potential!(b, i, i_f, μ⁻, μ⁺, pe, U, Δz)
    HallThruster.potential_source_term!(b, i, i_f, μ⁻, μ⁺, pe, U, Δz)
    #HallThruster.OVS_potential_source_term!(b, i)
end

function boundary_potential!(U, fluid, N, ϕ, pe, ne, B, A, b, Tev, νan, Δz, OVS) #if OVS, this sets to true
    ϕ_L = 300.0
    ϕ_R = 0.0
    HallThruster.boundary_conditions_potential!(U, fluid, N, pe, ne, B, A, b, Tev, νan, ϕ, ϕ_L, ϕ_R, Δz)
    #HallThruster.OVS_boundary_conditions_potential!(N, A, b, ϕ, ϕ_L, ϕ_R, Δz, OVS)
end

function IC!(U, z, fluids, L)
    ρ2 = 2.1801715574645586e-7 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0
    ρ1 = 5e-6/HallThruster.Xenon.m/0.004/u1
    @show ρ1
    U[1] = ρ1
    U[2] = ρ2
    U .= SA[ρ1, ρ2, ρ2*u1] #[ρ1, ρ1*u1, ρ1*E]
    return U
end

function IC_E!(E, U, z, L, fluid_ranges, fluids, i)
    Tev = 30 * exp(-(2 * (z - L/2) / 0.033)^2)
    #println("Tev in a row: ", Tev)
    ne = HallThruster.electron_density(U, fluid_ranges) / fluids[1].species.element.m
    #println("electron density in a row: ", ne)
    E[i] = 3/2*ne*Tev
end

function run_sim(end_time = 0.0002, n_save = 2)
    fluid = HallThruster.Xenon
    timestep = 0.1e-9 #0.9e-8

    ρ2 = 2.1801715574645586e-7 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0
    ρ1 = 5e-6/HallThruster.Xenon.m/0.004/u1
    T1 = 1000.0

    left_state = [ρ1, ρ2, ρ2 * u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2 * (u1 + 0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())
    saveat = if n_save == 1
        [end_time]
    else
        LinRange(0.0, end_time, n_save) |> collect
    end
    saved_values = SavedValues(Float64, NTuple{3, Vector{Float64}})
    callback = SavingCallback((U, tspan, integrator)->(integrator.p.cache.ϕ, integrator.p.cache.Tev, integrator.p.cache.ne), saved_values, saveat = saveat)

    sim = HallThruster.MultiFluidSimulation(
        grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = BCs,
        scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, identity, false),
        initial_condition = IC!,
        initial_condition_E = IC_E!,
        source_term! = source!,
        source_potential! = source_potential!,
        boundary_potential! = boundary_potential!,
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(u1, 300.0))
            HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time, #0.0002
        saveat = saveat, 
        timestepcontrol = (timestep, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = callback
    )

    @time sol = HallThruster.run_simulation(sim)

    p = plot() #plot(sol.u[end][1, :], yaxis = :log)
    plot!(p, sol.u[end][3, :] ./ sol.u[end][2, :])
    #plot!(p, sol.u[end][1, :]/HallThruster.Xenon.m)

    display(p)
    return sol, saved_values.saveval
end

function animate_solution(sol, saved_values)
    ϕ = Array{Union{Nothing, Vector{Float64}}}(nothing, length(sol.t))
    Tev = Array{Union{Nothing, Vector{Float64}}}(nothing, length(sol.t))
    ne = Array{Union{Nothing, Vector{Float64}}}(nothing, length(sol.t))
    for i in 1:length(sol.t)
        ϕ[i] = saved_values[i][1]
        Tev[i] = saved_values[i][2]
        ne[i] = saved_values[i][3]
    end
    mi = HallThruster.Xenon.m
    @gif for (u, t) in zip(sol.u, sol.t)
        p = plot(ylims = (1e13, 1e20))
        plot!(p, u[1, :] / mi, yaxis = :log)
        plot!(p, u[2, :] / mi)
    end
    @gif for (u, t) in zip(sol.u, sol.t)
        p = plot(ylims = (0, 3e4))
        plot!(p, u[3, :] ./ u[2, :])
    end
    @gif for (u, t) in zip(sol.u, sol.t)
        p = plot(ylims = (1e17, 1e21))
        plot!(p, u[4, :])
    end
    @gif for (ϕ, t) in zip(ϕ, sol.t)
        p = plot(ylims = (-100, 400))
        plot!(p, ϕ)
    end
    @gif for (Tev, t) in zip(Tev, sol.t)
        p = plot(ylims = (0, 30))
        plot!(p, Tev)
    end
    @gif for (ne, t) in zip(ne, sol.t)
        p = plot(ylims = (1e13, 1e20))
        plot!(p, ne)
    end
end