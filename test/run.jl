using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra

function source!(Q, U, params, i)
    HallThruster.apply_reactions!(Q, U, params, i)
    HallThruster.apply_ion_acceleration!(Q, U, params, i)
    #HallThruster.source_electron_energy!(Q, U, params, i)
    return Q
end

function source_potential!(b, U, s_consts)
    HallThruster.potential_source_term!(b, U, s_consts)
    #HallThruster.OVS_potential_source_term!(b, s_consts)
end

function boundary_potential!(A, b, U, bc_consts)
    ϕ_L = 300.0
    ϕ_R = 0.0
    HallThruster.boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
    #HallThruster.OVS_boundary_conditions_potential!((A, b, U, bc_consts, ϕ_L, ϕ_R)
end

function IC!(U, z, fluids, L)
    ρ2 = 2.1801715574645586e-7 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0
    ρ1 = 5e-6/0.004/u1
    Tev = 30 * exp(-(2 * (z - L/2) / 0.033)^2)
    ne = 2.1801715574645586e-7 / fluids[1].species.element.m
    U .= SA[ρ1, ρ2, ρ2*u1, 3/2*ne*Tev]
    return U
end

function run_sim(end_time = 0.0002, n_save = 2)
    fluid = HallThruster.Xenon
    timestep = 0.9e-8 #0.9e-8

    #fluid BCs #############################
    ρ2 = 2.1801715574645586e-7 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0 #150
    ρ1 = 5e-6/0.004/u1
    @show ρ1/HallThruster.Xenon.m
    left_state = [ρ1, ρ2, ρ2 * u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2 * (u1 + 0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    saveat = if n_save == 1
        [end_time]
    else
        LinRange(0.0, end_time, n_save) |> collect
    end
    #saved_values = SavedValues(Float64, NTuple{3, Vector{Float64}})
    callback = nothing #SavingCallback((U, tspan, integrator)->(integrator.p.cache.ϕ, integrator.p.cache.Tev, integrator.p.cache.ne), saved_values, saveat = saveat)

    sim = HallThruster.MultiFluidSimulation(
        grid = HallThruster.generate_grid(HallThruster.SPT_100, 10),
        boundary_conditions = BCs,
        scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, identity, false),
        initial_condition = IC!,
        source_term! = source!,
        source_potential! = source_potential!,
        boundary_potential! = boundary_potential!,
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(u1, 300.0))
            HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time, #0.0002
        saveat = saveat, 
        timestepcontrol = (timestep, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = callback,
        solve_energy = true
    )

    @time sol = HallThruster.run_simulation(sim)

    p = plot() #plot(sol.u[end][1, :], yaxis = :log)
    plot!(p, sol.u[end][3, :] ./ sol.u[end][2, :])
    #plot!(p, sol.u[end][1, :]/HallThruster.Xenon.m)

    display(p)
    return sol #, saved_values.saveval
end

function animate_solution(sol)
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
    @gif for (u, t) in zip(sol.u, sol.t) #nϵ
        p = plot(ylims = (1e13, 1e21))
        plot!(p, u[4, :], yaxis = :log)
    end
    @gif for (u, t) in zip(sol.u, sol.t) #Tev
        p = plot(ylims = (0, 40))
        plot!(p, u[5, :])
    end
    @gif for (u, t) in zip(sol.u, sol.t) #ne
        p = plot(ylims = (1e16, 1e20))
        plot!(p, u[6, :], yaxis = :log)
        plot!(p, u[7, :] ./ HallThruster.e)
    end
    @gif for (u, t) in zip(sol.u, sol.t) #pe
        p = plot(ylims = (1e13, 1e22))
        plot!(p, u[7, :] ./ HallThruster.e, yaxis = :log)
    end
    @gif for (u, t) in zip(sol.u, sol.t) #ϕ
        p = plot(ylims = (-100, 400))
        plot!(p, u[8, :])
    end
end
