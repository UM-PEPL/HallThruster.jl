using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase

function source!(Q, U, params, i)
    #HallThruster.apply_reactions!(Q, U, params, i)
    HallThruster.apply_ion_acceleration!(Q, U, params, i)
    HallThruster.source_electron_energy!(Q, U, params, i)
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
    u1 = 150
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

    condition(u,t,integrator) = t < 1
    function affect!(integrator)
        U, params = integrator.u, integrator.p
        
        fluids, fluid_ranges = params.fluids, params.fluid_ranges
        index = params.index

        B = params.cache.B

        z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
        ncells = size(U, 2) - 2

        ####################################################################
        #PREPROCESS
        #calculate useful quantities relevant for potential, electron energy and fluid solve
        L_ch = 0.025
        fluid = fluids[1].species.element

        @inbounds for i in 1:(ncells + 2)
            #update electron temperature from energy using old density
            if params.solve_energy
                U[index.Tev, i] = max(1, U[index.nϵ, i]/3*2/U[index.ne, i])
            end
            U[index.ne, i] = max(1e-10, HallThruster.electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
            #U[index.pe, i] = HallThruster.electron_pressure(U[index.ne, i], U[index.Tev, i]) #this would be real electron pressure, ie next step use for previous in energy convection update
            U[index.pe, i] = U[index.nϵ, i]/3*2*HallThruster.e #if using the same for pe and ne, might solve some instabilities
            U[index.grad_ϕ, i] = HallThruster.first_deriv_central_diff(U[index.ϕ, :], params.z_cell, i)
            U[index.ue, i] = HallThruster.electron_velocity(U, params, i)
            params.cache.νan[i] = HallThruster.get_v_an(z_cell[i], B[i], L_ch)
            params.cache.νc[i] = HallThruster.get_v_c(U[index.Tev, i], U[1, i]/fluid.m , U[index.ne, i], fluid.m)
            params.cache.μ[i] = HallThruster.cf_electron_transport(params.cache.νan[i], params.cache.νc[i], B[i])
        end
        
        #POTENTIAL #########################################################
        HallThruster.solve_potential!(U, params)
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(false,false))

    callback = cb #SavingCallback((U, tspan, integrator)->(integrator.p.cache.ϕ, integrator.p.cache.Tev, integrator.p.cache.ne), saved_values, saveat = saveat)

    sim = HallThruster.MultiFluidSimulation(
        grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
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

function animate_solution1(sol)
    mi = HallThruster.Xenon.m
    @gif for (u, t) in zip(sol.u, sol.t)
        p1 = plot(ylims = (1e17, 1e20))
        plot!(p1, u[1, :] / mi, yaxis = :log, title = "Neutral and ion densities [n/m^3]", label = ["nₙ" ""])
        plot!(p1, u[2, :] / mi, label = ["nᵢ" ""])
        p2 = plot(ylims = (-1300, 1300))
        plot!(p2, u[3, :] ./ u[2, :], title = "Ion velocity [m/s]", label = ["vᵢ" ""])
        p3 = plot(ylims = (0, 5))
        plot!(p3, u[4, :], title = "Internal electron energy [eV*n/m^3]", label = ["nϵ" ""])
        p4 = plot(ylims = (0, 120000))
        plot!(p4, u[5, :], title = "Electron temperature [eV]", label = ["Tev" ""])
        p5 = plot(ylims = (1e16, 1e20))
        plot!(p5, u[6, :], yaxis = :log, title = "Electron density and pressure", label = ["nₑ [n/m^3]" ""])
        plot!(p5, u[7, :] ./ HallThruster.e, label = ["pₑ [n*eV/m^3]" ""])
        p6 = plot(ylims = (-1e5, 1e5))
        plot!(p6, u[10, :], title = "Electron velocity", label = ["uₑ [m/s]" ""])
        p7 = plot(ylims = (-100, 400))
        plot!(p7, u[8, :], title = "Potential", label = ["ϕ [V]" ""])
        p8 = plot(ylims = (-80000, 1000))
        plot!(p8, u[9, :], title = "Electric field", label = ["E [V/m]" ""])

        plot!(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4), size = (2000, 1000))
    end
end

function animate_solution2(sol)
    mi = HallThruster.Xenon.m
    @gif for (u, t) in zip(sol.u, sol.t)
        p1 = plot(ylims = (0, 5000))
        plot!(p1, u[1, :], title = "Neutral and ion densities [n/m^3]", label = ["nₙ" ""])
        plot!(p1, u[2, :], label = ["nᵢ" ""])
        p2 = plot(ylims = (200, 500))
        plot!(p2, u[3, :] ./ u[2, :], title = "Ion velocity [m/s]", label = ["vᵢ" ""])
        p3 = plot(ylims = (1e13, 1e25))
        plot!(p3, u[4, :], yaxis = :log, title = "Internal electron energy [eV*n/m^3]", label = ["nϵ" ""])
        p4 = plot(ylims = (0, 40))
        plot!(p4, u[5, :], title = "Electron temperature [eV]", label = ["Tev" ""])
        p5 = plot(ylims = (1e16, 1e20))
        plot!(p5, u[6, :], yaxis = :log, title = "Electron density and pressure", label = ["nₑ [n/m^3]" ""])
        plot!(p5, u[7, :] ./ HallThruster.e, label = ["pₑ [n*eV/m^3]" ""])
        p6 = plot(ylims = (-1e5, 1e5))
        plot!(p6, u[10, :], title = "Electron velocity", label = ["uₑ [m/s]" ""])
        p7 = plot(ylims = (-100, 400))
        plot!(p7, u[8, :], title = "Potential", label = ["ϕ [V]" ""])
        p8 = plot(ylims = (-80000, 1000))
        plot!(p8, u[9, :], title = "Electric field", label = ["E [V/m]" ""])

        plot!(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4), size = (2000, 1000))
    end
end