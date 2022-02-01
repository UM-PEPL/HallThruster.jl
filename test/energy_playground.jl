#solve electron issues
#1) get rid of electron source terms, no heat conduction, and put the fluxes in explicit solve. Set electron velocity to a constant. See if this does indeed run and plot solution over time, set electron solve to false. 
#2) if above works, do MMS with that. 


using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase

function source!(Q, U, params, i)
    HallThruster.apply_reactions!(Q, U, params, i)
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
    ρ1 = 5e-6/0.004/abs(u1)
    u1 = -1000.0 + 80000*z
    #Tev = 200000 * exp(-(2 * (z - L/2) / 0.033)^2)
    Tev = 52220.331546975365
    ne = 2.1801715574645586e-7 / fluids[1].species.element.m
    U .= SA[ρ1, ρ2, ρ2*u1, 3/2*ne*Tev*HallThruster.kB]
    return U
end

function run_sim_energytest(end_time = 0.0002, n_save = 2)
    fluid = HallThruster.Xenon
    timestep = 0.9e-8

    #fluid BCs #############################
    ρ2 = 2.1801715574645586e-7 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0
    ρ1 = 5e-6/0.004/abs(u1)
    @show ρ1/HallThruster.Xenon.m
    left_state = [ρ1, ρ2, ρ2 * u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1*2, ρ2, ρ2 * (u1 + 0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet_ionbohm(left_state), HallThruster.Neumann_ionbohm())

    Tev = 52220.331546975365
    left_state_elec = 3/2*ρ2/HallThruster.Xenon.m*Tev*HallThruster.kB
    left_internal_energy = 3/2*Tev*HallThruster.kB
    right_state_elec = left_state_elec
    BCs_elec = (HallThruster.Dirichlet_energy_upd_ne(left_internal_energy), HallThruster.Dirichlet_energy_upd_ne(left_internal_energy))

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
            U[index.ne, i] = max(1e-10, HallThruster.electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
            if params.solve_energy
                U[index.Tev, i] = max(1, U[index.nϵ, i]/3*2/U[index.ne, i]/HallThruster.kB)
            end
            #U[index.ne, i] = max(1e-10, HallThruster.electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
            #U[index.pe, i] = HallThruster.electron_pressure(U[index.ne, i], U[index.Tev, i]) #this would be real electron pressure, ie next step use for previous in energy convection update
            U[index.pe, i] = U[index.nϵ, i]/3*2 #if using the same for pe and ne, might solve some instabilities
            U[index.grad_ϕ, i] = HallThruster.first_deriv_central_diff_pot(U[index.ϕ, :], params.z_cell, i)
            U[index.ue, i] = min(-150.0, HallThruster.electron_velocity(U, params, i))
            #@show U[index.ue, i]
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
        boundary_conditions = boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
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
