function test_ion_accel_source(fluxfn, reconstruct, end_time, dt)

    fluid = HallThruster.Xenon

    function source!(Q, U, params, i)
        #HallThruster.apply_reactions!(Q, U, params, i)
        HallThruster.apply_ion_acceleration!(Q, U, params, i)
        #HallThruster.source_electron_energy!(Q, U, params, i)
        return Q
    end

    function source_potential!(b, U, s_consts)
        HallThruster.potential_source_term!(b, U, s_consts)
        #HallThruster.OVS_potential_source_term!(b, s_consts)
    end
    
    function boundary_potential!(A, b, U, bc_consts)
        ϕ_L = 400.0
        ϕ_R = 0.0
        HallThruster.boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
        #HallThruster.OVS_boundary_conditions_potential!((A, b, U, bc_consts, ϕ_L, ϕ_R)
    end

    function IC!(U, z, fluids, L)
        ρ1 = 1e19 * fluid.m
        u1 = 300.0
        U .= [ρ1, ρ1, ρ1*u1, 0.0] #[ρ1, ρ1*u1, ρ1*E]
        return U
    end

    ρ1 = 1.0
    u1 = 300.0
    T1 = 300.0

    left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    left_state_elec = 0.0
    right_state_elec = left_state_elec
    BCs_elec = (HallThruster.Dirichlet_energy(left_state_elec), HallThruster.Dirichlet_energy(right_state_elec))

    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
        scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
        initial_condition = IC!,
        source_term! = source!,
        source_potential! = source_potential!,
        boundary_potential! = boundary_potential!, 
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, T1));
        HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(T1))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time,
        saveat = [end_time],
        timestepcontrol = (dt, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = nothing,
        solve_energy = false
    )

    sol = HallThruster.run_simulation(sim)

    #check if ion exit velocity corresponds to energy conservation
    println("ion exit velocity: $(sol.u[1][3, end]./sol.u[1][2, end])")
    @test sol.u[1][3, end]./sol.u[1][2, end] ≈ sqrt(2*400*HallThruster.e/fluid.m) atol = 1000

end

function test_ionization_source(fluxfn, reconstruct, end_time, dt)

    fluid = HallThruster.Xenon

    function source!(Q, U, params, i)
        HallThruster.apply_reactions!(Q, U, params, i)
        #HallThruster.apply_ion_acceleration!(Q, U, params, i)
        #HallThruster.source_electron_energy!(Q, U, params, i)
        return Q
    end

    function source_potential!(b, U, s_consts)
        HallThruster.potential_source_term!(b, U, s_consts)
        #HallThruster.OVS_potential_source_term!(b, s_consts)
    end
    
    function boundary_potential!(A, b, U, bc_consts)
        ϕ_L = 400.0
        ϕ_R = 0.0
        HallThruster.boundary_conditions_potential!(A, b, U, bc_consts, ϕ_L, ϕ_R)
        #HallThruster.OVS_boundary_conditions_potential!((A, b, U, bc_consts, ϕ_L, ϕ_R)
    end

    function IC!(U, z, fluids, L)
        ρ1 = 1e19 * fluid.m
        ρ2 = 0.01 * ρ1
        u1 = 300.0
        Tev = 10.0
        ne = ρ2 / fluid.m
        U .= SA[ρ1, ρ2, ρ2*u1, 3/2*ne*Tev]
        return U
    end

    ρ1 = 1e19 * fluid.m
    ρ2 = 0.01 * ρ1
    u1 = 300.0
    T1 = 300.0

    left_state = [ρ1, ρ2, ρ2*u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    left_state_elec = 0.0
    right_state_elec = left_state_elec
    BCs_elec = (HallThruster.Dirichlet_energy(left_state_elec), HallThruster.Dirichlet_energy(right_state_elec))

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
                #U[index.Tev, i] = max(1, U[index.nϵ, i]/3*2/U[index.ne, i])
            end
            U[index.ne, i] = max(1e-10, HallThruster.electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
            U[index.pe, i] = HallThruster.electron_pressure(U[index.ne, i], U[index.Tev, i]) #this would be real electron pressure, ie next step use for previous in energy convection update
            #U[index.pe, i] = U[index.nϵ, i]/3*2*HallThruster.e #if using the same for pe and ne, might solve some instabilities
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

    callback = cb


    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
        scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
        initial_condition = IC!,
        source_term! = source!,
        source_potential! = source_potential!,
        boundary_potential! = boundary_potential!, 
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, 300.0));
        HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time, #0.0002
        saveat = [end_time],
        timestepcontrol = (dt, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = callback,
        solve_energy = false
    )

    sol = HallThruster.run_simulation(sim)

    @show ρ1, ρ2
    @show sol.u[1][1, end]
    @show sol.u[1][2, end]

    #see if ion and neutral mass fractions changed
    @test sol.u[1][1, end]/ρ1 ≈ 0.01 atol = 0.02
    @test sol.u[1][2, end]/ρ1 ≈ 1 rtol = 0.2
    @test sol.u[1][3, end]/sol.u[1][2, end] ≈ 300 atol = 30

end