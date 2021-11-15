function test_ion_accel_source(fluxfn, reconstruct, end_time, dt)

    fluid = HallThruster.Xenon

    function source!(Q, U, params, i)
        #HallThruster.apply_reactions!(Q, U, params, i)
        HallThruster.apply_ion_acceleration!(Q, U, params, i)
        return Q
    end
    
    function IC!(U, z, fluids, L)
        ρ1 = 1.0
        u1 = 300.0
        U .= [ρ1, ρ1, ρ1*u1] #[ρ1, ρ1*u1, ρ1*E]
        return U
    end
    
    ρ1 = 1.0
    u1 = 300.0
    T1 = 300.0
    
    left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E] 
    right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())
    
    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, 100), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, reconstruct),
    initial_condition = IC!, 
    source_term! = source!, 
    fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, 300.0));
    HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
    #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
    end_time = end_time, 
    saveat = [end_time],
    timestepcontrol = (dt, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = nothing
    )
    
    sol = HallThruster.run_simulation(sim)

    #check if ion exit velocity corresponds to energy conservation 
    println("ion exit velocity: $(sol.u[1][3, end]./sol.u[1][2, end])")
    @test sol.u[1][3, end]./sol.u[1][2, end]  ≈ sqrt(2*400*HallThruster.e/fluid.m) atol = 10000

end