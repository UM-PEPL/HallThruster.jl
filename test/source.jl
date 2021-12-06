function test_ion_accel_source(fluxfn, reconstruct, end_time, dt)

    fluid = HallThruster.Xenon

    function source!(Q, U, params, ϕ, Tev, i)
        #HallThruster.apply_reactions!(Q, U, params, i)
        HallThruster.apply_ion_acceleration!(Q, U, params, ϕ, i)
        return Q
    end

    function source_potential!(b, i, i_f, μ⁻, μ⁺, pe, U, Δz)
        HallThruster.potential_source_term!(b, i, i_f, μ⁻, μ⁺, pe, U, Δz)
        #HallThruster.OVS_potential_source_term!(b, i)
    end
    
    function boundary_potential!(U, fluid, N, pe, ne, B, A, b, Tev, νan, Δz, OVS)
        ϕ_L = 400.0
        ϕ_R = 0.0
        HallThruster.boundary_conditions_potential!(U, fluid, N, pe, ne, B, A, b, Tev, νan, ϕ_L, ϕ_R, Δz)
        #HallThruster.OVS_boundary_conditions_potential!(N, A, b, ϕ_L, ϕ_R, Δz, OVS)
    end

    function IC!(U, z, fluids, L)
        ρ1 = 1e19 * fluid.m
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
        source_potential! = source_potential!,
        boundary_potential! = boundary_potential!, 
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, T1));
        HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(T1))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time,
        saveat = [end_time],
        timestepcontrol = (dt, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = nothing
    )

    sol = HallThruster.run_simulation(sim)

    #check if ion exit velocity corresponds to energy conservation
    println("ion exit velocity: $(sol.u[1][3, end]./sol.u[1][2, end])")
    @test sol.u[1][3, end]./sol.u[1][2, end] ≈ sqrt(2*400*HallThruster.e/fluid.m) atol = 1000

end

function test_ionization_source(fluxfn, reconstruct, end_time, dt)

    fluid = HallThruster.Xenon

    function source!(Q, U, params, ϕ, Tev, i)
        HallThruster.apply_reactions!(Q, U, params, Tev, i)
        #HallThruster.apply_ion_acceleration!(Q, U, params, i)
        return Q
    end

    function source_potential!(b, i, i_f, μ⁻, μ⁺, pe, U, Δz)
        HallThruster.potential_source_term!(b, i, i_f, μ⁻, μ⁺, pe, U, Δz)
        #HallThruster.OVS_potential_source_term!(b, i)
    end
    
    function boundary_potential!(U, fluid, N, pe, ne, B, A, b, Tev, νan, Δz, OVS)
        ϕ_L = 400.0
        ϕ_R = 0.0
        HallThruster.boundary_conditions_potential!(U, fluid, N, pe, ne, B, A, b, Tev, νan, ϕ_L, ϕ_R, Δz)
        #HallThruster.OVS_boundary_conditions_potential!(N, A, b, ϕ_L, ϕ_R, Δz, OVS)
    end

    function IC!(U, z, fluids, L)
        ρ1 = 1e19 * fluid.m
        ρ2 = 0.01 * ρ1
        u1 = 300.0
        U .= [ρ1, ρ2, ρ2*u1] #[ρ1, ρ1*u1, ρ1*E]
        return U
    end

    ρ1 = 1e19 * fluid.m
    ρ2 = 0.01 * ρ1
    u1 = 300.0
    T1 = 300.0

    left_state = [ρ1, ρ2, ρ2*u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = BCs,
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
        callback = nothing
    )

    sol = HallThruster.run_simulation(sim)

    @show ρ1, ρ2
    @show sol.u[1][1, end]
    @show sol.u[1][2, end]

    #see if ion and neutral mass fractions changed
    @test sol.u[1][1, end]/ρ1 ≈ 0.98 rtol = 0.2
    @test sol.u[1][2, end]/ρ2 ≈ 1.5 rtol = 0.2
    @test sol.u[1][3, end]/sol.u[1][2, end] ≈ 300 atol = 30

end