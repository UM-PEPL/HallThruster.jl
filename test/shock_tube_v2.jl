using HallThruster, Plots, SodShockTube, CairoMakie

function shock_tube(fluxfn, ncells, end_time)
    ρL, uL, pL = 1.0, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1

    gas = HallThruster.Air

    TL = pL / gas.R / ρL
    TR = pR / gas.R / ρR

    EL = gas.cv * TL
    ER = gas.cv * TR

    left_state = [ρL, ρL * uL, ρL * EL, ρL , ρL * uL, ρL * EL]
    right_state = [ρR, ρR * uR, ρR * ER, ρR, ρR * uR, ρR * ER]

    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    L = 1.0

    geometry = (
        domain = (0.0, L),
        channel_length = L/2,
        inner_radius = 0.0345,
        outer_radius = 0.05
    )

    function source!(Q, U, params, i)
       return Q
    end
        
    function IC!(U, z, fluids, L)
        ρL, uL, pL = 1.0, 0.0, 1.0
        ρR, uR, pR = 0.125, 0.0, 0.1
        gas1 = fluids[1].species.element
        gas2 = fluids[2].species.element
        if z < L/2
            ρ1 = ρL
            T1 = pL / ρL * gas1.m / HallThruster.kB
            ρ2 = ρL
            T2 = pL / ρL * gas2.m / HallThruster.kB
        else
            ρ1 = ρR
            T1 = pR / ρR * gas1.m / HallThruster.kB
            ρ2 = ρR
            T2 = pR / ρR * gas2.m / HallThruster.kB

        end
        u1 = 0.0
        u2 = 0.0
        U .= [ρ1, ρ1*u1, ρ1*T1*gas1.cv, ρ2, ρ2*u2, ρ2*T2*gas2.cv, 0.0]
        return U
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

    #simulation input, need to somehow do initial conditions now as well, do as with boundary conditions in thomas code
    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(geometry, ncells), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, true),
    initial_condition = IC!, 
    source_term! = source!, 
    source_potential! = source_potential!,
    boundary_potential! = boundary_potential!,
    fluids = [HallThruster.Fluid(HallThruster.Species(HallThruster.Air, 0), HallThruster.EulerEquations());
    HallThruster.Fluid(HallThruster.Species(HallThruster.Air, 0), HallThruster.EulerEquations())], 
    end_time = end_time, 
    saveat = [0, end_time],
    timestepcontrol = (1e-6, true), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
    callback = nothing,
    solve_energy = false
    )
    sol = HallThruster.run_simulation(sim)
end

function run_shock_tube(ncells, time)
    #analytical
    problem = ShockTubeProblem(
        geometry = (0.0, 1.0, 0.5), # left edge, right edge, initial shock location
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = time,
        γ = 1.4
    );
    xs = LinRange(0.0, 1.0, ncells+2); # x locations at which to solve
    positions, regions, values = SodShockTube.solve(problem, xs);

    #discretised
    @time sol = shock_tube(HallThruster.HLLE!, ncells, time)
    ρ_d = sol.u[2][1, :]
    ρu_d = sol.u[2][2, :]
    ρE_d = sol.u[2][3, :]
    zs = LinRange(0, 1, ncells+2)
    @show sol.t[end]
    
    #plot the discrepancies
    f = Figure(resolution = (1000, 1000))
    ax_ρ = Axis(f[1,1], xlabel = "x", ylabel = "ρ", title = "Density")
    ax_u = Axis(f[2,1], xlabel = "x", ylabel = "u", title = "Velocity")
    ax_p = Axis(f[1,2], xlabel = "x", ylabel = "p", title = "Pressure")
    ax_E = Axis(f[2,2], xlabel = "x", ylabel = "E", title = "Stagnation/Internal Energy")

    opts = (;linewidth = 4)

    lines!(ax_ρ, values.x, values.ρ; opts...)
    lines!(ax_ρ, values.x, ρ_d; opts...)
    lines!(ax_u, values.x, values.u; opts...)
    lines!(ax_u, values.x, ρu_d./ρ_d; opts...)
    lines!(ax_p, values.x, values.p; opts...)
    lines!(ax_p, values.x, @.((1.4-1)*(ρE_d-0.5*ρu_d*ρu_d/ρ_d)); opts...)
    lines!(ax_E, values.x, values.e; opts...) #stagnation energy
    lines!(ax_E, values.x, @.(ρE_d + 0.5*ρu_d*ρu_d/ρ_d/ρ_d*(1 -ρ_d)); opts...) #stagnation energy
    lines!(ax_E, values.x, @.((values.e - 0.5*values.u*values.u)/values.ρ); opts...) #internal energy
    lines!(ax_E, values.x, @.((ρE_d - 0.5*ρu_d*ρu_d/ρ_d)/ρ_d); opts...) #internal energy
    display(f)
    #save("shocktube.png", f) 
end