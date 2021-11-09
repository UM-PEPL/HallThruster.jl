using HallThruster, Plots

function shock_tube_v2(fluxfn, ncells, end_time)
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

    function source(z)
       return 0.0
    end
        
    function IC!(U, z, fluids, L)
        ρL, uL, pL = 1.0, 0.0, 1.0
        ρR, uR, pR = 0.125, 0.0, 0.1
        gas1 = fluids[1].species.element
        gas2 = fluids[2].species.element
        if z < L/2
            n1 = ρL / gas1.m
            T1 = pL / ρL * gas1.m / HallThruster.kB
            n2 = ρL / gas2.m
            T2 = pL / ρL * gas2.m / HallThruster.kB
        else
            n1 = ρR / gas1.m
            T1 = pR / ρR * gas1.m / HallThruster.kB
            n2 = ρR / gas2.m
            T2 = pR / ρR * gas2.m / HallThruster.kB

        end
        u1 = 0.0
        u2 = 0.0
        U .= [n1*gas1.m, n1*u1*gas1.m, n1*gas1.m*T1*gas1.cv, n2*gas2.m, n2*u2*gas2.m, n2*gas2.m*T2*gas2.cv]
        return U
    end


    #simulation input, need to somehow do initial conditions now as well, do as with boundary conditions in thomas code
    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(geometry, ncells), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, identity, false),
    initial_condition = IC!, 
    source_term! = source, 
    fluids = [HallThruster.Fluid(HallThruster.Species(HallThruster.Air, 0), HallThruster.EulerEquations());
    HallThruster.Fluid(HallThruster.Species(HallThruster.Air, 0), HallThruster.EulerEquations())], 
    end_time = end_time, 
    saveat = [0, end_time]
    )
    sol = HallThruster.run_simulation_v2(sim)
end

function run_shock_tube_v2(ncells, time)
    @time sol = shock_tube_v2(HallThruster.HLLE!, ncells, time)
    zs = LinRange(0, 1, ncells+2)
    p = plot(zs, sol.u[1][1, :])
    plot!(zs, sol.u[end][4, :])
    display(p)
    @show sol.t[end]
    sol

end
