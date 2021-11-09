using HallThruster, Plots

function shock_tube_v2(fluxfn, ncells, end_time)
    ρL, uL, pL = 1.0, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1

    gas = HallThruster.Air

    TL = pL / gas.R / ρL
    TR = pR / gas.R / ρR

    EL = gas.cv * TL
    ER = gas.cv * TR

    left_state = [ρL, ρL * uL, ρL * EL]
    right_state = [ρR, ρR * uR, ρR * ER]

    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    L = 1.0

    geometry = (
        domain = (0.0, L),
        channel_length = L/2,
        inner_radius = 0.0345,
        outer_radius = 0.05
    )

    ni_func = x -> if x < L/2
        ρL / gas.m
    else
        ρR / gas.m
    end

    ui_func = x -> 0.0

    Ti_func = x -> if x < L/2
        pL / ρL * gas.m / HallThruster.kB
    else
        pR / ρR * gas.m / HallThruster.kB
    end

    function source(z)
       return 0.0
    end
    
    IC = [ni_func, ui_func, Ti_func]

    #simulation input, need to somehow do initial conditions now as well, do as with boundary conditions in thomas code
    sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(geometry, ncells), 
    boundary_conditions = BCs,
    scheme = HallThruster.HyperbolicScheme(fluxfn, HallThruster.minmod, true),
    initial_condition = IC, 
    source_term! = source, 
    fluids = [HallThruster.Fluid(HallThruster.Species(HallThruster.Air, 0), HallThruster.EulerEquations())],  
    end_time = end_time, 
    saveat = [0, end_time]
    )
    sol = HallThruster.run_simulation_v2(sim)
end

function run_shock_tube_v2(ncells, time)
    @time sol = shock_tube_v2(HallThruster.HLLE!, ncells, time)
    zs = LinRange(0, 1, ncells+2)
    p = plot(zs, sol.u[1][1, :])
    plot!(zs, sol.u[end][1, :])
    display(p)
    @show sol.t[end]
    sol

end
