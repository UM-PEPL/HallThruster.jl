using HallThruster

function shock_tube(fluxfn, ncells, end_time)
    ρL, uL, pL = 1.0, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1

    gas = HallThruster.Air


    L = 1.0

    geometry = (
        domain = (0.0, L),
        channel_length = L/2,
        inner_radius = 0.0345,
        outer_radius = 0.05
    )

    un = uL
    Tn = 0.003484320557491289

    Te_func = x -> 0.0
    ϕ_func = x -> 0.0
    ni_func = x -> if x < L/2
        ρL / gas.m
    else
        ρR / gas.m
    end


    Ti_func = x -> if x < L/2
        pL / ρL * gas.m / HallThruster.kB
    else
        pR / ρR * gas.m / HallThruster.kB
    end

    simulation = (
        ncells = ncells,
        propellant = HallThruster.Air,
        ncharge = 1,
        MMS = false,
        mms! = nothing,
        cb = nothing,#DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass),  #abstol, reltol, test
        geometry = geometry,
        neutral_temperature = Tn,
        neutral_velocity = un,
        ion_temperature = 0.0,
        initial_nn_mms = nothing,
        initial_Te = Te_func,
        initial_ϕ = ϕ_func,
        initial_ni = ni_func,
        initial_Ti = Ti_func,
        solve_Te = false,
        solve_ne = false,
        inlet_mdot = 5e-6,
        saveat = [0, end_time],
        tspan = (0., end_time),
        dt = 1e-8, #5e-8
        scheme = (
            flux_function = fluxfn,
            limiter = identity,
            reconstruct = false
        ),
    )

    sol = HallThruster.run_simulation(simulation)
end

function run_shock_tube(ncells, time)
    @time sol = shock_tube(HallThruster.HLLE!, ncells, time)
    zs = LinRange(0, 1, ncells+2)
    p = plot(zs, sol.u[1][2, :])
    plot!(zs, sol.u[end][2, :])
    display(p)
    @show sol.t[end]
    sol

end
