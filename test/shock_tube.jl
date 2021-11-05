using HallThruster

function shock_tube(fluxfn, ncells)
    ρL, uL, pL = 1.0, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1

    gas = HallThruster.Air

    un = 0.0
    Tn = 0.0
    L = 1.0

    Te_func = x -> 0.0
    ϕ_func = x -> 0.0
    ni_func = x -> if x < L
        ρL / gas.m
    else
        ρR / gas.m
    end

    Ti_func = x -> if x < L
        pL / nL / HallThruster.kB
    else
        pR / nR / HallThruster.kB
    end

    end_time = 2.0

    simulation = (
        ncells = ncells,
        propellant = HallThruster.Xenon,
        ncharge = 1.0,
        MMS = false,
        mms! = nothing,
        cb = DiffEqCallbacks.TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass),  #abstol, reltol, test
        geometry = HallThruster.SPT_100,
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
        dt = 200e-8, #5e-8
        scheme = (
            flux_function = fluxfn,
            limiter = identity,
            reconstruct = reconstruct
        ),
    )

    sol = run_simulation(simulation)
end