Returns(x) = _ -> x

function freestream_preservation(;CFL, fluxfn, reconstruct)
    Te_func = Returns(0.0)
    ϕ_func = Returns(0.0)
    ni_func = Returns(1e6)

    end_time = 1e-3

    geometry = HallThruster.SPT_100
    ncells = 100
    un = 300.0
    dx = (geometry.domain[2] - geometry.domain[1]) / ncells
    dt = CFL * dx / un

    # Test no electric field, no ion temperature, no electron temperature, constant ion density
    simulation = (
        ncells = ncells,
        propellant = HallThruster.Xenon,
        ncharge = 1,
        geometry = geometry,
        neutral_temperature = 500.,
        neutral_velocity = un,
        ion_temperature = 500.0,
        initial_Te = Returns(0.0),
        initial_ϕ = Returns(0.0),
        initial_ni = Returns(1e6),
        solve_Te = false,
        solve_ne = false,
        inlet_mdot = 5e-6,
        tspan = (0., end_time),
        dt = dt,
        scheme = (
            flux_function = fluxfn,
            limiter = identity,
            reconstruct = reconstruct
        ),
        saveat = [0, end_time],
        save_everystep = false,
    )
    @show CFL, dt

    sol = HallThruster.run_simulation(simulation)
    return sol
end

function test_preservation(CFL)
    for (fluxfn, name) in zip([HallThruster.upwind!, HallThruster.HLLE!], ["upwind", "HLLE"])
        println(name)
        sol = freestream_preservation(CFL = CFL, fluxfn = fluxfn, reconstruct = false);
        @test sol.retcode == :Success
        @test sol.u[end] ≈ sol.u[1]
    end
end
