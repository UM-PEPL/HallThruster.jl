using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, DifferentialEquations, Statistics, Plots

Te_func = z -> 30 * exp(-(2(z - SPT_100.channel_length) / 0.033)^2)
ϕ_func = z -> 300 * (1 - 1/(1 + exp(-1000 * (z - SPT_100.channel_length))))
ni_func = z -> 2000.0

end_time = 1e-5

const SPT_100 = (
    domain = (0.0, 0.05),
    channel_length = 0.025,
    inner_radius = 0.0345,
    outer_radius = 0.05
)

simulation = (
    ncells = 20,
    propellant = HallThruster.Xenon,
    ncharge = 1,
    geometry = SPT_100,
    neutral_temperature = 500.,
    neutral_velocity = 300.,
    ion_temperature = 0.0,
    initial_Te = Te_func,
    initial_ϕ = ϕ_func,
    initial_ni = ni_func,
    solve_Te = false,
    solve_ne = false,
    inlet_mdot = 5e-6,
    saveat = [end_time],
    tspan = (0., end_time),
    dt = 5e-8,
    scheme = (
        flux_function = HallThruster.upwind!,
        limiter = identity,
        reconstruct = false
    ),
)

sol = HallThruster.run_simulation(simulation)