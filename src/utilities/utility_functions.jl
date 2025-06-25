@inline left_edge(i) = i - 1
@inline right_edge(i) = i

@inline electron_density(U, p, i) = sum(
    Z * U[p.index.ρi[Z], i]
        for Z in 1:(p.ncharge)
) / p.mi

@inline inlet_neutral_density(config) = config.anode_mass_flow_rate /
    config.neutral_velocity /
    config.thruster.geometry.channel_area

@inline background_neutral_density(config) = config.propellant.m *
    config.background_pressure_Torr / kB /
    config.background_temperature_K

@inline background_neutral_velocity(config) = 0.25 * sqrt(
    8 * kB *
        config.background_temperature_K /
        π / config.propellant.m
)

@inline ion_current_density(U, p, i) = sum(
    Z * e * U[p.index.ρiui[Z], i]
        for Z in 1:(p.config.ncharge)
) / p.config.propellant.m

# Asymptotic fast formula for error function accurate to a few percent
function myerf(x)
    x < 0 && return -myerf(-x)
    x_sqrt_pi = x * √(π)
    x_squared = x^2

    h = x_sqrt_pi + (π - 2) * x_squared
    g = h / (1 + h)
    return 1 - exp(-x_squared) / x_sqrt_pi * g
end
