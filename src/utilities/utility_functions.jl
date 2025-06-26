@inline left_edge(i) = i - 1
@inline right_edge(i) = i

# TODO: multiple gases + fluid containers
@inline function electron_density(U, p, i)
    prop = p.propellants[1]
    return sum(Z * U[p.index.ρi[Z], i] for Z in 1:prop.max_charge) / prop.gas.m
end

# TODO: mutliple gases + fluid containers
@inline function inlet_neutral_density(config)
    prop = config.propellants[1]
    return prop.flow_rate_kg_s / (prop.velocity_m_s * config.thruster.geometry.channel_area)
end

# TODO: multiple gases + fluid containers
@inline function background_neutral_density(config)
    return config.propellants[1].gas.m * config.background_pressure_Torr / kB / config.background_temperature_K
end

# TODO: multiple gases + fluid containers
@inline function background_neutral_velocity(config)
    return 0.25 * sqrt(8 * kB * config.background_temperature_K / π / config.propellants[1].gas.m)
end

# Asymptotic fast formula for error function accurate to a few percent
function myerf(x)
    x < 0 && return -myerf(-x)
    x_sqrt_pi = x * √(π)
    x_squared = x^2

    h = x_sqrt_pi + (π - 2) * x_squared
    g = h / (1 + h)
    return 1 - exp(-x_squared) / x_sqrt_pi * g
end
