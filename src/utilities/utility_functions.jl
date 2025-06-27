"""
    @public
Replacement for the new `public` keyword in v1.11, which is not present in the LTS version v1.10.
When a new LTS version is declared, this can be safely removed and replaced with the bare `public` keyword.
"""
macro public(ex)
    return if VERSION >= v"1.11.0-DEV.469"
        args = ex isa Symbol ? (ex,) : Base.isexpr(ex, :tuple) ? ex.args : error("something informative")
        esc(Expr(:public, args...))
    else
        nothing
    end
end

"""
Placeholder to allow Unitful and DynamicalQuantities
"""
units(::Any) = 1

"""
Convert a number to Float64, first converting it to the appropriate unit, if required
"""
convert_to_float64(number::Number, @nospecialize unit::Any) = Float64(number)

#============================================
    Basic statistics functions
============================================#
mean(x) = sum(x) / length(x)

function var(x)
    μ = mean(x)
    return mean((_x - μ)^2 for _x in x)
end

std(x) = sqrt(var(x))

#============================================
    Solver utility functions
============================================#

@inline left_edge(i) = i - 1
@inline right_edge(i) = i

@inline function inlet_neutral_density(propellant, channel_area)
    return propellant.flow_rate_kg_s / (propellant.velocity_m_s * channel_area)
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
