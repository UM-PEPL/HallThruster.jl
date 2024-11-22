"""
    AnomalousTransportModel
The abstract supertype of all types of anomalous transport models.
Subtype this to define your own model.
"""
abstract type AnomalousTransportModel end
#=============================================================================
 Begin definition of built-in models
==============================================================================#

"""
    NoAnom <: AnomalousTransportModel
No anomalous collision frequency included in simulation
"""
struct NoAnom <: AnomalousTransportModel end

function (::NoAnom)(νan, params)
    for i in eachindex(νan)
        νan[i] = 0.0
    end
    return νan
end

"""
    Bohm(c) <: AnomalousTransportModel
Model where the anomalous collision frequency scales with the electron cyclotron frequency ωce times some scaling factor c
"""
struct Bohm <: AnomalousTransportModel
    c::Float64
end

function (model::Bohm)(νan, params)
    for i in eachindex(νan)
        B = params.cache.B[i]
        ωce = e * B / me
        νan = model.c * ωce
    end

    return νan
end

"""
    TwoZoneBohm(c1, c2) <: AnomalousTransportModel
Model where the anomalous collision frequency has two values: c1 * ωce inside the channel and c2 * ωce outside of the channel.
Takes two arguments: c1 and c2. The transition between these values can be smoothed by the user-provided transition function.
"""
struct TwoZoneBohm <: AnomalousTransportModel
    coeffs::NTuple{2, Float64}
    TwoZoneBohm(c1, c2) = new((c1, c2))
    TwoZoneBohm(t) = new(t)
end

function (model::TwoZoneBohm)(νan, params)
    (; grid, config) = params
    L_ch = config.thruster.geometry.channel_length
    c1, c2 = model.coeffs

    for i in eachindex(νan)
        B = params.cache.B[i]
        ωce = e * B / me

        β = linear_transition(
            grid.cell_centers[i], L_ch, config.transition_length, c1, c2,)
        νan[i] = β * ωce
    end

    return νan
end

"""
    MultiLogBohm(zs, cs) <: AnomalousTransportModel
Model similar to that employed in Hall2De, where the mobility is Bohm-like (i.e. `νan(z) = c(z) * ωce(z)`) and z is in meters.

The function `c(z)` is defined by a sequence of nodes `(z, c)` provided by the user. At `z = z[1]`, `c(z) = c[1]`, and so forth.

At `z[i] < z < z[i+1]`, `log(c)` is defined by linearly interpolating between `log(c[i])` and `log(c[i+1])`.

For `z < z[1]`, `c = c[1]` and for `z > z[end]`, `c(z) = c[end]`.

The user may also provide a single array of [z[1], z[2], ..., z[end], c[1], c[2], ..., c[end]]. The number of z values must be equal to the number of c values.
"""
struct MultiLogBohm <: AnomalousTransportModel
    zs::Vector{Float64}
    cs::Vector{Float64}
    function MultiLogBohm(zs, cs)
        if length(zs) != length(cs)
            throw(ArgumentError("Number of z values must be equal to number of c values"))
        end
        return new(zs, cs)
    end
end

function MultiLogBohm(coeffs)
    N = length(coeffs) ÷ 2
    zs = coeffs[1:N]
    cs = coeffs[(N + 1):end]
    return MultiLogBohm(zs, cs)
end

function (model::MultiLogBohm)(νan, params)
    (; grid) = params
    (; B) = params.cache

    for i in eachindex(νan)
        ωce = e * B[i] / me
        c = HallThruster.interpolate(
            grid.cell_centers[i], model.zs, model.cs, use_log = true,)
        νan[i] = c * ωce
    end

    return νan
end

"""
    ShiftedMultiBohm(zs, cs, z0, dz, alpha, pstar) <: AnomalousTransportModel
A version of MultiLogBohm where the coefficients are shifted depending on the background pressure.
"""
Base.@kwdef struct ShiftedMultiBohm <: HallThruster.AnomalousTransportModel
    zs::Vector{Float64}
    cs::Vector{Float64}
    z0::Float64
    dz::Float64
    pstar::Float64
    alpha::Float64
end

function (model::ShiftedMultiBohm)(νan, params)
    (; zs, cs, z0, dz, alpha, pstar) = model
    (; config, grid) = params

    pb = config.background_pressure

    torr_to_pa = 133.322

    p_ratio = pb / (pstar * torr_to_pa)
    zstar = z0 + dz / (1 + (alpha - 1)^(2 * p_ratio - 1))

    for i in eachindex(νan)
        B = params.cache.B[i]
        ωce = e * B / me
        c = HallThruster.interpolate(grid.cell_centers[i] - zstar, zs, cs, use_log = true)
        νan[i] = c * ωce
    end
end

"""
    ShiftedGaussianBohm(trough_location, trough_width, trough_depth, z0, dz, alpha, pstar) <: AnomalousTransportModel
Model in which the anomalous collision frequency is Bohm-like (ν_an ~ ω_ce), except in a Gaussian-shaped region defined by
the parameters trough_location, trough_width, trough_max, and trough_min where the collision frequency is lower.
The location of the trough is based on the background pressure and the user-provided coefficients.

# Arguments
- `trough_min`: the minimum Hall parameter
- `trough_max`: the maximum Hall parameter
- `trough_location`: the axial position (in meters) of the mean of the Gaussian trough
- `trough_width`: the standard deviation (in meters) of the Gaussian trough
- `z0`: the furthest upstream displacement permitted at high back-pressures, relative to `trough_location`
- `dz`: the maximum allowable amount of axial displacement
- `pstar`: the background pressure at which the shift upstream halts/plateaus
- `alpha`: the slope of the pressure response curve, with a higher value corresponding to a steeper pressure response
"""
Base.@kwdef struct ShiftedGaussianBohm <: HallThruster.AnomalousTransportModel
    trough_min::Float64
    trough_max::Float64
    trough_location::Float64
    trough_width::Float64
    z0::Float64
    dz::Float64
    pstar::Float64
    alpha::Float64
end

function (model::ShiftedGaussianBohm)(νan, params)
    (; grid) = params
    (; trough_location, trough_width, trough_min, trough_max, z0, dz, alpha, pstar) = model

    # do not recompute after a certain point - profile is meant to be fixed in time
    if (params.iteration[] > 5)
        return νan
    end

    pb = params.config.background_pressure
    torr_to_pa = 133.322

    p_ratio = pb / (pstar * torr_to_pa)
    zstar = z0 + dz / (1 + (alpha - 1)^(2 * p_ratio - 1))

    B_interp = LinearInterpolation(params.z_cell, params.cache.B)

    for i in eachindex(νan)
        z = grid.cell_centers[i] - zstar
        B = B_interp(z)
        ωce = e * B / me
        μ = trough_location
        c = trough_max * (1 - (1 - trough_min) * exp(-0.5 * ((z - μ) / trough_width)^2))
        νan[i] = c * ωce
    end

    return νan
end

"""
    ShiftedTwoZoneBohm(coeffs, z0, dz, alpha, pstar) <: AnomalousTransportModel
Model where the anomalous collision frequency has two values: c1 * ωce before some transition location and c2 * ωce after.
Takes two arguments: c1 and c2. The transition between these values can be smoothed by the user-provided transition function.
The location of the transition is based on the background pressure and the user-provided coefficients.
"""
Base.@kwdef struct ShiftedTwoZoneBohm <: HallThruster.AnomalousTransportModel
    coeffs::NTuple{2, Float64}
    z0::Float64
    dz::Float64
    pstar::Float64
    alpha::Float64
end

function (model::ShiftedTwoZoneBohm)(νan, params)
    (; coeffs, z0, dz, alpha, pstar) = model
    (; grid, config) = params

    L_ch = config.thruster.geometry.channel_length
    pb = config.background_pressure

    torr_to_pa = 133.322

    p_ratio = pb / (pstar * torr_to_pa)
    zstar = L_ch + z0 + dz / (1 + (alpha - 1)^(2 * p_ratio - 1))

    c1, c2 = coeffs

    for i in eachindex(νan)
        B = params.cache.B[i]
        ωce = HallThruster.e * B / HallThruster.me

        β = linear_transition(
            grid.cell_centers, zstar, config.transition_length, c1, c2,)
        νan[i] = β * ωce
    end
end

"""
    num_anom_variables(::AnomalousTransportModel)::Int

The number of variable arrays that should be allocated for the provided anomalous
transport model. These arrays are used to save state beyond the anomalous
collision frequency, and are useful for defining more complex anomalous transport
models. If not defined by the user, this defaults to zero.
"""
num_anom_variables(::AnomalousTransportModel)::Int = 0

"""
    allocate_anom_variables(::AnomalusTransportModel, ncells)
Allocate arrays for anomalous transport state variables. `ncells` is the length
of the arrays to be allocated. These anomalous transport variables are then stored
in params.cache.anom_variables
"""
function allocate_anom_variables(model::AnomalousTransportModel, ncells)
    [zeros(ncells) for _ in 1:num_anom_variables(model)]
end
