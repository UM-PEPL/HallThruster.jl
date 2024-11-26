"""
    AnomalousTransportModel
The abstract supertype of all types of anomalous transport models.
Subtype this to define your own model.
"""
abstract type AnomalousTransportModel end

#=============================================================================
 Serialization
==============================================================================#
"""
$(SIGNATURES)
Returns a NamedTuple mapping symbols to transport models for all built-in models.
"""
@inline function anom_models()
    return (; NoAnom, Bohm, TwoZoneBohm,
        MultiLogBohm, GaussianBohm,
        LogisticPressureShift,)
end

function Serialization.SType(::Type{T}) where {T <: AnomalousTransportModel}
    Serialization.TaggedUnion()
end
Serialization.options(::Type{T}) where {T <: AnomalousTransportModel} = anom_models()

#=============================================================================
 Begin definition of built-in models
==============================================================================#

"""
    NoAnom <: AnomalousTransportModel
No anomalous collision frequency included in simulation
"""
struct NoAnom <: AnomalousTransportModel end

function (::NoAnom)(νan, ::Any)
    for i in eachindex(νan)
        νan[i] = 0.0
    end
    return νan
end

"""
    Bohm(c) <: AnomalousTransportModel
Model where the anomalous collision frequency scales with the electron cyclotron frequency ωce times some scaling factor c
"""
@kwdef struct Bohm <: AnomalousTransportModel
    c::Float64
end

function (model::Bohm)(νan, params)
    (; config, cache, grid) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    z_shift = pressure_shift(config.anom_model, params)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (_, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(z)
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
@kwdef struct TwoZoneBohm <: AnomalousTransportModel
    c1::Float64
    c2::Float64
end

function (model::TwoZoneBohm)(νan, params)
    (; c1, c2) = model
    (; config, cache, grid) = params
    (; B) = cache

    L_trans = config.transition_length
    L_ch = config.thruster.geometry.channel_length

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    z_shift = pressure_shift(config.anom_model, params)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(z)
        ωce = e * B / me
        c = linear_transition(z, L_ch, L_trans, c1, c2)
        νan[i] = c * ωce
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
@kwdef struct MultiLogBohm <: AnomalousTransportModel
    zs::Vector{Float64}
    cs::Vector{Float64}
    function MultiLogBohm(zs, cs)
        if length(zs) != length(cs)
            throw(ArgumentError("Number of z values must be equal to number of c values"))
        end
        return new(zs, cs)
    end
end

function (model::MultiLogBohm)(νan, params)
    (; grid, config) = params
    (; B) = params.cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    z_shift = pressure_shift(config.anom_model, params)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(z)
        ωce = e * B / me
        c = HallThruster.interpolate(z, model.zs, model.cs, use_log = true)
        νan[i] = c * ωce
    end

    return νan
end

"""
    GaussianBohm(hall_min, hall_max, center, width) <: AnomalousTransportModel
Model in which the anomalous collision frequency is Bohm-like (ν_an ~ ω_ce), 
except in a Gaussian-shaped region defined centered on z = `center`,
where the collision frequency is lower.

# Arguments
- `hall_min`: the minimum Hall parameter
- `hall_max`: the maximum Hall parameter
- `center`: the axial position (in meters) of the mean of the Gaussian trough
- `width`: the standard deviation (in meters) of the Gaussian trough
"""
@kwdef struct GaussianBohm <: AnomalousTransportModel
    hall_min::Float64
    hall_max::Float64
    center::Float64
    width::Float64
end

function (model::GaussianBohm)(νan, params)
    (; hall_min, hall_max, center, width) = model
    (; config, cache, grid) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    z_shift = pressure_shift(config.anom_model, params)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(z)
        ωce = e * B / me
        c = hall_max * (1 - (1 - hall_min) * exp(-0.5 * ((z - center) / width)^2))
        νan[i] = c * ωce
    end

    return νan
end

abstract type PressureShift <: AnomalousTransportModel end

"""
    LogisticPressureShift(model, z0, dz, pstar, alpha)
A wrapper model that allows a transport profile to shift axially in response
to changes in background pressure. The displacement/shift of the transport profile
follows a logistic curve.

# Arguments
- model: An anomalous transport model
- z0: the center of the shift at 0 background pressure
- dz: the total pressure shift across (0, Inf) background pressure
- pstar: the "turning point" pressure
- alpha: the slope of the pressure-displacement response curve
"""
@kwdef struct LogisticPressureShift{A <: AnomalousTransportModel} <: PressureShift
    model::A
    z0::Float64
    dz::Float64
    pstar::Float64
    alpha::Float64
end

function (model::LogisticPressureShift)(args...; kwargs...)
    return model.model(args...; kwargs...)
end

pressure_shift(model::AnomalousTransportModel, ::Any) = 0.0

function pressure_shift(model::LogisticPressureShift, params)
    (; z0, dz, alpha, pstar) = model
    pb = params.config.background_pressure_Torr

    torr_to_pa = 133.322
    p_ratio = pb / (pstar * torr_to_pa)
    zstar = z0 + dz / (1 + (alpha - 1)^(2 * p_ratio - 1))
    return zstar
end

"""
    num_anom_variables(::AnomalousTransportModel)::Int

The number of variable arrays that should be allocated for the provided anomalous
transport model. These arrays are used to save state beyond the anomalous
collision frequency, and are useful for defining more complex anomalous transport
models. If not defined by the user, this defaults to zero.
"""
num_anom_variables(::AnomalousTransportModel)::Int = 0
