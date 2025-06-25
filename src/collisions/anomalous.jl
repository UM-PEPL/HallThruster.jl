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
    return (;
        NoAnom,
        Bohm,
        TwoZoneBohm,
        MultiLogBohm,
        GaussianBohm,
        ScaledGaussianBohm,
        LogisticPressureShift,
        SimpleLogisticShift,
    )
end

function Serialization.SType(::Type{T}) where {T <: AnomalousTransportModel}
    return Serialization.TaggedUnion()
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

function (::NoAnom)(νan, @nospecialize(_x1), @nospecialize(_x2))
    for i in eachindex(νan)
        νan[i] = 0.0
    end
    return νan
end

"""
    Bohm(c) <: AnomalousTransportModel
Model where the anomalous collision frequency scales with the electron cyclotron frequency ωce times some scaling factor c

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct Bohm <: AnomalousTransportModel
    c::Float64
end

function (model::Bohm)(νan, params, config)
    (; cache, grid, thruster) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
    z_shift = pressure_shift(config.anom_model, config.background_pressure_Torr, L_ch)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (_, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(zc)
        ωce = e * B / me
        νan = model.c * ωce
    end

    return νan
end

"""
    TwoZoneBohm(c1, c2) <: AnomalousTransportModel
Model where the anomalous collision frequency has two values: c1 * ωce inside the channel and c2 * ωce outside of the channel.
Takes two arguments: c1 and c2. The transition between these values is smoothed over `config.transition_length`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct TwoZoneBohm <: AnomalousTransportModel
    c1::Float64
    c2::Float64
end

function (model::TwoZoneBohm)(νan, params, config)
    (; c1, c2) = model
    (; cache, grid, thruster) = params
    (; B) = cache

    L_trans = config.transition_length
    L_ch = config.thruster.geometry.channel_length

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
    z_shift = pressure_shift(config.anom_model, config.background_pressure_Torr, L_ch)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(zc)
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

# Fields
$(TYPEDFIELDS)
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

function (model::MultiLogBohm)(νan, params, config)
    (; grid, thruster) = params
    (; B) = params.cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
    z_shift = pressure_shift(config.anom_model, config.background_pressure_Torr, L_ch)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(zc)
        ωce = e * B / me
        c = HallThruster.interpolate(z, model.zs, model.cs, use_log = true)
        νan[i] = c * ωce
    end

    return νan
end

"""
    GaussianBohm(hall_min, hall_max, center, width) <: AnomalousTransportModel
Model in which the anomalous collision frequency is Bohm-like (`νan ~ ω_ce`),
except in a Gaussian-shaped region defined centered on z = `center`,
where the collision frequency is lower.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct GaussianBohm <: AnomalousTransportModel
    """the minimum inverse Hall parameter"""
    hall_min::Float64
    """the maximum inverse Hall parameter"""
    hall_max::Float64
    """the axial position (in meters) of the mean of the Gaussian trough"""
    center::Float64
    """the standard deviation (in meters) of the Gaussian trough"""
    width::Float64
end

function (model::GaussianBohm)(νan, params, config)
    (; hall_min, hall_max, center, width) = model
    (; cache, grid, thruster) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
    z_shift = pressure_shift(config.anom_model, config.background_pressure_Torr, L_ch)
    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(zc)
        ωce = e * B / me
        c = hall_max * (1 - (1 - hall_min) * exp(-0.5 * ((z - center) / width)^2))
        νan[i] = c * ωce
    end

    return νan
end

"""
    ScaledGaussianBohm(anom_scale, barrier_scale, width, center) <: AnomalousTransportModel
Model in which the anomalous collision frequency is Bohm-like (`νan ~ ω_ce`),
except in a Gaussian-shaped region defined centered on z = `center`,
where the collision frequency is lower.
Reparameterized version of the `GaussianBohm` model to make parameters non-dimensional and closer to O(1)

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct ScaledGaussianBohm <: AnomalousTransportModel
    """the maximum inverse hall parameter, should be in [0, 1]"""
    anom_scale::Float64 = 0.0625
    """the factor by which transport is reduced by the baseline value at the center of the trough, should be in [0,1]. """
    barrier_scale::Float64 = 0.9
    """the standard deviation of the Gaussian trough, in channel lengths"""
    width::Float64
    """the axial position of the mean of the Gaussian trough, in channel lengths"""
    center::Float64
end

function (model::ScaledGaussianBohm)(νan, params, config)
    (; anom_scale, barrier_scale, width, center) = model
    (; cache, grid, thruster) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
    z_shift = pressure_shift(config.anom_model, config.background_pressure_Torr, L_ch)
    B_interp = LinearInterpolation(grid.cell_centers, B)
    mean = L_ch * center
    std = L_ch * width

    for (i, zc) in enumerate(grid.cell_centers)
        z = zc - z_shift
        B = B_interp(zc)
        ωce = e * B / me
        c = anom_scale * (1 - barrier_scale * exp(-0.5 * ((z - mean) / (std))^2))
        νan[i] = c * ωce
    end

    return νan
end

abstract type PressureShift <: AnomalousTransportModel end

pressure_shift(model::AnomalousTransportModel, ::Any, ::Any) = 0.0

function (model::PressureShift)(args...; kwargs...)
    return model.model(args...; kwargs...)
end

"""
    LogisticPressureShift(model, z0, dz, pstar, alpha)
A wrapper model that allows a transport profile to shift axially in response to changes in background pressure.
The displacement/shift of the transport profile follows a logistic curve.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct LogisticPressureShift{A <: AnomalousTransportModel} <: PressureShift
    """
    An anomalous transport model
    """
    model::A
    """
    the center of the shift at 0 background pressure
    """
    z0::Float64
    """
    the total pressure shift across (0, Inf) background pressure
    """
    dz::Float64
    """
    the "turning point" pressure
    """
    pstar::Float64
    """
    the slope of the pressure-displacement response curve
    """
    alpha::Float64
end

function pressure_shift(model::LogisticPressureShift, pB::Float64, channel_length::Float64)
    (; z0, dz, alpha, pstar) = model
    p_ratio = pB / pstar
    zstar = z0 + dz / (1 + (alpha - 1)^(2 * p_ratio - 1))
    return channel_length * zstar
end

"""
    SimpleLogisticShift(model, z0, dz, pstar, alpha)
A wrapper model that allows a transport profile to shift axially in response to changes in background pressure.
As with LogisticPressureShift, the displacement/shift of the transport profile follows a logistic curve.
However, the parameterization is different, so that the shift is zero when
the background pressure is zero.
As such, it does not have a z0 parameter.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct SimpleLogisticShift{A <: AnomalousTransportModel} <: PressureShift
    """
    An AnomalousTransportModel
    """
    model::A
    """
    The maximum displacement in response to increasing pressure, scaled by the discharge channel length.
    This should be positive.
    """
    shift_length::Float64
    """
    The pressure at the midpoint of the shift, in Torr.
    Defaults to 25e-6 Torr, which gives good fits for the H9 and SPT-100.
    """
    midpoint_pressure::Float64 = 25.0e-6
    """
    The slope of the pressure response curve.
    Defaults to 2, which gives good fits for the H9 and SPT-100.
    """
    slope::Float64 = 2.0
end

function pressure_shift(model::SimpleLogisticShift, pB::Float64, channel_length::Float64)
    (; shift_length, midpoint_pressure, slope) = model
    p_ratio = pB / midpoint_pressure
    zstar = shift_length * (inv(1 + exp(-slope * (p_ratio - 1))) - inv(1 + exp(slope)))
    return -channel_length * zstar
end

"""
    num_anom_variables(::AnomalousTransportModel)::Int

The number of variable arrays that should be allocated for the provided anomalous
transport model. These arrays are used to save state beyond the anomalous
collision frequency, and are useful for defining more complex anomalous transport
models. If not defined by the user, this defaults to zero.
"""
num_anom_variables(::AnomalousTransportModel)::Int = 0
