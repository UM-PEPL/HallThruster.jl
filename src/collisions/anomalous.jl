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
        StepTroughBohm1,
        StepTroughBohm2,
    )
end

function Serialization.SType(::Type{T}) where {T <: AnomalousTransportModel}
    return Serialization.TaggedUnion()
end
Serialization.options(::Type{T}) where {T <: AnomalousTransportModel} = anom_models()

#=============================================================================
  Validation functions
==============================================================================#
function check_finite(var, var_name)
    isfinite(var) || throw(ArgumentError("`$(var_name)` must be finite. Got: $(var)."))
    return nothing
end
function check_positive(var, var_name)
    (isfinite(var) && var > 0) || throw(ArgumentError("`$(var_name)` must be positive. Got: $(var)."))
    return nothing
end
function check_nonnegative(var, var_name)
    (isfinite(var) && var >= 0) || throw(ArgumentError("`$(var_name)` must be nonnegative. Got: $(var)."))
    return nothing
end
function check_in_interval(var, var_name, lb, ub; upper_inclusive = true)
    valid_upper = upper_inclusive ? var <= ub : var < ub
    if !(isfinite(var) && var >= lb && valid_upper)
        right_bracket = upper_inclusive ? "]" : ")"
        throw(ArgumentError("`$(var_name)` must be in the interval [$(lb), $(ub)$(right_bracket). Got: $(var)."))
    end
    return nothing
end

macro check_finite(var)
    return :(check_finite($(esc(var)), $(string(var))))
end

macro check_positive(var)
    return :(check_positive($(esc(var)), $(string(var))))
end

macro check_nonnegative(var)
    return :(check_nonnegative($(esc(var)), $(string(var))))
end

macro check_in_interval(var, lb, ub, upper_inclusive = true)
    return :(
        check_in_interval(
            $(esc(var)), $(string(var)), $(esc(lb)), $(esc(ub));
            upper_inclusive = $(esc(upper_inclusive)),
        )
    )
end

#=============================================================================
 Definition of built-in models
==============================================================================#

"""
    NoAnom <: AnomalousTransportModel
No anomalous collision frequency included in simulation
"""
struct NoAnom <: AnomalousTransportModel end

function (::NoAnom)(νan, @nospecialize(_params), @nospecialize(_z_shift::Float64 = 0.0))
    @inbounds for i in eachindex(νan)
        νan[i] = 0.0
    end
    return νan
end

"""
    Bohm(c) <: AnomalousTransportModel
Model where the anomalous collision frequency scales with the electron cyclotron frequency
as `νan = c * ωce`.

# Fields
$(TYPEDFIELDS)
"""
struct Bohm <: AnomalousTransportModel
    """Nonnegative inverse Hall parameter."""
    c::Float64
    function Bohm(c)
        c = Float64(c)
        @check_nonnegative c
        return new(c)
    end
end
Bohm(; c) = Bohm(c)

function (model::Bohm)(νan, params, ::Float64 = 0.0)
    (; cache, grid) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    B_interp = LinearInterpolation(grid.cell_centers, B)

    for (i, zc) in enumerate(grid.cell_centers)
        B = B_interp(zc)
        ωce = e * B / me
        νan[i] = model.c * ωce
    end

    return νan
end

"""
    TwoZoneBohm(c1, c2) <: AnomalousTransportModel
Model where the anomalous collision frequency has two values: `c1 * ωce` inside
the channel and `c2 * ωce` outside it. The transition between these values is
smoothed over `params.transition_length`.

# Fields
$(TYPEDFIELDS)
"""
struct TwoZoneBohm <: AnomalousTransportModel
    """Nonnegative inverse Hall parameter inside the channel."""
    c1::Float64
    """Nonnegative inverse Hall parameter outside the channel."""
    c2::Float64
    function TwoZoneBohm(c1, c2)
        c1 = Float64(c1)
        c2 = Float64(c2)
        @check_nonnegative c1
        @check_nonnegative c2
        return new(c1, c2)
    end
end
TwoZoneBohm(; c1, c2) = TwoZoneBohm(c1, c2)

function (model::TwoZoneBohm)(νan, params, z_shift::Float64 = 0.0)
    (; c1, c2) = model
    (; cache, grid, thruster) = params
    (; B) = cache

    L_trans = params.transition_length
    L_ch = params.thruster.geometry.channel_length

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
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
Model similar to that employed in Hall2De, where the anomalous collision frequency
is Bohm-like (i.e. `νan(z) = c(z) * ωce(z)`) and `z` is in meters.

The function `c(z)` is defined by a sequence of nodes `(z, c)` provided by the user. At `z = z[1]`, `c(z) = c[1]`, and so forth.

At `z[i] < z < z[i+1]`, `log(c)` is defined by linearly interpolating between `log(c[i])` and `log(c[i+1])`.

For `z < z[1]`, `c = c[1]` and for `z > z[end]`, `c(z) = c[end]`.

The `zs` values must be finite and strictly increasing. The `cs` values must be
finite and positive because their logarithms are interpolated. Both arrays must
be nonempty and have the same length.

# Fields
$(TYPEDFIELDS)
"""
struct MultiLogBohm <: AnomalousTransportModel
    """Finite, strictly increasing axial node positions in meters."""
    zs::Vector{Float64}
    """Positive inverse Hall parameters at the axial nodes."""
    cs::Vector{Float64}
    function MultiLogBohm(zs, cs)
        zs = Float64.(zs)
        cs = Float64.(cs)
        isempty(zs) && throw(ArgumentError("`zs` and `cs` must be nonempty."))
        length(zs) == length(cs) || throw(ArgumentError("Number of z values must be equal to number of c values."))
        all(isfinite, zs) || throw(ArgumentError("All `zs` values must be finite."))
        all(>(0), diff(zs)) || throw(ArgumentError("`zs` values must be strictly increasing."))
        all(c -> isfinite(c) && c > 0, cs) || throw(ArgumentError("All `cs` values must be positive and finite."))
        return new(zs, cs)
    end
end
MultiLogBohm(; zs, cs) = MultiLogBohm(zs, cs)

function (model::MultiLogBohm)(νan, params, z_shift::Float64 = 0.0)
    (; grid) = params
    (; B) = params.cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

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
Model in which the anomalous collision frequency is Bohm-like (`νan ~ ωce`),
with a Gaussian trough centered at `center`. The inverse Hall parameter is
`hall_max` far from the trough and `hall_min * hall_max` at its center.

# Fields
$(TYPEDFIELDS)
"""
struct GaussianBohm <: AnomalousTransportModel
    """Fraction of `hall_max` retained at the trough center; must be in [0, 1]."""
    hall_min::Float64
    """Inverse Hall parameter far from the trough; must be positive."""
    hall_max::Float64
    """Finite axial position of the center of the Gaussian trough, in meters."""
    center::Float64
    """Positive standard deviation of the Gaussian trough, in meters."""
    width::Float64
    function GaussianBohm(hall_min, hall_max, center, width)
        hall_min = Float64(hall_min)
        hall_max = Float64(hall_max)
        center = Float64(center)
        width = Float64(width)
        @check_in_interval hall_min 0 1
        @check_positive hall_max
        @check_finite center
        @check_positive width
        return new(hall_min, hall_max, center, width)
    end
end
GaussianBohm(; hall_min, hall_max, center, width) = GaussianBohm(hall_min, hall_max, center, width)

function (model::GaussianBohm)(νan, params, z_shift::Float64 = 0.0)
    (; hall_min, hall_max, center, width) = model
    (; cache, grid) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

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
struct ScaledGaussianBohm <: AnomalousTransportModel
    """The maximum inverse hall parameter. Must be positive."""
    anom_scale::Float64
    """The factor by which transport is reduced by the baseline value at the center of the trough, must be in [0,1]. """
    barrier_scale::Float64
    """The standard deviation of the Gaussian trough, in channel lengths. Must be positive."""
    width::Float64
    """The axial position of the mean of the Gaussian trough, in channel lengths. Must be positive."""
    center::Float64
    function ScaledGaussianBohm(anom_scale, barrier_scale, width, center)
        anom_scale = Float64(anom_scale)
        barrier_scale = Float64(barrier_scale)
        width = Float64(width)
        center = Float64(center)
        @check_positive anom_scale
        @check_in_interval barrier_scale 0 1
        @check_positive width
        @check_positive center
        return new(anom_scale, barrier_scale, width, center)
    end
end
function ScaledGaussianBohm(; anom_scale = 0.0625, barrier_scale = 0.9, width, center)
    return ScaledGaussianBohm(anom_scale, barrier_scale, width, center)
end

function (model::ScaledGaussianBohm)(νan, params, z_shift::Float64 = 0.0)
    (; anom_scale, barrier_scale, width, center) = model
    (; cache, grid, thruster) = params
    (; B) = cache

    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end

    L_ch = thruster.geometry.channel_length
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

"""
    StepTroughBohm1(anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent) <: AnomalousTransportModel
Model in which the anomalous collision frequency is Bohm-like (`νan ~ ω_ce`),
with a shape function defined by a product of a logistic "step" function and an inverted generalized Gaussian near the exit plane.
Parameterization is similar to the ScaledGaussianBohm, with all quantities chosen to be O(1).


# Fields
$(TYPEDFIELDS)
"""
struct StepTroughBohm1 <: AnomalousTransportModel
    """The maximum inverse Hall parameter; must be positive."""
    anom_scale::Float64
    """The axial position of the co-located center of the generalized Gaussian trough and logistic step, in channel lengths. Must be positive."""
    anom_center::Float64
    """The size of the logistic step. Zero removes the step; one makes its upstream limit zero. Must be in [0, 1]."""
    step_scale::Float64
    """Dimensionless logistic sharpness parameter. The step sharpens as this approaches zero. Must be in (0, 1)."""
    step_width::Float64
    """Fraction of the baseline transport retained at the trough center. Must be in [0, 1]."""
    trough_floor::Float64
    """Positive generalized-Gaussian width relative to `anom_center`. Must be positive."""
    trough_width::Float64
    """Scaled generalized-Gaussian exponent. It maps 0 to 1 and 0.5 to 2, and approaches infinity as it approaches 1. Must be in [0, 1)."""
    trough_exponent::Float64

    function StepTroughBohm1(anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent)
        anom_scale = Float64(anom_scale)
        anom_center = Float64(anom_center)
        step_scale = Float64(step_scale)
        step_width = Float64(step_width)
        trough_floor = Float64(trough_floor)
        trough_width = Float64(trough_width)
        trough_exponent = Float64(trough_exponent)
        @check_positive anom_scale
        @check_positive anom_center
        @check_in_interval step_scale 0 1
        @check_in_interval step_width 0 1 false
        @check_positive step_width
        @check_in_interval trough_floor 0 1
        @check_positive trough_width
        @check_in_interval trough_exponent 0 1 false
        return new(anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent)
    end
end

function StepTroughBohm1(; anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent)
    return StepTroughBohm1(anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent)
end

function (model::StepTroughBohm1)(νan::Vector{Float64}, z::Vector{Float64}, B::Vector{Float64}, L_ch::Float64 = 1.0, z_shift::Float64 = 0.0)
    (; anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent) = model

    @inbounds for i in eachindex(νan)
        z0 = (z[i] - z_shift) / L_ch
        z_aux = (z0 / anom_center) - 1

        # logistic part
        step = (1 - step_scale) + step_scale / (1 + exp(-z_aux * (1 - step_width) / step_width))

        # exponential part
        trough = 1 - (1 - trough_floor) * exp(-abs(z_aux / trough_width)^(1 / (1 - trough_exponent)))

        # result
        inverse_hall = anom_scale * step * trough
        ωce = e * B[i] / me
        νan[i] = inverse_hall * ωce
    end
    return νan
end

function (model::StepTroughBohm1)(νan::Vector{Float64}, params, z_shift::Float64 = 0.0)
    (; cache, grid, thruster) = params
    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end
    return model(νan, grid.cell_centers, cache.B, thruster.geometry.channel_length, z_shift)
end

"""
    StepTroughBohm2(anom_scale, anom_center, step_scale, step_width, trough_floor, trough_width, trough_exponent) <: AnomalousTransportModel
Similar concept to StepTroughBohm1, but the step is a smoothstep and the trough is a bump function with compact support.
These ensure that the trough has a finite extent in space.


# Fields
$(TYPEDFIELDS)
"""
struct StepTroughBohm2 <: AnomalousTransportModel
    """The maximum inverse Hall parameter in the plume. Must be positive."""
    anom_scale::Float64
    """The axial position of the co-located center of the generalized Gaussian trough and logistic step, in channel length. Must be positive."""
    anom_center::Float64
    """Combined width parameter for the step and trough. Must be positive."""
    anom_width::Float64
    """The ratio of the anomalous collision frequency at the anode to that downstream. Must be positive."""
    anode_scale::Float64
    """Fraction of the baseline transport retained at the trough center. Must be in [0, 1)."""
    trough_floor::Float64
    """Controls how round the shoulders of the trough are. Must be in [0, 1]."""
    trough_roundness::Float64
    """Scaled generalized-Gaussian exponent. It maps 0 to 1 and 0.5 to 2, and approaches infinity as it approaches 1. Must be in [0, 1)."""
    trough_exponent::Float64

    function StepTroughBohm2(anom_scale, anom_center, anom_width, anode_scale, trough_floor, trough_roundness, trough_exponent)
        anom_scale = Float64(anom_scale)
        anom_center = Float64(anom_center)
        anom_width = Float64(anom_width)
        anode_scale = Float64(anode_scale)
        trough_floor = Float64(trough_floor)
        trough_roundness = Float64(trough_roundness)
        trough_exponent = Float64(trough_exponent)
        @check_positive anom_scale
        @check_positive anom_center
        @check_positive anode_scale
        @check_positive anom_width
        @check_positive anom_width
        @check_in_interval trough_floor 0 1
        @check_in_interval trough_roundness 0 1 false
        @check_in_interval trough_exponent 0 1 false
        return new(anom_scale, anom_center, anom_width, anode_scale, trough_floor, trough_roundness, trough_exponent)
    end
end

function StepTroughBohm2(; anom_scale, anom_center, anom_width, anode_scale, trough_floor, trough_roundness, trough_exponent)
    return StepTroughBohm2(anom_scale, anom_center, anom_width, anode_scale, trough_floor, trough_roundness, trough_exponent)
end

function smootherstep(x::T) where T <: AbstractFloat
    return x <= 0 ? zero(T) : x >= 1 ? oneunit(T) : 6 * x^5 - 15 * x^4 + 10 * x^3
end

function bump_function(x::T, roundness, exponent) where T <: AbstractFloat
    p = 1 / (1 - exponent)
    r = 1 / (1 - roundness)
    return abs(x) >= 1 ? zero(T) : (1 - abs(x)^p)^r
end

function (model::StepTroughBohm2)(νan::Vector{Float64}, z::Vector{Float64}, B::Vector{Float64}, L_ch::Float64 = 1.0, z_shift::Float64 = 0.0)
    S = model.anom_scale
    L = model.anom_center
    w = model.anom_width
    s = model.anode_scale
    a = 1 - model.trough_floor
    p = model.trough_exponent
    r = model.trough_roundness

    @inbounds for i in eachindex(νan)
        z0 = (z[i] - z_shift) / L_ch
        z_norm = (z0 - L) / w

        # Bump (trough)
        bump = bump_function(z_norm, model.trough_exponent, model.trough_roundness)

        # Smootherstep part
        step = S * (s + (1 - s) * smootherstep(0.5 * (z_norm + 1)))

        # result
        inverse_hall = step * (1 - a * bump)
        ωce = e * B[i] / me
        νan[i] = inverse_hall * ωce
    end
    return νan
end

function (model::StepTroughBohm2)(νan::Vector{Float64}, params, z_shift::Float64 = 0.0)
    (; cache, grid, thruster) = params
    # Profile is fixed in time, do not update after 5 iterations
    if (params.iteration[] > 5)
        return νan
    end
    return model(νan, grid.cell_centers, cache.B, thruster.geometry.channel_length, z_shift)
end


#=============================================================================
 Begin definition of built-in models
==============================================================================#

abstract type PressureShift <: AnomalousTransportModel end

pressure_shift(model::AnomalousTransportModel, ::Any, ::Any) = 0.0

function (model::PressureShift)(νan, params, _z::Float64 = 0.0)
    z_shift = pressure_shift(model, params.background_pressure_Torr, params.thruster.geometry.channel_length)
    return model.model(νan, params, z_shift)
end

"""
    LogisticPressureShift(model, z0, dz, pstar, alpha)
A wrapper model that allows a transport profile to shift axially in response to changes in background pressure.
The displacement/shift of the transport profile follows a logistic curve.

# Fields
$(TYPEDFIELDS)
"""
struct LogisticPressureShift{A <: AnomalousTransportModel} <: PressureShift
    """
    An anomalous transport model
    """
    model::A
    """
    Dimensionless shift offset, scaled by the channel length. Must be finite.
    """
    z0::Float64
    """
    The shift amplitude, scaled by the channel length. Must be finite.
    """
    dz::Float64
    """
    The positive pressure scale in Torr.
    """
    pstar::Float64
    """
    Shape parameter for the pressure-displacement response curve; must be greater than 1.
    """
    alpha::Float64
    function LogisticPressureShift(model::A, z0, dz, pstar, alpha) where {A <: AnomalousTransportModel}
        z0 = Float64(z0)
        dz = Float64(dz)
        pstar = Float64(pstar)
        alpha = Float64(alpha)
        @check_finite z0
        @check_finite dz
        @check_positive pstar
        (isfinite(alpha) && alpha > 1) || throw(ArgumentError("`alpha` must be finite and greater than 1. Got: $(alpha)."))
        return new{A}(model, z0, dz, pstar, alpha)
    end
end
function LogisticPressureShift(; model, z0, dz, pstar, alpha)
    return LogisticPressureShift(model, z0, dz, pstar, alpha)
end

function pressure_shift(model::LogisticPressureShift, pB::Float64, channel_length::Float64)
    (; z0, dz, alpha, pstar) = model
    p_ratio = pB / pstar
    zstar = z0 + dz / (1 + (alpha - 1)^(2 * p_ratio - 1))
    return channel_length * zstar
end

"""
    SimpleLogisticShift(model, shift_length, midpoint_pressure, slope)
A wrapper model that allows a transport profile to shift axially in response to changes in background pressure.
As with LogisticPressureShift, the displacement/shift of the transport profile follows a logistic curve.
However, the parameterization is different, so that the shift is zero when
the background pressure is zero.
As such, it does not have a z0 parameter.

# Fields
$(TYPEDFIELDS)
"""
struct SimpleLogisticShift{A <: AnomalousTransportModel} <: PressureShift
    """
    An AnomalousTransportModel
    """
    model::A
    """
    Scale of the upstream displacement in response to increasing pressure, relative to
    the discharge channel length. The asymptotic displacement magnitude is
    `shift_length / (1 + exp(-slope))`. Must be positive.
    """
    shift_length::Float64
    """
    The pressure at the midpoint of the shift, in Torr.
    Defaults to 25e-6 Torr, which gives good fits for the H9 and SPT-100.
    """
    midpoint_pressure::Float64
    """
    The slope of the pressure response curve.
    Defaults to 2, which gives good fits for the H9 and SPT-100.
    """
    slope::Float64
    function SimpleLogisticShift(model::A, shift_length, midpoint_pressure, slope) where {A <: AnomalousTransportModel}
        shift_length = Float64(shift_length)
        midpoint_pressure = Float64(midpoint_pressure)
        slope = Float64(slope)
        @check_positive shift_length
        @check_positive midpoint_pressure
        @check_positive slope
        return new{A}(model, shift_length, midpoint_pressure, slope)
    end
end
function SimpleLogisticShift(; model, shift_length, midpoint_pressure = 25.0e-6, slope = 2.0)
    return SimpleLogisticShift(model, shift_length, midpoint_pressure, slope)
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
