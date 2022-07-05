"""
    AnomalousTransportModel
The abstract supertype of all types of anomalous transport model. Currently, there are two subtypes: `FixedAnomModel` and `ZeroEquationModel`, which
correspond, respectively, to models which are initialized and then never updated, and models which depend algebraically on local plasma properties.
"""
abstract type AnomalousTransportModel end

"""
    FixedAnomModel <: AnomalousTransportModel
An anomalous collision frequency model which is fixed at the start of a simulation and then never updated again.
This is the most common type of model, and includes the ubiquitous multi-zone bohm-like collision frequency models.

Users defining their own `FixedAnomModel` must implement the `initialize_anom!(νan, model, U, params)` method for their type.
"""
abstract type FixedAnomModel <: AnomalousTransportModel end

"""
    ZeroEquationModel <: AnomalousTransportModel
An anomalous collision frequency which depends algebraically on the local plasma properties and is updated at each iteration.

Users defining their own `ZeroEquationModel` must implement the `evaluate_anom(model, U, params)` method for their type.
"""
abstract type ZeroEquationModel <: AnomalousTransportModel end

"""
    initialize_anom!(νan, model::FixedAnomModel, U, params)
Initialize the anomalous collision frequency profile for models in which the anomalous mobility does not change with time. Must be
implemented by the user for specific subtypes of FixedAnomModel. Should update argument νan in-place with the desired value of anomalous
collision frequency and then return the mutated array.
"""
function initialize_anom!(νan, model::FixedAnomModel, U, params)
    throw(ArgumentError(
        "Function `initialize_anom!(νan, model, U, params)` is not implemented for model of type $(typeof(model)). This function must be implemented by the user for `FixedAnomModel` objects."
    ))
end

"""
    evaluate_anom(model::FixedAnomModel, U, params, i)
By default, this does nothing for FixedAnomModel objects, as the anomalous collision frequency is not updated at each iteration.
"""
function evaluate_anom(model::FixedAnomModel, U, params, i)
    return params.cache.νan[i]
end


"""
    initialize_anom!(νan, model::ZeroEquationModel, U, params)
Default anomalous collision frequency initialization for `ZeroEquationModel`s. Just calls evalute_anom at each grid location.
"""
function initialize_anom!(νan, model::ZeroEquationModel, U, params)
    for i in eachindex(z_cell)
        νan[i] = evaluate_anom(model, U, params, i)
    end
    return νan
end

"""
    evaluate_anom(model::ZeroEquationModel, U, params, i)
Return the anomalous collision frequency at index `i` for a provided ZeroEquationModel. Must be implemented by the user.
"""
function evaluate_anom(model::ZeroEquationModel, U, params, i)
    throw(ArgumentError(
        "Function `evaluate_anom(model, U, params, i)` is not implemented for model of type $(typeof(model)). This function must be implemented by the user for `ZeroEquationModel` objects."
    ))
end


#=============================================================================
 Begin definition of built-in models of type FixedAnomModel
==============================================================================#

"""
    NoAnom <: FixedAnomModel
No anomalous collision frequency included in simulation
"""
struct NoAnom <: FixedAnomModel end

function initialize_anom!(νan, model::NoAnom, U, params)
    νan .= 0.0
    return νan
end

"""
    Bohm(c) <: FixedAnomModel
Model where the anomalous collision frequency scales with the electron cyclotron frequency ωce times some scaling factor c
"""
struct Bohm <: FixedAnomModel
    c::Float64
end

function initialize_anom!(νan, model::Bohm, U, params)
    for i in eachindex(νan)
        B = params.cache.B[i]
        ωce = e * B / me
        νan[i] = model.c * ωce
    end

    return νan
end

"""
    TwoZoneBohm(c1, c2) <: FixedAnomModel
Model where the anomalous collision frequency has two values: c1 * ωce inside the channel and c2 * ωce outside of the channel.
Takes two arguments: c1 and c2. The transition between these values can be smoothed by the user-provided transition function.
"""
struct TwoZoneBohm <: FixedAnomModel
    coeffs::NTuple{2, Float64}
    TwoZoneBohm(c1, c2) = new((c1, c2))
end

function initialize_anom!(νan, model::TwoZoneBohm, U, params)
    c1, c2 = model.coeffs

    for i in eachindex(νan)
        B = params.cache.B[i]
        ωce = e * B / me

        β = params.config.transition_function(params.z_cell[i], params.L_ch, c1, c2)
        νan[i] = β * ωce
    end

    return νan
end

"""
    MultiLogBohm(zs, cs) <: FixedAnomModel
Model similar to that employed in Hall2De, where the mobility is Bohm-like (i.e. `νan(z) = c(z) * ωce(z)`) and z is in meters.

The function `c(z)` is defined by a sequence of nodes `(z, c)` provided by the user. At `z = z[1]`, `c(z) = c[1]`, and so forth.

At `z[i] < z < z[i+1]`, `log(c)` is defined by linearly interpolating between `log(c[i])` and `log(c[i+1])`.

For `z < z[1]`, `c = c[1]` and for `z > z[end]`, `c(z) = c[end]`.

The user may also provide a single array of [z[1], z[2], ..., z[end], c[1], c[2], ..., c[end]]. The number of z values must be equal to the number of c values.
"""
struct MultiLogBohm <: FixedAnomModel
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
    cs = coeffs[N+1:end]
    return MultiLogBohm(zs, cs)
end

function initialize_anom!(νan, model::MultiLogBohm, U, params)
    (;z_cell) = params
    (;B) = params.cache

    for i in eachindex(νan)
        ωce = e * B[i] / me
        c = HallThruster.interpolate(z_cell[i], model.zs, model.cs, use_log = true)
        νan[i] = c * ωce
    end

    return νan
end

#=============================================================================
 Begin definition of built-in models of type ZeroEquationModel
==============================================================================#

# None yet, but users may define their own.
