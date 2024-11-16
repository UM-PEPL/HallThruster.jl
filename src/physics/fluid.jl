@enum ConservationLawType begin
    _ContinuityOnly
    _IsothermalEuler
    _EulerEquations
end

Base.@kwdef struct Fluid
    species::Species
    type::ConservationLawType
    nvars::Int
    u::Union{Float64,Nothing} = nothing
    T::Union{Float64,Nothing} = nothing
    a::Union{Float64,Nothing} = nothing
end

@inline nvars(f::Fluid) = f.nvars

function Fluid(s; u=nothing, T=nothing)
    if isnothing(u) && isnothing(T)
        return EulerEquations(s)
    elseif isnothing(u)
        return Fluid(s, T)
    else
        return Fluid(s, u, T)
    end
end

function ContinuityOnly(s; u, T)
    a = √(s.element.γ * kB * T / s.element.m)
    return Fluid(;
        species=s, type=_ContinuityOnly, nvars=1, u=Float64(u), T=Float64(T), a=Float64(a)
    )
end
ContinuityOnly(s, u, T) = ContinuityOnly(s; u, T)
Fluid(s, u, T) = ContinuityOnly(s, u, T)

function IsothermalEuler(s; T)
    a = √(s.element.γ * kB * T / s.element.m)
    return Fluid(; species=s, type=_IsothermalEuler, nvars=2, T=Float64(T), a=Float64(a))
end
IsothermalEuler(s, T) = IsothermalEuler(s; T)
Fluid(s, T) = IsothermalEuler(s, T)

function EulerEquations(s)
    return Fluid(; species=s, type=_EulerEquations, nvars=3)
end

function ranges(fluids)
    fluid_ranges = fill(1:1, length(fluids))
    start_ind = 1
    last_ind = 1
    for (i, f) in enumerate(fluids)
        nf = nvars(f)
        last_ind = start_ind - 1 + nf
        fluid_ranges[i] = start_ind:last_ind
        start_ind = last_ind + 1
    end

    return fluid_ranges
end
