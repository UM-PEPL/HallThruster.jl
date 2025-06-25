@inline γ(f::Fluid) = f.species.element.γ
@inline m(f::Fluid) = f.species.element.m
@inline R(f::Fluid) = f.species.element.R
@inline cp(f::Fluid) = f.species.element.cp
@inline cv(f::Fluid) = f.species.element.cv

@inline number_density(U, f::Fluid) = density(U, f) / m(f)
@inline density(U, f::Fluid) = U[1]

@inline velocity(U::NTuple{1, T}, f::Fluid) where {T} = f.u
@inline velocity(U::NTuple{2, T}, f::Fluid) where {T} = U[2] / U[1]
@inline velocity(U::NTuple{3, T}, f::Fluid) where {T} = U[2] / U[1]

@inline temperature(U::NTuple{1, T}, f::Fluid) where {T} = f.T
@inline temperature(U::NTuple{2, T}, f::Fluid) where {T} = f.T

@inline function temperature(U::NTuple{3, T}, f::Fluid) where {T}
    return (γ(f) - 1) * (U[3] - 0.5 * U[2]^2 / U[1]) / (U[1] * R(f))
end

@inline pressure(U::NTuple{1, T}, f::Fluid) where {T} = U[1] * R(f) * f.T
@inline pressure(U::NTuple{2, T}, f::Fluid) where {T} = U[1] * R(f) * f.T
@inline pressure(U::NTuple{3, T}, f::Fluid) where {T} = (γ(f) - 1) *
    (U[3] - 0.5 * U[2]^2 / U[1])

@inline sound_speed(U::NTuple{1, T}, f::Fluid) where {T} = f.a
@inline sound_speed(U::NTuple{2, T}, f::Fluid) where {T} = f.a
@inline sound_speed(U::NTuple{3, T}, f::Fluid) where {T} = sqrt(
    γ(f) * R(f) *
        temperature(U, f)
)
