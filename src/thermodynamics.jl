@inline γ(f::Fluid) = f.species.element.γ
@inline m(f::Fluid) = f.species.element.m
@inline R(f::Fluid) = f.species.element.R
@inline cp(f::Fluid) = f.species.element.cp
@inline cv(f::Fluid) = f.species.element.cv

@inline number_density(U, f::Fluid) = density(U, f) / m(f)
@inline density(U, f::Fluid) = U[1]

function compute_primitive_continuity(U, f::Fluid)
    ρ = U[1]
    u = f.conservation_laws.u
    p = ρ * R(f) * f.conservation_laws.T
    return ρ, u, p
end

function compute_primitive_isothermal(U, f::Fluid)
    ρ = U[1]
    u = U[2] / U[1]
    p = U[1] * R(f) * f.conservation_laws.T
    return ρ, u, p
end

function compute_primitive_euler(U, f::Fluid)
    γ = fluid.species.element.γ
    ρ, ρu, ρE = U
    u = ρu / ρ
    p = (γ - 1) * (ρE - 0.5 * ρu * u)
    return ρ, u, p
end

function compute_primitive(U, f::Fluid)
    ρ = U[1]
    if f.conservation_laws.type == _ContinuityOnly
        u = f.conservation_laws.u
    else
        u = U[2] / U[1]
    end

    if f.conservation_laws.type == _EulerEquations
        ρu = U[2]
        ρE = U[3]
        p = (γ - 1) * (ρE - 0.5 * ρu * u)
    else
        p = U[1] * R(f) * f.conservation_laws.T
    end

    return ρ, u, p
end

function compute_conservative(ρ, u, p, γ)
    ρE = p / (γ - 1) + 0.5 * ρ * u^2
    return ρ, ρ * u, ρE
end

@inline velocity(U::SVector{1, T}, f::Fluid) where T = f.conservation_laws.u
@inline velocity(U::SVector{2, T}, f::Fluid) where T = U[2] / U[1]
@inline velocity(U::SVector{3, T}, f::Fluid) where T = U[2] / U[1]

@inline function velocity(U, f::Fluid)
    if f.conservation_laws.type == _ContinuityOnly
        return f.conservation_laws.u
    else
        return U[2] / U[1]
    end
end

@inline function temperature(U, f::Fluid)
    if f.conservation_laws.type == _EulerEquations
        (γ(f) - 1) * (U[3] - 0.5 * U[2]^2 / U[1]) / U[1] / R(f)
    else
        return f.conservation_laws.T
    end
end

@inline function pressure(U, f::Fluid)
    if f.conservation_laws.type == _EulerEquations
        return (γ(f) - 1) * (U[3] - 0.5 * U[2]^2 / U[1])
    else
        return U[1] * R(f) * f.conservation_laws.T
    end
end

function stagnation_energy(U, f::Fluid)
    if f.conservation_laws.type == _EulerEquations
        return U[3] + (0.5 * (U[2])^2 / (U[1])^2) * (1 - U[1])
    else
        return 0.5 * velocity(U, f)^2 + static_energy(U, f)
    end
end

function static_energy(U, f::Fluid)
    if f.conservation_laws.type == _EulerEquations
        return U[3] / U[1] - 0.5 * (U[2] / U[1])^2
    else
        return cv(f) * temperature(U, f)
    end
end

function static_enthalpy(U, f::Fluid)
    if f.conservation_laws.type == _EulerEquations
        return U[3] / U[1] - 0.5 * (U[2] / U[1])^2 + pressure(U, f) / density(U, f)
    else
        return cp(f) * temperature(U, f)
    end
end

@inline sound_speed(U, f::Fluid) = sqrt(γ(f) * R(f) * temperature(U, f))
@inline mach_number(U, f::Fluid) = velocity(U, f) / sound_speed(U, f)

@inline stagnation_enthalpy(U, f) = stagnation_energy(U, f) + pressure(U, f) / density(U, f)
@inline function critical_sound_speed(U, f)
    return 2 * (γ(f) - 1) / (γ(f) + 1) * stagnation_enthalpy(U, f)
end
