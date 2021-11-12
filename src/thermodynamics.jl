γ(f::Fluid) = f.species.element.γ
m(f::Fluid) = f.species.element.m
R(f::Fluid) = f.species.element.R
cp(f::Fluid) = f.species.element.cp
cv(f::Fluid) = f.species.element.cv

number_density(U, f::Fluid) = density(U, f) / m(f)
density(U, f::Fluid) = U[1]

function velocity(U, f::Fluid)
    if f.conservation_laws.type == :ContinuityOnly
        return f.conservation_laws.u
    else
        return U[2] / U[1]
    end
end

function temperature(U, f::Fluid)
    if f.conservation_laws.type == :EulerEquations
        pressure(U, f) / density(U, f) / R(f)
    else
        return f.conservation_laws.T
    end
end

function pressure(U, f::Fluid)
    if f.conservation_laws.type == :EulerEquations
        return (γ(f)-1) * (U[3] - 0.5 * U[2]^2/U[1])
    else
        return density(U, f) * R(f) * temperature(U, f)
    end
end

function stagnation_energy(U, f::Fluid)
    if f.conservation_laws.type == :EulerEquations
        return  U[3] + (0.5*(U[2])^2/(U[1])^2)*(1 - U[1])
    else
        return 0.5 * velocity(U, f)^2 + static_energy(U, f)
    end
end

function static_energy(U, f::Fluid)
    if f.conservation_laws.type == :EulerEquations
        return U[3]/U[1] - 0.5 * (U[2]/U[1])^2
    else
        return cv(f) * temperature(U, f)
    end
end

function static_enthalpy(U, f::Fluid)
    if f.conservation_laws.type == :EulerEquations
        return U[3]/U[1] - 0.5 * (U[2]/U[1])^2 + pressure(U, f) / density(U, f)
    else
        return cp(f) * temperature(U, f)
    end
end

sound_speed(U, f::Fluid) = sqrt(γ(f) * R(f) * temperature(U, f))
mach_number(U, f::Fluid) = velocity(U, f) / sound_speed(U, f)

stagnation_enthalpy(U, f) =
    stagnation_energy(U, f) + pressure(U, f) / density(U, f)
critical_sound_speed(U, f) = let γ = γ(f)
    2 * (γ - 1) / (γ + 1) * stagnation_enthalpy(U, f)
end
