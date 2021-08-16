γ(f::Fluid) = f.species.element.γ
m(f::Fluid) = f.species.element.m
R(f::Fluid) = f.species.element.R
cp(f::Fluid) = f.species.element.cp
cv(f::Fluid) = f.species.element.cv

number_density(U, f::Fluid) = U[1]
density(U, f::Fluid) = number_density(U, f) * m(f)

velocity(U, f::ContinuityFluid) = f.conservation_laws.u
velocity(U, f::Fluid) = U[2] / U[1]

temperature(U, f::Fluid) = f.conservation_laws.T
temperature(U, f::EulerFluid) =
    pressure(U, f) / number_density(U, f) / kB

pressure(U, f::Fluid) = number_density(U, f) * kB * temperature(U, f)
pressure(U, f::EulerFluid) =
    m(f) * (γ(f)-1) * (U[3] - 0.5 * U[2]^2/U[1])

stagnation_energy(U, f::Fluid) =
    0.5 * velocity(U, f)^2 + static_energy(U, f)
stagnation_energy(U, f::EulerFluid) = U[3] / U[1]

static_energy(U, f::Fluid) = cv(f) * temperature(U, f)
static_energy(U, f::EulerFluid) =  U[3]/U[1] - 0.5 * (U[2]/U[1])^2

static_enthalpy(U, f::Fluid) = cp(f) * temperature(U, f)
static_enthalpy(U, f::EulerFluid) =
    U[3]/U[1] - 0.5 * (U[2]/U[1])^2 + pressure(U, f) / density(U, f)

sound_speed(U, f::Fluid) = sqrt(γ(f) * R(f) * temperature(U, f))
mach_number(U, f::Fluid) = velocity(U, f) / sound_speed(U, f)

stagnation_enthalpy(U, f) =
    stagnation_energy(U, f) + pressure(U, f) / density(U, f)
critical_sound_speed(U, f) = let γ = γ(f)
    2 * (γ - 1) / (γ + 1) * stagnation_enthalpy(U, f)
end