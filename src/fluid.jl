struct Fluid{C<:ConservationLawSystem}
    species::Species
    conservation_laws::C
end

const ContinuityFluid = Fluid{ContinuityOnly}
const IsothermalFluid = Fluid{IsothermalEuler}
const EulerFluid = Fluid{EulerEquations}

nvars(::Type{Fluid{C}}) where C = nvars(C)

