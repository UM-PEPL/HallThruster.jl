struct Fluid{C<:ConservationLawSystem}
    species::Species
    conservation_laws::C
end

const ContinuityFluid = Fluid{ContinuityOnly}
const IsothermalFluid = Fluid{IsothermalEuler}
const EulerFluid = Fluid{EulerEquations}

nvars(::Type{Fluid{C}}) where C = nvars(C)

function ranges(fluids)
	fluid_ranges = fill(1:1, length(fluids))

	start_ind = 1
	last_ind = 1
	for (i, f) in enumerate(fluids)
		nf = nvars(f.conservation_laws |> typeof)
		last_ind = start_ind - 1 + nf
		fluid_ranges[i] = start_ind : last_ind
		start_ind = last_ind + 1
	end

	return fluid_ranges
end