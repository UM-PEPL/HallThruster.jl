struct Fluid
    species::Species
    conservation_laws::ConservationLawSystem
end

nvars(f::Fluid) = nvars(f.conservation_laws)

function ranges(fluids)
	fluid_ranges = fill(1:1, length(fluids))

	start_ind = 1
	last_ind = 1
	for (i, f) in enumerate(fluids)
		nf = nvars(f)
		last_ind = start_ind - 1 + nf
		fluid_ranges[i] = start_ind : last_ind
		start_ind = last_ind + 1
	end

	return fluid_ranges
end