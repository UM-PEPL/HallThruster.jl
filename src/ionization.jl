Base.@kwdef struct IonizationReaction{I}
	reactant::Species
	product::Species
	rate_coeff::I
end

function Base.show(io::IO, i::IonizationReaction)
	electron_input = "e-"
	electron_output = string(i.product.Z - i.reactant.Z + 1) * "e-"
	reactant_str = string(i.reactant)
	product_str = string(i.product)
	rxn_str = electron_input * " + " * reactant_str * " -> "
	rxn_str *= electron_output * " + " * product_str
	print(io, rxn_str)
end

function load_ionization_reaction(reactant, product)
	rates_file = rate_coeff_filename(reactant, product, "ionization")
	rates_file = joinpath(REACTION_FOLDER, rates_file)
	try
		rates = DataFrame(CSV.File(rates_file))
		Te = rates[!, 1]
		k = rates[!, 2]
		rate_coeff = LinearInterpolation(Te, k, extrapolation_bc = Flat())
		return IonizationReaction(reactant, product, rate_coeff)
	catch e
		return nothing
	end
end

function rate_coeff_filename(reactant, product, reaction_type)
	return join([reaction_type, repr(reactant), repr(product)], "_") * ".dat"
end

function load_ionization_reactions(species)
	species_sorted = sort(species, by = x -> x.Z)
	reactions = IonizationReaction[]
	for i in 1:length(species)
		for j in i+1:length(species)
			reaction = load_ionization_reaction(species_sorted[i], species_sorted[j])
			if !isnothing(reaction)
				push!(reactions, reaction)
			end
		end
	end

	return reactions
end