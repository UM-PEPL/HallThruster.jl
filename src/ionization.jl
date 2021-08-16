Base.@kwdef struct IonizationReaction{I}
	reactant::Species
	product::Species
	rate_coeff::I
end

function Base.show(io::IO, i::IonizationReaction)
	electron_input = "e-"
	electron_output = string(i.product.Z - i.reactant.Z) * "e-"
	reactant_str = string(i.reactant)
	product_str = string(i.product)
	rxn_str = electron_input * " + " * reactant_str * " -> "
	rxn_str *= electron_output * " + " * product_str
	print(io, rxn_str)
end
