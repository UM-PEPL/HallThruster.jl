Base.@kwdef struct IonizationReaction{I}
	reactant::Species
	product::Species
	rate_coeff::I
end

struct LinearInterpolation{X<:Number, Y<:Number}
	xs::Vector{X}
	ys::Vector{Y}
	function LinearInterpolation(x, y)
		if length(x) != length(y)
			throw(ArgumentError("x and y must have same length"))
		else
			return new{typeof(x[1]), typeof(y[1])}(x, y)
		end
	end
end

function (itp::LinearInterpolation)(x::T) where T
	xs, ys = itp.xs, itp.ys
	if x ≤ xs[1]
		return ys[1]/oneunit(T)
	elseif x ≥ xs[end]
		return ys[end]/oneunit(T)
	end
	i = find_left_index(x, xs)
	return ys[i] + (ys[i+1] - ys[i]) * (x - xs[i])/(xs[i+1] - xs[i])
end

function find_left_index(value, array)
	N = length(array)

	if value ≥ array[end]
		return N
	elseif value < array[1]
		return 0
	elseif value == array[1]
		return 1
	end

	left = 1
	right = N
	while true
		mid = (left + right) ÷ 2
		if value > array[mid+1]
			left = mid
		elseif value < array[mid]
			right = mid
		else
			return mid
		end
	end
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
		rate_coeff = LinearInterpolation(Te, k)
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
	reactions = IonizationReaction{LinearInterpolation{Float64, Float64}}[]
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