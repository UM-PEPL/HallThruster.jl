Base.@kwdef struct IonizationReaction{I}
    reactant::Species
    product::Species
    rate_coeff::I
end

Base.@kwdef struct LandmarkTable{L, IL}
    loss_coeff::L
    rate_coeff::IL
end

function Base.show(io::IO, i::IonizationReaction)
    electron_input = "e-"
    electron_output = string(i.product.Z - i.reactant.Z + 1) * "e-"
    reactant_str = string(i.reactant)
    product_str = string(i.product)
    rxn_str = electron_input * " + " * reactant_str * " -> "
    rxn_str *= electron_output * " + " * product_str
    return print(io, rxn_str)
end

function load_ionization_reaction(reactant, product)
    rates_file = rate_coeff_filename(reactant, product, "ionization")
    rates_file = joinpath(REACTION_FOLDER, rates_file)
    rates = DataFrame(CSV.File(rates_file))
    ϵ = rates[!, 1]
    k = rates[!, 2]
    rate_coeff = LinearInterpolation(ϵ, k)
    return IonizationReaction(reactant, product, rate_coeff)
end

function load_landmark()
    rates_file = joinpath(LANDMARK_FOLDER, "landmark_rates.csv")
    rates = DataFrame(CSV.File(rates_file))
    rates = DataFrame(CSV.File(rates_file))
    ϵ = rates[!, 1]
    k_ionization = rates[!, 2]
    k_loss = rates[!, 3]
    rate_coeff = LinearInterpolation(ϵ, k_ionization)
    loss_coeff = LinearInterpolation(ϵ, k_loss)
    return LandmarkTable(loss_coeff, rate_coeff)
end

function rate_coeff_filename(reactant, product, reaction_type)
    return join([reaction_type, repr(reactant), repr(product)], "_") * ".dat"
end

function load_ionization_reactions(species)
    species_sorted = sort(species; by=x -> x.Z)
    reactions = IonizationReaction{LinearInterpolation{Float64,Float64}}[]
    for i in 1:length(species)
        for j in (i + 1):length(species)
            reaction = load_ionization_reaction(species_sorted[i], species_sorted[j])
            if !isnothing(reaction)
                push!(reactions, reaction)
            end
        end
    end

    return reactions
end