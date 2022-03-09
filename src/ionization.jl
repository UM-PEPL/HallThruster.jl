Base.@kwdef struct IonizationReaction{I}
    reactant::Species
    product::Species
    rate_coeff::I
end

Base.@kwdef struct LandmarkTable{L, IL}
    loss_coeff::L
    rate_coeff::IL
end

struct LinearInterpolation{X<:Number,Y<:Number}
    xs::Vector{X}
    ys::Vector{Y}
    function LinearInterpolation(x, y)
        if length(x) != length(y)
            throw(ArgumentError("x and y must have same length"))
        else
            return new{typeof(x[1]),typeof(y[1])}(x, y)
        end
    end
end

function (itp::LinearInterpolation)(x::T) where {T}
    xs, ys = itp.xs, itp.ys
    if x ≤ xs[1]
        return ys[1] / oneunit(T)
    elseif x ≥ xs[end]
        return ys[end] / oneunit(T)
    end
    i = find_left_index(x, xs)
    return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])
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
        if value > array[mid + 1]
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
    return print(io, rxn_str)
end

function load_ionization_reaction(reactant, product)
    rates_file = rate_coeff_filename(reactant, product, "ionization")
    rates_file = joinpath(REACTION_FOLDER, rates_file)
    rates = DataFrame(CSV.File(rates_file))
    Te = rates[!, 1]
    k = rates[!, 2]
    rate_coeff = LinearInterpolation(Te, k)
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


@inline biexponential(x, c1, c2, c3, c4, c5) = c1 * (exp(c2 / (x - c5)) - c4 * exp(c2 * c3 / (x - c5)))


function ionization_fits_Xe(ncharge::Int)

    Xe0_Xe1 = IonizationReaction(
        Xenon(0), Xenon(1), ϵ -> biexponential(ϵ, 3.9e-13, 48.0, 0.0, 0.0, 6.0)
    )

    Xe0_Xe2 = IonizationReaction(
        Xenon(0), Xenon(2), ϵ -> biexponential(ϵ, 3.8e-14, 57.0, 10.0, 0.7, 0.0)
    )

    Xe0_Xe3 = IonizationReaction(
        Xenon(0), Xenon(3), ϵ -> biexponential(ϵ, 1.7e-14, 120, 6.0, 0.5, 0.0)
    )

    Xe1_Xe2 = IonizationReaction(
        Xenon(1), Xenon(2), ϵ -> biexponential(ϵ, 1.48e-13, 35.0, 11.0, 0.45, 0.0)
    )

    Xe1_Xe3 = IonizationReaction(
        Xenon(1), Xenon(3), ϵ -> biexponential(ϵ, 4e-14, 100, 4.0, 0.8, 0.0)
    )

    Xe2_Xe3 = IonizationReaction(
        Xenon(2), Xenon(3), ϵ -> biexponential(ϵ, 1.11e-13, 43.0, 7.0, 0.68, 0.0)
    )

    if ncharge == 1
        return [Xe0_Xe1]
    elseif ncharge == 2
        return [Xe0_Xe1, Xe0_Xe2, Xe1_Xe2]
    elseif ncharge == 3
        return [Xe0_Xe1, Xe0_Xe2, Xe0_Xe3, Xe1_Xe2, Xe1_Xe3, Xe2_Xe3]
    else
        throw(ArgumentError("ncharge must be 1, 2, or 3"))
    end
end

loss_coeff_fit(ϵ) = 6e-12 * exp(-39.0 / (ϵ + 3.0))