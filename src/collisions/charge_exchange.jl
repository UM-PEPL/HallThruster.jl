Base.@kwdef struct ChargeExchangeReaction{I} <: Reaction
    reactant::Species
    product::Species
    rate_coeff::I
end

abstract type ChargeExchangeModel <: ReactionModel end

struct NoChargeExchange <: ReactionModel end

load_reactions(::NoChargeExchange, species) = ChargeExchangeReaction{Nothing}[]


struct ChargeExchangeFit <: ChargeExchangeModel end

supported_species(::ChargeExchangeFit) = [Xenon, Krypton]

"""
    Charge exchange fits from Hause, Prince and Bemish, 2013
"""
function load_reactions(::ChargeExchangeFit, species; kwargs...)

    CEX_fit(A, B, mi, ui) = let E = 0.5 * mi * ui^2 / e
        1e-20 * max(0.0, abs(ui) * (A - B * log(E)))
    end

    A_Xe_I,  B_Xe_I  = 87.3, 13.6
    A_Xe_II, B_Xe_II = 45.7, 8.9
    A_Kr_I,  B_Kr_I  = 80.7, 14.7
    A_Kr_II, B_Kr_II = 44.6, 9.8

    [
        let reactant = s
            product = if reactant.element == Xenon
                Xenon(0)
            elseif reactant.element == Krypton
                Krypton(0)
            else
                error("Unsupported gas $(s.element) provided for charge exchange reaction. This should be unreachable.")
            end

            A, B = if reactant == Xenon(1)
                A_Xe_I, B_Xe_I
            elseif reactant == Xenon(2)
                A_Xe_II, B_Xe_II
            elseif reactant == Krypton(1)
                A_Kr_I, B_Kr_I
            elseif reactant == Krypton(2)
                A_Kr_II, B_Kr_II
            else
                error("Unsupported charge state provided for charge exchange reaction. This should be unreachable.")
            end

            ChargeExchangeReaction(reactant, product, CEX_fit $ (A, B, reactant.element.m))
        end

        for s in species
        if s == Xenon(1) || s == Xenon(2) || s == Krypton(1) || s == Krypton(2)
    ]
end


