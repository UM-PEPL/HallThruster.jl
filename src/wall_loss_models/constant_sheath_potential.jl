Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = 20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

freq_electron_wall(::ConstantSheathPotential, params, i) = 1e7

function wall_power_loss!(Q, model::ConstantSheathPotential, params)
    (;z_cell, L_ch, config, cache, ncells) = params
    (;ϵ) = cache
    (;sheath_potential, inner_loss_coeff, outer_loss_coeff) = model

    @inbounds for i in 2:ncells-1
        αϵ = linear_transition(z_cell[i], L_ch, config.transition_length, inner_loss_coeff, outer_loss_coeff)
        Q[i] = 1e7 * αϵ * ϵ[i] * exp(-sheath_potential / ϵ[i])
    end

    return nothing
end
