Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = 20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

freq_electron_wall(::ConstantSheathPotential, params, i) = 1e7

function wall_power_loss(model::ConstantSheathPotential, params, i)
    (;z_cell, L_ch, config, cache) = params

    ne = cache.ne[i]
    nϵ = cache.nϵ[i]
    ϵ = nϵ / ne

    (;sheath_potential, inner_loss_coeff, outer_loss_coeff) = model
    αϵ = config.transition_function(z_cell[i], L_ch, inner_loss_coeff, outer_loss_coeff)

    W = 1e7 * αϵ * ϵ * exp(-sheath_potential / ϵ)

    return W
end
