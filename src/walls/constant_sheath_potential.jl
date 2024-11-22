Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = 20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

freq_electron_wall(::ConstantSheathPotential, params, i) = 1e7

function wall_power_loss!(Q, model::ConstantSheathPotential, params)
    (; config, cache, grid) = params
    (; ϵ) = cache
    (; sheath_potential, inner_loss_coeff, outer_loss_coeff) = model
    L_ch = config.thruster.geometry.channel_length

    @inbounds for i in 2:(length(grid.cell_centers) - 1)
        αϵ = linear_transition(grid.cell_centers[i], L_ch, config.transition_length,
            inner_loss_coeff, outer_loss_coeff,)
        Q[i] = 1e7 * αϵ * ϵ[i] * exp(-sheath_potential / ϵ[i])
    end

    return nothing
end
