Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = 20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

function freq_electron_wall(::ConstantSheathPotential, U, params, i)
    (;z_cell, L_ch) = params
    αϵ = params.config.transition_function(z_cell[i], L_ch, 1.0, 0.0)
    return αϵ * 1e7
end

function wall_electron_current(model::ConstantSheathPotential, U, params, i)
    (;z_edge, cache, A_ch) = params
    (;ne) = cache
    νew = freq_electron_wall(model, U, params, i)
    Δz = z_edge[right_edge(i)] - z_edge[left_edge(i)]
    V_cell = A_ch * Δz
    return e * νew * V_cell * ne[i]
end

function wall_power_loss(model::ConstantSheathPotential, U, params, i)
    (;z_cell, index, L_ch, config, cache) = params

    ne = cache.ne[i]
    ϵ = U[index.nϵ, i] / ne

    (;sheath_potential, inner_loss_coeff, outer_loss_coeff) = model
    αϵ = config.transition_function(z_cell[i], L_ch, inner_loss_coeff, outer_loss_coeff)

    W = 1e7 * αϵ * ϵ * exp(-sheath_potential / ϵ)

    return W
end