struct NoWallLosses <: WallLossModel end

freq_electron_wall(model::NoWallLosses, params, i) = 0.0
wall_electron_current(model::NoWallLosses, params, i) = 0.0
function wall_power_loss(Q, model::NoWallLosses, params)
  Q .= 0.0
end
