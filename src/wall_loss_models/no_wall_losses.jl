struct NoWallLosses <: WallLossModel end

freq_electron_wall(model::NoWallLosses, U, params, i) = 0.0
wall_electron_current(model::NoWallLosses, U, params, i) = 0.0
wall_power_loss(model::NoWallLosses, U, params, i) = 0.0

