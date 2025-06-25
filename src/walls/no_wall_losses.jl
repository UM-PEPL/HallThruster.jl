struct NoWallLosses <: WallLossModel end

function freq_electron_wall(::NoWallLosses, _params, _i)
    @nospecialize _params, _i
    return 0.0
end

function wall_electron_current(::NoWallLosses, _params, _i)
    @nospecialize _params, _i
    return 0.0
end

function wall_power_loss!(Q, ::NoWallLosses, _params)
    @nospecialize _params
    return Q .= 0.0
end

wall_loss_scale(::NoWallLosses) = 0.0
