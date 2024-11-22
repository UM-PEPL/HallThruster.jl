"""
$(TYPEDEF)

Abstract type for wall loss models in the electron energy equation. Types included with HallThruster.jl are:

- `NoWallLosses`
    No electron energy losses to the wall.
- `ConstantSheathPotential(sheath_potential, inner_loss_coeff, outer_loss_coeff)`
    Power to the walls scales with `α * exp(-sheath_potential))`, where `α = inner_loss_coeff` inside the channel and `α = outer_loss_coeff` outside.
- `WallSheath(material::WallMaterial)`
    Power to the walls is computed self-consistently using a sheath model, dependent on the secondary electron emission yield of the provided material.
    See `WallMaterial` for material options.

Users implementing their own `WallLossModel` will need to implement at least three methods
    1) `freq_electron_wall(model, params, i)`: Compute the electron-wall momentum transfer collision frequency in cell `i`
    2) `wall_power_loss!(Q, model, params)`: Compute the electron power lost to the walls in array Q

A third method, `wall_electron_current(model, params, i)`, will compute the electron current to the walls in cell `i`. If left unimplemented,
it defaults to Ie,w = e ne νew V_cell where V_cell is the cell volume.

A fourth method, `wall_ion_current(model, params, i, Z)`, for computing the current of ions of charge Z to the walls in cell `i`, may also be implemented.
If left unimplemented, it will default to computing the current assuming Ie,w = Ii,w.
"""
abstract type WallLossModel end

#=============================================================================
 Serialization
==============================================================================#

wall_loss_models() = (; NoWallLosses, ConstantSheathPotential, WallSheath)
Serialization.SType(::Type{T}) where {T <: WallLossModel} = Serialization.TaggedUnion()
Serialization.options(::Type{T}) where {T <: WallLossModel} = wall_loss_models()

#=============================================================================
 Placeholder definitions
==============================================================================#

function freq_electron_wall(model::WallLossModel, ::Any, ::Any)
    error("freq_electron_wall not implemented for wall loss model of type $(typeof(model)). See documentation for WallLossModel for a list of required methods")
end

function wall_power_loss(::Any, model::WallLossModel, ::Any)
    error("wall_power_loss not implemented for wall loss model of type $(typeof(model)). See documentation for WallLossModel for a list of required methods")
end

function wall_electron_current(::WallLossModel, params, i)
    (; grid, config, cache) = params
    (; ne, nu_wall) = cache
    A_ch = config.thruster.geometry.channel_area
    V_cell = A_ch * grid.dz_cell[i]
    return e * nu_wall[i] * V_cell * ne[i]
end

function wall_ion_current(model::WallLossModel, params, i, Z)
    (; ne, ni) = params.cache

    Iew = wall_electron_current(model, params, i)
    return Z * ni[Z, i] / ne[i] * Iew * (1 - params.cache.γ_SEE[i])
end
