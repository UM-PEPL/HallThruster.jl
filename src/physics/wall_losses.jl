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
    1) `freq_electron_wall(model, params, config, i)`: Compute the electron-wall momentum transfer collision frequency in cell `i`
    2) `wall_power_loss!(Q, model, params, config)`: Compute the electron power lost to the walls in array Q

A third method, `wall_electron_current(model, params, i)`, will compute the electron current to the walls in cell `i`. If left unimplemented,
it defaults to Ie,w = e ne νew V_cell where V_cell is the cell volume.
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

function freq_electron_wall(::T, _params, _i) where {T <: WallLossModel}
    @nospecialize _params, _i
    error("freq_electron_wall not implemented for wall loss model of type $(nameof(T)). 
          See documentation for WallLossModel for a list of required methods")
end

function wall_power_loss(_Q, ::T, _params) where {T <: WallLossModel}
    @nospecialize _Q, _params
    error("wall_power_loss not implemented for wall loss model of type $(nameof(T)). 
          See documentation for WallLossModel for a list of required methods")
end

function wall_electron_current(::WallLossModel, params, i)
    (; grid, cache, grid, thruster) = params
    (; ne, νew_momentum) = cache
    A_ch = thruster.geometry.channel_area
    V_cell = A_ch * grid.dz_cell[i]
    return e * νew_momentum[i] * V_cell * ne[i]
end

#=============================================================================
# NoWallLosses
==============================================================================#

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

#=============================================================================
# ConstantSheathPotential
==============================================================================#

Base.@kwdef struct ConstantSheathPotential <: WallLossModel
    sheath_potential::Float64 = 20.0
    inner_loss_coeff::Float64
    outer_loss_coeff::Float64
end

function freq_electron_wall(::ConstantSheathPotential, _params, _i)
    @nospecialize(_params, _i)
    return 1.0e7
end

function wall_power_loss!(Q, model::ConstantSheathPotential, params)
    (; cache, grid, transition_length, thruster) = params
    (; ϵ) = cache
    (; sheath_potential, inner_loss_coeff, outer_loss_coeff) = model
    L_ch = thruster.geometry.channel_length

    @inbounds for i in 2:(length(grid.cell_centers) - 1)
        αϵ = linear_transition(
            grid.cell_centers[i], L_ch, transition_length,
            inner_loss_coeff, outer_loss_coeff,
        )
        Q[i] = 1.0e7 * αϵ * ϵ[i] * exp(-sheath_potential / ϵ[i])
    end

    return nothing
end

function wall_loss_scale(::ConstantSheathPotential)
    return 1.0
end

#=============================================================================
# WallMaterial
==============================================================================#

@kwdef struct WallMaterial
    name::String
    ϵ_star::Float64
    σ₀::Float64
end

@inline @fastmath function SEE_yield(material::WallMaterial, Tev, γ_max)
    (; σ₀, ϵ_star) = material
    γ = min(γ_max, σ₀ + 1.5 * Tev / ϵ_star * (1 - σ₀))
    return γ
end

const Alumina = WallMaterial(name = "Alumina", ϵ_star = 22, σ₀ = 0.57)
const BoronNitride = WallMaterial(name = "Boron Nitride", ϵ_star = 45, σ₀ = 0.24)
const SiliconDioxide = WallMaterial(name = "Silicon Dioxide", ϵ_star = 18, σ₀ = 0.5)
const BNSiO2 = WallMaterial(name = "Boron Nitride Silicon Dioxide", ϵ_star = 40, σ₀ = 0.54)
const SiliconCarbide = WallMaterial(name = "Silicon Carbide", ϵ_star = 43, σ₀ = 0.69)

#=============================================================================
 Serialization
==============================================================================#
const wall_materials = (; Alumina, BoronNitride, SiliconDioxide, BNSiO2, SiliconCarbide)
Serialization.SType(::Type{WallMaterial}) = Serialization.EnumType()
Serialization.options(::Type{WallMaterial}) = wall_materials

#=============================================================================
# WallSheath
==============================================================================#

function wall_electron_temperature(params, transition_length, i)
    (; cache, grid, thruster) = params

    shielded = thruster.shielded

    Tev = cache.Tev[i]

    Tev_channel = shielded * cache.Tev[1] + !shielded * Tev
    Tev_plume = Tev

    L_ch = thruster.geometry.channel_length

    Tev = linear_transition(
        grid.cell_centers[i], L_ch, transition_length, Tev_channel, Tev_plume,
    )

    return Tev
end

"""
    sheath_potential(Tev, γ, mi))
compute wall sheath to be used for radiative losses and loss to wall.
Goebel Katz equ. 7.3-29, 7.3-44. Assumed nₑuₑ/nᵢuᵢ ≈ 0.5
Sheath potentials are positive by convention in HallThruster.jl.
"""
@inline @fastmath sheath_potential(Tev, γ, mi) = Tev * log((1 - γ) * sqrt(mi / π / me / 2))

Base.@kwdef struct WallSheath <: WallLossModel
    material::WallMaterial
    loss_scale::Float64 = 1.0
    function WallSheath(material::WallMaterial, α::Float64 = 0.15)
        return new(material, α)
    end
end

# compute the edge-to-center density ratio
# (c.f https://iopscience.iop.org/article/10.1088/0963-0252/24/2/025017)
# this is approximate, but deviations only become noticable when the knudsen number becomes small, which is not true in our case
@inline edge_to_center_density_ratio() = 0.86 / sqrt(3);

# TODO: reorganize this into something that operates on arrays
function freq_electron_wall(model::WallSheath, params, i)
    (; cache, thruster, transition_length) = params
    geometry = thruster.geometry

    # compute difference in radii
    Δr = geometry.outer_radius - geometry.inner_radius
    # compute electron wall temperature
    Tev = wall_electron_temperature(params, transition_length, i)

    # calculate and store SEE coefficient
    # use number-averaged mass here
    γ_SEE_max = 1 - 8.3 * sqrt(me / cache.m_eff[i])
    γ = SEE_yield(model.material, Tev, γ_SEE_max)
    cache.γ_SEE[i] = γ

    # compute the ion current to the walls
    h = edge_to_center_density_ratio()
    j_iw = 0.0
    for fluid in params.fluid_containers.isothermal
        Z = fluid.species.Z
        inv_mi = inv(fluid.species.element.m)
        niw = h * fluid.density[i] * inv_mi
        j_iw += Z * model.loss_scale * niw * sqrt(Z * e * Tev * inv_mi)
    end

    # compute electron wall collision frequency
    νew = j_iw / (Δr * (1 - γ)) / cache.ne[i]

    return νew
end

function wall_power_loss!(Q, ::WallSheath, params)
    (; cache, grid, thruster, transition_length, plume_loss_scale) = params
    L_ch = thruster.geometry.channel_length

    @inbounds for i in 2:(length(grid.cell_centers) - 1)
        Tev = wall_electron_temperature(params, transition_length, i)

        # space charge limited SEE coefficient
        γ = params.cache.γ_SEE[i]

        # Space charge-limited sheath potential
        ϕ_s = sheath_potential(Tev, γ, cache.m_eff[i])

        # Compute electron wall collision frequency with transition function for energy wall collisions in plume
        νew = cache.radial_loss_frequency[i] * linear_transition(
            grid.cell_centers[i], L_ch, transition_length, 1.0, plume_loss_scale,
        )

        # Compute wall power loss rate
        Q[i] = νew * (2 * Tev + (1 - γ) * ϕ_s)
    end

    return nothing
end

@inline function wall_loss_scale(m::WallSheath)
    return m.loss_scale
end