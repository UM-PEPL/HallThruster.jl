"""
$(TYPEDEF)

Struct containing secondary electron emission (SEE) yield fit coefficients for a material. Used in the `WallSheath` wall loss model

# Fields

$(TYPEDFIELDS)

The SEE yield is computed as a function of electron temperature in eV (`Tev`) as `Γ * a * Tev^b`.

# Available materials

- `IdealDielectric` (a material with zero SEE)
- `Alumina`
- `BoronNitride`
- `BNSiO2`
- `StainlessSteel`

"""
Base.@kwdef struct WallMaterial
    name::String = ""
    a::Float64
    b::Float64
    Γ::Float64
end

"""
    SEE_yield(material::WallMaterial, Tev, mi)
fit function for SEE with different wall materials
Goebel Katz equ. 7.3-30
"""
@inline function SEE_yield(material::WallMaterial, Tev, mi)
    (;a, b, Γ) = material
    γ_max = 1 - 8.3 * sqrt(me/mi) # Space-charge limited SEE coefficient
    γ = min(γ_max,  Γ * a * Tev^b)
    return γ
end

const IdealDielectric = WallMaterial(name = "IdealDielectric", a = 0.0,   b = 0.0,   Γ = 0.0)
const Alumina         = WallMaterial(name = "Alumina", a = 0.145, b = 0.650, Γ = 1.49)
const BoronNitride    = WallMaterial(name = "BoronNitride", a = 0.150, b = 0.549, Γ = 1.38)
const BNSiO2          = WallMaterial(name = "BNSiO2", a = 0.123, b = 0.528, Γ = 1.36)
const StainlessSteel  = WallMaterial(name = "StainlessSteel", a = 0.040, b = 0.610, Γ = 1.44)

Base.@kwdef struct WallSheath <: WallLossModel
    material::WallMaterial
    α::Float64 = 0.15 # scaling coefficient
    function WallSheath(material::WallMaterial, α::Float64 = 0.15)
        return new(material, α)
    end
end



function wall_power_loss(model::WallSheath, U, params, i)
    (;config) = params
    mi = config.propellant.m
    Tev = wall_electron_temperature(U, params, i)

    # space charge limited SEE coefficient
    γ = SEE_yield(model.material, Tev, mi)

    # Space charge-limited sheath potential
    ϕ_s = sheath_potential(Tev, γ, mi)

    νew = freq_electron_wall(model, U, params, i)

    W = νew * (2Tev + (1 - γ) * ϕ_s)

    return W
end

function wall_electron_current(model::WallSheath, U, params, i)
    (;config) = params
    mi = config.propellant.m
    Tev = wall_electron_temperature(U, params, i)
    γ = SEE_yield(model.material, Tev, mi)

    return inv(1 - γ) * sum(wall_ion_current(model, U, params, i, Z) for Z in 1:config.ncharge)
end

function wall_ion_current(model::WallSheath, U, params, i, Z)
    (;z_edge, config, index, z_cell, L_ch) = params
    (;propellant, thruster, transition_function) = config
    (;α) = model

    mi = propellant.m
    ni = U[index.ρi[Z], i] / mi
    Tev = params.cache.Tev[i]

    u_bohm = sqrt(Z * e * Tev / mi)

    if firstindex(z_cell) < i < lastindex(z_cell)
        Δz = z_edge[right_edge(i)] - z_edge[left_edge(i)]
    elseif i == firstindex(z_cell)
        Δz = z_edge[begin+1] - z_edge[begin]
    elseif i == lastindex(z_cell)
        Δz = z_edge[end] - z_edge[end-1]
    end

    in_channel = transition_function(z_cell[i], L_ch, 1.0, 0.0)

    Iiw = in_channel * α * Z * e * ni * u_bohm * channel_perimeter(thruster) * Δz

    return Iiw
end

function freq_electron_wall(model::WallSheath, U, params, i)
    (;A_ch, z_edge, z_cell) = params

    if firstindex(z_cell) < i < lastindex(z_cell)
        Δz = z_edge[right_edge(i)] - z_edge[left_edge(i)]
    elseif i == firstindex(z_cell)
        Δz = z_edge[begin+1] - z_edge[begin]
    elseif i == lastindex(z_cell)
        Δz = z_edge[end] - z_edge[end-1]
    end

    V_cell = A_ch * Δz

    Iew = wall_electron_current(model, U, params, i)

    ne = params.cache.ne[i]

    νew = Iew / e / ne / V_cell

    return νew
end

function wall_electron_temperature(U, params, i)
    (;cache, config, z_cell) = params

    Tev_channel = if config.thruster.shielded
        cache.Tev[1]
    else
        cache.Tev[i]
    end

    Tev_plume = cache.Tev[i]

    L_ch = config.thruster.geometry.channel_length

    Tev = config.transition_function(z_cell[i], L_ch, Tev_channel, Tev_plume)

    return Tev
end

"""
    sheath_potential(Tev, γ, mi))
compute wall sheath to be used for radiative losses and loss to wall.
Goebel Katz equ. 7.3-29, 7.3-44. Assumed nₑuₑ/nᵢuᵢ ≈ 0.5
Sheath potentials are positive by convention in HallThruster.jl.
"""
sheath_potential(Tev, γ, mi) = Tev*log((1 - γ) * sqrt(mi/π/me/2))
