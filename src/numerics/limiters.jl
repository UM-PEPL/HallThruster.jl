"""
    $(TYPEDEF)
Defines a slope limiter for a second-order fluid solve.
See `HallThruster.slope_limiters` for a list of available limiters.
"""
struct SlopeLimiter{L}
    limiter::L
end

@inline check_r(r) = isfinite(r) && r >= 0

function (limiter::SlopeLimiter)(r)
    return check_r(r) * limiter.limiter(r)
end

__piecewise_constant(::Any) = 0.0
__no_limiter(::Any) = 1.0
__van_leer(r) = (4r / (r + 1)^2)
__van_albada(r) = (2r / (r^2 + 1))
__minmod(r) = min(2 / (1 + r), 2r / (1 + r))
__koren(r) = max(0, min(2r, min((1 + 2r) / 3, 2))) * 2 / (r + 1)

const piecewise_constant = SlopeLimiter(__piecewise_constant)
const no_limiter = SlopeLimiter(__no_limiter)
const van_leer = SlopeLimiter(__van_leer)
const van_albada = SlopeLimiter(__van_albada)
const minmod = SlopeLimiter(__minmod)
const koren = SlopeLimiter(__koren)

function stage_limiter!(U, params, config)
    (; grid, index, min_Te, cache) = params
    (; ncharge, propellant) = config
    stage_limiter!(U, grid.cell_centers, cache.nϵ, index, min_Te, ncharge, propellant.m)
end

function stage_limiter!(U, z_cell, nϵ, index, min_Te, ncharge, mi)
    min_density = MIN_NUMBER_DENSITY * mi
    @inbounds for i in eachindex(z_cell)
        U[index.ρn, i] = max(U[index.ρn, i], min_density)

        for Z in 1:ncharge
            density_floor = max(U[index.ρi[Z], i], min_density)
            velocity = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
            U[index.ρi[Z], i] = density_floor
            U[index.ρiui[Z], i] = density_floor * velocity
        end
        nϵ[i] = max(nϵ[i], 1.5 * MIN_NUMBER_DENSITY * min_Te)
    end
end

#=============================================================================
 Serialization
==============================================================================#
const slope_limiters = (;
    piecewise_constant, no_limiter, van_leer,
    van_albada, minmod, koren,
)
Serialization.SType(::Type{T}) where {T <: SlopeLimiter} = Serialization.Enum()
Serialization.options(::Type{T}) where {T <: SlopeLimiter} = slope_limiters
