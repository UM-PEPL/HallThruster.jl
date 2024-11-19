"""
    $(TYPEDFIELDS)
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

const slope_limiters = (;
    piecewise_constant,
    no_limiter,
    van_leer,
    van_albada,
    minmod,
    koren
)

@__register_stringtype(SlopeLimiter, slope_limiters)

function stage_limiter!(U, params)
    (; ncells, cache, config, index) = params
    (; nϵ) = cache
    min_density = config.min_number_density * config.propellant.m
    @inbounds for i in 1:ncells
        U[index.ρn, i] = max(U[index.ρn, i], min_density)

        for Z in 1:(config.ncharge)
            density_floor = max(U[index.ρi[Z], i], min_density)
            velocity = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
            U[index.ρi[Z], i] = density_floor
            U[index.ρiui[Z], i] = density_floor * velocity
        end
        nϵ[i] = max(nϵ[i], config.min_number_density * config.min_electron_temperature)
    end
end
