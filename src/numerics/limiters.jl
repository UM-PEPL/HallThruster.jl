struct SlopeLimiter{F}
    limiter_func::F
end

function (L::SlopeLimiter)(r)
    if isinf(r) || isnan(r) || r <= 0.0
        ϕ = 0.0
    else
        ϕ = L.limiter_func(r)
    end
    return ϕ
end

const piecewise_constant_SL = SlopeLimiter(Returns(0.0))
const van_leer_SL = SlopeLimiter(r -> 4r / (r + 1)^2)
const van_albada_SL = SlopeLimiter(r -> 2r / (r^2 + 1))
const no_limiter_SL = SlopeLimiter(Returns(1.0))
const minmod_SL = SlopeLimiter(r -> min(2 / (1 + r), 2r / (1 + r)))
const koren_SL = SlopeLimiter(r -> max(0, min(2r, min((1 + 2r) / 3, 2))) * 2 / (r+1))
const osher_SL(β) = SlopeLimiter(r -> (max(0, min(r, β))) * 2 / (r+1))

function stage_limiter!(U, integrator, p, t)
    min_density = p.config.min_number_density * p.config.propellant.m
    @inbounds for j in 1:size(U, 2)
        U[p.index.ρn, j] = max(U[p.index.ρn, j], min_density)
        for Z in p.config.ncharge
            density_floor = max(U[p.index.ρi[Z], j], min_density)
            velocity = U[p.index.ρiui[Z], j] / U[p.index.ρi[Z], j]
            U[p.index.ρi[Z], j] = density_floor
            U[p.index.ρiui[Z], j] = density_floor * velocity
        end
        U[p.index.nϵ, j] = max(U[p.index.nϵ, j], p.config.min_number_density * p.config.min_electron_temperature)
    end
end