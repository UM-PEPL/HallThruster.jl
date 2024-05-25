struct SlopeLimiter{F}
    limiter_func::F
end

function (L::SlopeLimiter)(r)
    valid = !(isinf(r) || isnan(r) || r <= 0.0)
    ϕ = valid * L.limiter_func(r)
    return ϕ
end

const piecewise_constant = SlopeLimiter(Returns(0.0))
const superbee = SlopeLimiter(r -> max(0.0, min(2r, 1), min(r, 2)))
const van_leer = SlopeLimiter(r -> 4r / (r + 1)^2)
const van_albada = SlopeLimiter(r -> 2r / (r^2 + 1))
const no_limiter = SlopeLimiter(Returns(1.0))
const minmod = SlopeLimiter(r -> min(2 / (1 + r), 2r / (1 + r)))
const koren = SlopeLimiter(r -> max(0, min(2r, min((1 + 2r) / 3, 2))) * 2 / (r+1))
const osher(β) = SlopeLimiter(r -> (max(0, min(r, β))) * 2 / (r+1))

function stage_limiter!(U, params)
    (;ncells, cache, config, index) = params
    (;nϵ) = cache
    min_density = config.min_number_density * config.propellant.m
    @inbounds for i in 1:ncells
        U[index.ρn, i] = max(U[index.ρn, i], min_density)

        for Z in 1:config.ncharge
            density_floor = max(U[index.ρi[Z], i], min_density)
            velocity = U[index.ρiui[Z], i] / U[index.ρi[Z], i]
            U[index.ρi[Z], i] = density_floor
            U[index.ρiui[Z], i] = density_floor * velocity
        end
        nϵ[i] = max(nϵ[i], config.min_number_density * config.min_electron_temperature)
    end
end
