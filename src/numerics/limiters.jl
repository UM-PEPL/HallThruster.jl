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

function stage_limiter!(U, p)
    ncells = size(U, 2)
    (;nϵ) = p.cache
    min_density = p.config.min_number_density * p.config.propellant.m
    @inbounds for i in 1:ncells
        U[p.index.ρn, i] = max(U[p.index.ρn, i], min_density)

        for Z in 1:p.config.ncharge
            density_floor = max(U[p.index.ρi[Z], i], min_density)
            velocity = U[p.index.ρiui[Z], i] / U[p.index.ρi[Z], i]
            U[p.index.ρi[Z], i] = density_floor
            U[p.index.ρiui[Z], i] = density_floor * velocity
        end
        nϵ[i] = max(nϵ[i], p.config.min_number_density * p.config.min_electron_temperature)
    end
end
