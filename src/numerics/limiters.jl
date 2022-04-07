struct FluxLimiter{F}
    limiter_func::F
end

function (L::FluxLimiter)(r)
    if isinf(r) || isnan(r) || r <= 0.0
        ϕ = 0.0
    else
        ϕ = L.limiter_func(r)
    end
    return ϕ
end

const no_limiter = FluxLimiter(identity)
const koren = FluxLimiter(r -> max(0, min(2r, min((1 + 2r) / 3, 2))))
const minmod = FluxLimiter(r -> max(0, min(1, r)))
const osher = FluxLimiter(r -> 1.5(r^2 + r) / (r^2 + r + 1))
const superbee = FluxLimiter(r -> max(0, min(2r, 1), min(r, 2)))
const van_albada = FluxLimiter(r -> (r^2 + r) / (r^2 + 1))
const van_leer = FluxLimiter(r -> (r + abs(r)) / (1 + abs(r)))

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