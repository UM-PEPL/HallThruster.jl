@inline check_r(r) = isfinite(r) && r >= 0

piecewise_constant(::Any) = 0.0
no_limiter(::Any) = 1.0
van_leer(r) = check_r(r) * (4r / (r + 1)^2)
van_albada(r) = check_r(r) * (2r / (r^2 + 1))
minmod(r) = check_r(r) * min(2 / (1 + r), 2r / (1 + r))
koren(r) = check_r(r) * max(0, min(2r, min((1 + 2r) / 3, 2))) * 2 / (r + 1)

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
