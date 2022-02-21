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
const koren = FluxLimiter(r -> max(0, min(2r, min((1 + 2r) / 3, r))))
const minmod = FluxLimiter(r -> max(0, min(1, r)))
const osher = FluxLimiter(r -> 1.5(r^2 + r) / (r^2 + r + 1))
const superbee = FluxLimiter(r -> max(0, min(2r, 1), min(r, 2)))
const van_albada = FluxLimiter(r -> (r^2 + r) / (r^2 + 1))
const van_albada_2 = FluxLimiter(r -> 2r / (1 + r^2))
const van_leer = FluxLimiter(r -> (r + abs(r)) / (1 + abs(r)))