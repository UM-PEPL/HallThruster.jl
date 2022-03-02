@inline function rk_step!(U, RHS!, params, rk_coeffs, cache, t, Δt)
    (;a, b, c) = rk_coeffs
    (;k, dU′, U′) = cache

    s = length(k)
    ℓ = length(U)

    RHS!(k[1], U, params, t)

    @inbounds for i in 2:s
        t′ = t + c[i-1]*Δt
        for j in 1:ℓ
            dU′[j] = 0.0
            for m in 1:i-1
                dU′[j] += a[i, m] * k[i][j]
            end
            U′[j] = U[j] + Δt * dU′[j]
        end
        RHS!(k[i], U′, params, t′)
    end

    @inbounds for j in 1:ℓ
        dU′[j] = 0.0
        for i in 1:s
            dU′[j] += b[i] * k[i][j]
        end
        U[j] += Δt * dU′[j]
    end

    return nothing
end

function simulate!(U, params, RHS!, rk_coeffs, Δt, tspan, saveat, ::SavedValues{T1, T2}, save_func, callback) where {T1, T2}
    s = size(rk_coeffs.a, 1)
    k = [zeros(size(U)) for i in 1:s]
    dU′ = zeros(size(U))
    U′ = zeros(size(U))
    cache = (;k, dU′, U′)

    nsave = length(saveat)
    saved_U = [zeros(size(U)) for i in 1:nsave]
    ts = zeros(length(saveat))
    savevals = Vector{T2}(undef, nsave)

    integrator = (cache = cache, p = params)

    t = tspan[1]
    iter = 0
    save_iter = 1
    while t < tspan[2]
        if t ≥ saveat[save_iter]
            @. saved_U[save_iter] = U
            ts[save_iter] = t
            savevals[save_iter] = save_func(U, t, integrator)
            save_iter += 1
        end
        iter += 1
        rk_step!(U, RHS!, params, rk_coeffs, cache, t, Δt)
        callback(U, params)
        if params.implicit_energy > 0
            energy_crank_nicholson!(U, params)
        end
        t += Δt
    end

    ts[end] = t
    @. saved_U[end] = U
    savevals[end] = save_func(U, t, integrator)

    return ts, saved_U, savevals
end