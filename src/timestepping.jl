@inline function rk_step!(U, RHS!, params, rk_coeffs, cache, t, Δt)
    (;a, b, c) = rk_coeffs
    (;k, dU′, U′) = cache

    s = length(k)
    N = length(U)

    RHS!(k[1], U, params, t)

    @inbounds for i in 2:s
        t′ = t + c[i-1]*Δt
        for j in 1:N
            dU′[j] = 0.0
            for m in 1:i-1
                dU′[j] += a[i, m] * k[i][j]
            end
            U′[j] = U[j] + Δt * dU′[j]
        end
        RHS!(k[i], U′, params, t′)
    end

    @inbounds for j in 1:N
        dU′[j] = 0.0
        for i in 1:s
            dU′[j] += b[i] * k[i][j]
        end
        U[j] += Δt * dU′[j]
    end

    return nothing
end

# NOTE: currently assumes tridiagonal jacobian
function crank_nicholson_step!(x, RHS!, jacobian!, params, cache, t, Δt)
    s = params.num_newton_iterations

    (;k_impl, J, A, b, x₀, f₀, fᵢ) = cache

    N = length(x)
    k = k_impl

    # compute initial value of f(x₀)
    RHS!(f₀, x, params, t)

    x₀ .= x

    for i in 1:s
        jacobian!(J, fᵢ, x, params, t)
        RHS!(fᵢ, x, params, t + Δt)

        # A = I - Δt/2 J
        # b = -xᵢ + x₀ + Δt/2 * (fₙ + fₙ₊₁)
        A[1,1] = 1 - Δt/2 * J[1,1]
        A[1,2] = -Δt/2 * J[1, 2]
        b[1] = -x[1] + x₀[1] + Δt/2 * (f₀[1] + fᵢ[1])
        for j in 2:N-1
            A[j, j] = 1 - Δt/2 * J[j, j]
            A[j, j+1] = -Δt/2 * J[j, j+1]
            A[j, j-1] = -Δt/2 * J[j, j-1]
            b[j] = -x[j] + x₀[j] + Δt/2 * (f₀[j] + fᵢ[j])
        end
        A[end, end-1] = -Δt/2 * J[end-1, end-1]
        A[end, end] = 1 - Δt/2 * J[end, end]
        b[end] = -x[end] + x₀[end] + Δt/2 * (f₀[end] + fᵢ[end])

        # kᵢ = A \ b
        tridiagonal_solve!(k[i], A, b)

        # update xᵢ
        x .+= k[i]
    end

    return nothing
end

function simulate!(U, params, RHS!, RHS_impl!, rk_coeffs, Δt, tspan, saveat, ::SavedValues{T1, T2}, save_func, callback) where {T1, T2}
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
    implicit = RHS_impl! !== nothing
    index = params.index

    if implicit
        nϵ = U[index.nϵ, :]
        ncells = length(nϵ)
        dnϵ = zeros(ncells)
        cfg = ForwardDiff.JacobianConfig(RHS_impl!, dnϵ, nϵ)

        jac! = @closure (J, dU, U, params, t) -> let
            f = @closure (dU, U) -> RHS_impl!(dU, U, params, t)
            ForwardDiff.jacobian!(J, f, dU, U)
        end
        A = Tridiagonal(zeros(ncells-1), zeros(ncells), zeros(ncells-1))
        J = copy(A)
        b = zeros(ncells)
        k_impl = [zeros(ncells) for i in 1:params.num_newton_iterations]

        cache_impl = (;k_impl, A, b, J, x₀ = zeros(ncells), f₀ = dnϵ, fᵢ = zeros(ncells))
    end

    while t < tspan[2]
        if t ≥ saveat[save_iter]
            @. saved_U[save_iter] = U
            ts[save_iter] = t
            savevals[save_iter] = save_func(U, t, integrator)
            save_iter += 1
        end
        iter += 1
        rk_step!(U, RHS!, params, rk_coeffs, cache, t, Δt)
        if implicit
            @. @views nϵ = U[index.nϵ, :]
            crank_nicholson_step!(nϵ, RHS_impl!, jac!, params, cache_impl, t, Δt)
        end
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