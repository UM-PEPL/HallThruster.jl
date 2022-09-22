# Perform one step of the Strong-stability-preserving RK22 algorithm
function ssprk22_step!(u, f, params, t)
    (;dt, cache) = params
    (;k, u1) = cache

    # First step of SSPRK22
    f(k, u, params, t)
    @. u1 = u + dt * k
    stage_limiter!(u1, nothing, params, t)

    # Second step of SSPRK22
    f(k, u1, params, t+dt)
    @. u = (u + u1 + dt * k) / 2
    stage_limiter!(u, nothing, params, t)

    return nothing
end

struct MyODEProblem{F, U, T, P}
    f::F
    u::U
    tspan::T
    params::P
end

function mysolve(prob::MyODEProblem; saveat, dt)
    (;f, u, tspan, params) = prob
    i = 1
    save_ind = 2
    t = tspan[1]

    retcode = :success

    fields_to_save = (
        :μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen,
        :νei, :νew, :νiz, :νex, :νe, :Id, :ni, :ui, :ji, :niui, :nn, :nn_tot
    )

    first_saveval = NamedTuple{fields_to_save}(params.cache)
    u_save = [deepcopy(u) for _ in saveat]
    savevals = [first_saveval for _ in saveat]

    (nvars, ncells) = size(u)

    while t < tspan[2]
        i += 1
        t += dt
        ssprk22_step!(u, f, params, t)
        update_values!(u, params, t)

        # Check for NaNs and terminate if necessary
        nandetected = false
        infdetected = false

        @inbounds for j in 1:ncells, i in 1:nvars
            if isnan(u[i, j])
                println("NaN detected in variable $i in cell $j at time $(integrator.t)")
                nandetected = true
                retcode = :NaNDetected
                break
            elseif isinf(u[i, j])
                println("Inf detected in variable $i in cell $j at time $(integrator.t)")
                infdetected = true
                retcode = :InfDetected
                break
            end
        end

        # Save values
        if t > saveat[save_ind]
            u_save[save_ind] .= u
            savevals[save_ind] = NamedTuple{fields_to_save}(params.cache)
            save_ind += 1
        end

        if nandetected || infdetected
            break
        end
    end

    return Solution(saveat, u_save, savevals, retcode, params)
end
