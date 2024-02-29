# Perform one step of the Strong-stability-preserving RK22 algorithm
function ssprk22_step!(U, f, params, t, dt)
    (;k, u1) = params.cache

    # First step of SSPRK22
    f(k, U, params, t)
    @. u1 = U + dt * k
    stage_limiter!(u1, params)

    # Second step of SSPRK22
    f(k, u1, params, t+dt)
    @. U = (U + u1 + dt * k) / 2
    stage_limiter!(U, params)

    return nothing
end

struct Solution{T, U, P, S}
    t::T
    u::U
    savevals::S
    retcode::Symbol
    params::P
end

function Solution(sol::S, params::P, savevals::SV) where {S, P, SV}
    return Solution(sol.t, sol.u, savevals, sol.retcode, params)
end

function Base.show(io::IO, mime::MIME"text/plain", sol::Solution)
    println(io, "Hall thruster solution with $(length(sol.u)) saved frames")
    println(io, "Retcode: $(string(sol.retcode))")
    print(io, "End time: $(sol.t[end]) seconds")
end

@inline _saved_fields_vector() = (
    :μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen,
    :νei, :νew_energy, :νew_momentum, :νiz, :νex, :νe, :Id, :ji, :nn_tot,
    :anom_multiplier, :ohmic_heating, :wall_losses, :inelastic_losses, :Vs,
    :channel_area, :inner_radius, :outer_radius, :dA_dz, :tanδ, :anom_variables,
    :dt
)

@inline _saved_fields_matrix() = (:ni, :ui, :niui, :nn)
@inline saved_fields() = (_saved_fields_vector()..., _saved_fields_matrix()...)

function solve(U, params, tspan; saveat)
    i = 1
    save_ind = 2
    t = tspan[1]

    retcode = :success

    fields_to_save = saved_fields()

    first_saveval = NamedTuple{fields_to_save}(params.cache)
    u_save = [deepcopy(U) for _ in saveat]
    savevals = [deepcopy(first_saveval) for _ in saveat]

    (nvars, ncells) = size(U)

    while t < tspan[2]
        i += 1

        params.dt[] = params.adaptive ? params.cache.dt[] : params.dt[]
        params.dt[] = clamp(params.dt[], params.dtmin, params.dtmax)

        t += params.dt[]

        # Update heavy species
        ssprk22_step!(U, update_heavy_species!, params, t, params.dt[])

        # Update electron quantities
        update_electrons!(U, params, t)

        # Update plume geometry
        update_plume_geometry!(U, params)

        # Check for NaNs and terminate if necessary
        nandetected = false
        infdetected = false

        @inbounds for j in 1:ncells, i in 1:nvars
            if isnan(U[i, j])
                @warn("NaN detected in variable $i in cell $j at time $(t)")
                nandetected = true
                retcode = :NaNDetected
                break
            elseif isinf(U[i, j])
                @warn("Inf detected in variable $i in cell $j at time $(t)")
                infdetected = true
                retcode = :InfDetected
                break
            end
        end

        # Save values at designated intervals
        # TODO interpolate these to be exact and make a bit more elegant
        if t > saveat[save_ind]
            u_save[save_ind] .= U

            # save vector fields
            for field in _saved_fields_vector()
                if field == :anom_variables
                    for i in 1:num_anom_variables(params.config.anom_model)
                        savevals[save_ind][field][i] .= params.cache[field][i]
                    end
                else
                    cached_field::Vector{Float64} = params.cache[field]
                    sv::Vector{Float64} = savevals[save_ind][field]
                    sv .= cached_field
                end
            end

            # save matrix fields
            for field in _saved_fields_matrix()
                cached_field::Matrix{Float64} = params.cache[field]
                sv::Matrix{Float64} = savevals[save_ind][field]
                sv .= cached_field
            end

            save_ind += 1
        end

        if nandetected || infdetected
            break
        end
    end

    ind = min(save_ind, length(savevals)+1)-1

    return Solution(saveat[1:ind], u_save[1:ind], savevals[1:ind], retcode, params)
end
