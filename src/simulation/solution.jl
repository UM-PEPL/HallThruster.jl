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
    :νei, :radial_loss_frequency, :νew_momentum, :νiz, :νex, :νe, :Id, :ji, :nn,
    :anom_multiplier, :ohmic_heating, :wall_losses, :inelastic_losses, :Vs,
    :channel_area, :inner_radius, :outer_radius, :dA_dz, :tanδ, :anom_variables,
    :dt
)

@inline _saved_fields_matrix() = (:ni, :ui, :niui)
@inline saved_fields() = (_saved_fields_vector()..., _saved_fields_matrix()...)

function solve(U, params, tspan; saveat)

    t = tspan[1]

    retcode = :success

    fields_to_save = saved_fields()

    first_saveval = NamedTuple{fields_to_save}(params.cache)
    u_save = [deepcopy(U) for _ in saveat]
    savevals = [deepcopy(first_saveval) for _ in saveat]

    (nvars, ncells) = size(U)
    iteration = params.iteration

    iteration[] = 1

    # Yield to signals only every few iterations
    yield_interval = 1
    next_yield = yield_interval
    save_ind = 2

    while t < tspan[2]

        if params.adaptive
            params.dt[] = clamp(params.cache.dt[], params.dtmin, params.dtmax)
        end

        t += params.dt[]

        # update heavy species quantities
        integrate_heavy_species!(U, params, params.dt[])
        update_heavy_species!(U, params)

        # Check for NaNs in heavy species solve and terminate if necessary
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

        if nandetected || infdetected
            break
        end

        # Update electron quantities
        update_electrons!(params, t)

        # Update plume geometry
        update_plume_geometry!(U, params)

        # Update the current iteration
        iteration[] += 1

        # Allow for system interrupts
        if iteration[] == next_yield
            yield()
            iteration[] += yield_interval
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
    end

    ind = min(save_ind, length(savevals)+1)-1

    return Solution(saveat[1:ind], u_save[1:ind], savevals[1:ind], retcode, params)
end
