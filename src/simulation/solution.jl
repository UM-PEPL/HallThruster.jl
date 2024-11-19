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

function Base.show(io::IO, ::MIME"text/plain", sol::Solution)
    println(io, "Hall thruster solution with $(length(sol.u)) saved frames")
    println(io, "Retcode: $(string(sol.retcode))")
    return print(io, "End time: $(sol.t[end]) seconds")
end

@inline _saved_fields_vector() = (:μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen,
    :νei, :radial_loss_frequency, :νew_momentum, :νiz, :νex,
    :νe, :Id, :ji, :nn,
    :anom_multiplier, :ohmic_heating, :wall_losses,
    :inelastic_losses, :Vs,
    :channel_area, :inner_radius, :outer_radius, :dA_dz,
    :tanδ, :anom_variables,
    :dt)

@inline _saved_fields_matrix() = (:ni, :ui, :niui)
@inline saved_fields() = (_saved_fields_vector()..., _saved_fields_matrix()...)

function solve(U, params, tspan; saveat)
    fields_to_save = saved_fields()

    first_saveval = NamedTuple{fields_to_save}(params.cache)
    u_save = [deepcopy(U) for _ in saveat]
    savevals = [deepcopy(first_saveval) for _ in saveat]

    iteration = params.iteration
    t = tspan[1]
    iteration[] = 1

    # Yield to signals only every few iterations
    yield_interval = 100
    save_ind = 2
    small_step_count = 0
    uniform_steps = false
    retcode = :success

    try
        while t < tspan[2]
            # compute new timestep 
            if params.adaptive
                if uniform_steps
                    params.dt[] = params.dtbase
                    #println(small_step_count)
                    small_step_count -= 1
                else
                    params.dt[] = clamp(params.cache.dt[], params.dtmin, params.dtmax)
                end
            end

            # Count number of minimal timesteps in a row
            if params.dt[] == params.dtmin
                small_step_count += 1
            elseif !uniform_steps
                small_step_count = 0
            end

            if small_step_count >= params.max_small_steps
                uniform_steps = true
            elseif small_step_count == 0
                uniform_steps = false
            end

            t += params.dt[]

            # update heavy species quantities
            integrate_heavy_species!(U, params, params.dt[])
            update_heavy_species!(U, params)

            # Check for NaNs in heavy species solve and terminate if necessary
            if any(!isfinite, U)
                @warn "Nan or Inf detected in heavy species solve at time $(t)"
                retcode = :failure
                break
            end

            # Update electron quantities
            update_electrons!(params, t)

            # Update plume geometry
            update_plume_geometry!(U, params)

            # Update the current iteration
            iteration[] += 1

            # Allow for system interrupts
            if iteration[] % yield_interval == 0
                yield()
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
    catch e
        # catch and print error, but return a normal retcode
        errstring = sprint(showerror, e, catch_backtrace())
        println(errstring)
        retcode = :error
    end

    ind = min(save_ind, length(savevals) + 1) - 1

    return Solution(saveat[1:ind], u_save[1:ind], savevals[1:ind], retcode, params)
end
