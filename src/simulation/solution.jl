struct Solution{T, U, P, C, S}
    t::T
    u::U
    savevals::S
    params::P
    config::C
    retcode::Symbol
    error::String
end

function Solution(
        sol::S, params::P, config::C, savevals::SV, error::String = "",) where {S, P, C, SV}
    return Solution(sol.t, sol.u, savevals, params, config, sol.retcode, error)
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
    :dt,)

@inline _saved_fields_matrix() = (:ni, :ui, :niui)
@inline saved_fields() = (_saved_fields_vector()..., _saved_fields_matrix()...)

function solve(U, params, config, tspan; saveat)
    # Initialie starting time and iterations
    iteration = params.iteration
    t = tspan[1]
    iteration[] = 1

    # Yield to signals only every few iterations
    yield_interval = 100

    # Error handling
    errstring = ""
    retcode = :success

    # Frame saving setup
    save_ind = 2
    fields_to_save = saved_fields()
    first_saveval = NamedTuple{fields_to_save}(params.cache)
    u_save = [deepcopy(U) for _ in saveat]
    savevals = [deepcopy(first_saveval) for _ in saveat]

    # Parameters for adaptive timestep escape hatch
    small_step_count = 0
    uniform_steps = false

    sim = params.simulation

    try
        while t < tspan[2]
            # compute new timestep 
            if sim.adaptive
                if uniform_steps
                    params.dt[] = sim.dt
                    small_step_count -= 1
                else
                    params.dt[] = clamp(params.cache.dt[], sim.min_dt, sim.max_dt)
                end
            end

            t += params.dt[]

            #====
            Count how many times we've taken the minimum allowable timestep.
            If we exceed the threshold, then start taking longer uniform steps for a bit.
            This helps break out of cases where adaptive timestepping gets stuck ---
            either by resolving the situation or by causing the simulation to fail fast
            ====#
            if params.dt[] == sim.min_dt
                small_step_count += 1
            elseif !uniform_steps
                small_step_count = 0
            end

            if small_step_count >= sim.max_small_steps
                uniform_steps = true
            elseif small_step_count == 0
                uniform_steps = false
            end

            # update heavy species quantities
            integrate_heavy_species!(U, params, config, params.dt[])
            update_heavy_species!(U, params, config)

            # Check for NaNs or Infs in heavy species solve and terminate if necessary
            if any(!isfinite, U)
                if sim.show_errors
                    @warn("NaN or Inf detected in heavy species solver at time $(t)")
                end
                retcode = :failure
                break
            end

            # Update electron quantities
            update_electrons!(params, config, t)

            # Update plume geometry
            if config.solve_plume
                update_plume_geometry!(params)
            end

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
                        for i in 1:num_anom_variables(config.anom_model)
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
        errstring = sprint(showerror, e, catch_backtrace())
        retcode = :error
        if sim.print_errors
            @warn "Error detected in solution: $(errstring)"
        end
    end

    ind = min(save_ind, length(savevals) + 1) - 1

    return Solution(
        saveat[1:ind], u_save[1:ind], savevals[1:ind], params, config, retcode, errstring,
    )
end
