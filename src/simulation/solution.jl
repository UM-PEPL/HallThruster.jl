struct Solution{T, U, P, S}
    t::T
    u::U
    savevals::S
    params::P
    retcode::Symbol
    error::String
end

function Solution(sol::S, params::P, savevals::SV, error::String = "") where {S, P, SV}
    return Solution(sol.t, sol.u, savevals, params, sol.retcode, error)
end

function Base.show(io::IO, ::MIME"text/plain", sol::Solution)
    println(io, "Hall thruster solution with $(length(sol.u)) saved frames")
    println(io, "Retcode: $(string(sol.retcode))")
    return print(io, "End time: $(sol.t[end]) seconds")
end

function Base.getindex(sol::Solution, field::Symbol, charge::Int = 1)
    if charge > sol.params.ncharge && field in [:ni, :ui, :niui]
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. \
             Maximum charge state in provided solution is $(sol.params.config.ncharge).",
        ))
    end

    if field == :ni
        return [saved[:ni][charge, :] for saved in sol.savevals]
    elseif field == :ui
        return [saved[:ui][charge, :] for saved in sol.savevals]
    elseif field == :niui
        return [saved[:niui][charge, :] for saved in sol.savevals]
    elseif field == :B
        return [sol.params.cache.B]
    elseif field == :cyclotron_frequency
        return [e * sol[:B, charge][1] / me]
    else
        return [getproperty(saved, field) for saved in sol.savevals]
    end
end

@inline _saved_fields_vector() = (
    :mobility,
    :Tev,
    :potential,
    :electric_field,
    :ne,
    :pe,
    :ue,
    :grad_pe,
    :nu_anom,
    :nu_class,
    :nu_en,
    :nu_ei,
    :radial_loss_frequency,
    :nu_wall,
    :nu_iz,
    :nu_ex,
    :nu_e,
    :Id,
    :ji,
    :nn,
    :anom_multiplier,
    :ohmic_heating,
    :wall_losses,
    :inelastic_losses,
    :Vs,
    :channel_area,
    :inner_radius,
    :outer_radius,
    :dA_dz,
    :tan_div_angle,
    :anom_variables,
    :dt
)

@inline _saved_fields_matrix() = (:ni, :ui, :niui)
@inline saved_fields() = (_saved_fields_vector()..., _saved_fields_matrix()...)

function solve(U, params, tspan; saveat, show_errors = true)
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

    try
        while t < tspan[2]
            # compute new timestep 
            if params.adaptive
                if uniform_steps
                    params.dt[] = params.dtbase
                    small_step_count -= 1
                else
                    params.dt[] = clamp(params.cache.dt[], params.dtmin, params.dtmax)
                end
            end

            t += params.dt[]

            #====
            Count how many times we've taken the minimum allowable timestep.
            If we exceed the threshold, then start taking longer uniform steps for a bit.
            This helps break out of cases where adaptive timestepping gets stuck ---
            either by resolving the situation or by causing the simulation to fail fast
            ====#
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

            # update heavy species quantities
            integrate_heavy_species!(U, params, params.dt[])
            update_heavy_species!(U, params)

            # Check for NaNs or Infs in heavy species solve and terminate if necessary
            if any(!isfinite, U)
                if show_errors
                    @warn("NaN or Inf detected in heavy species solver at time $(t)")
                end
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
        errstring = sprint(showerror, e, catch_backtrace())
        retcode = :error
        if show_errors
            @warn "Error detected in solution: $(errstring)"
        end
    end

    ind = min(save_ind, length(savevals) + 1) - 1

    return Solution(
        saveat[1:ind], u_save[1:ind], savevals[1:ind], params, retcode, errstring
    )
end
