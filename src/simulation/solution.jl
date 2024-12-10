struct Solution{T, P, C, S}
    t::T
    frames::S
    params::P
    config::C
    retcode::Symbol
    error::String
end

function Base.show(io::IO, ::MIME"text/plain", sol::Solution)
	return print(io, "Hall thruster solution with $(length(sol.frames)) saved frames (retcode: $(string(sol.retcode)), end time: $(sol.t[end]) seconds)")
end

function Base.getindex(sol::Solution, frames::Integer)
    return Solution(
		[sol.t[frames]],
		[sol.frames[frames]],
        sol.params,
        sol.config,
        sol.retcode,
        sol.error,
    )
end

function Base.getindex(sol::Solution, frames::AbstractVector)
    return Solution(
		sol.t[frames],
		sol.frames[frames],
        sol.params,
        sol.config,
        sol.retcode,
        sol.error,
    )
end

function Base.getindex(sol::Solution, field::Symbol, charge::Integer = 1)
    if charge > sol.config.ncharge && field in [:ni, :ui, :niui]
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. Maximum charge state in provided solution is $(sol.config.ncharge)."))
    end

    if field == :ni
        return [saved[:ni][charge, :] for saved in sol.frames]
    elseif field == :ui
        return [saved[:ui][charge, :] for saved in sol.frames]
    elseif field == :niui
        return [saved[:niui][charge, :] for saved in sol.frames]
    elseif field == :B
        return sol.params.cache.B
    elseif field == :ωce || field == :cyclotron_freq || field == :omega_ce
        return e * sol[:B] / me
    elseif field == :E
        return -sol[:∇ϕ]
    elseif field == :V || field == :potential
        return sol[:ϕ]
    elseif field == :nu_an
        return sol[:νan]
    elseif field == :nu_ei
        return sol[:νei]
    elseif field == :nu_en
        return sol[:νen]
    elseif field == :mobility
        return sol[:μ]
    elseif field == :thermal_conductivity
        return sol[:κ]
	elseif field == :z || field == :cell_centers
		return sol.params.grid.cell_centers
    else
        return [getproperty(saved, field) for saved in sol.frames]
    end
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
    frames = [deepcopy(first_saveval) for _ in saveat]

    # Parameters for adaptive timestep escape hatch
    small_step_count = 0
    uniform_steps = false

    sim = params.simulation

    # Extract stuff from config
    (; source_neutrals, source_ion_continuity, source_ion_momentum, scheme) = config
    sources = (; source_neutrals, source_ion_continuity, source_ion_momentum)

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
            integrate_heavy_species!(U, params, scheme, sources, params.dt[])
            update_heavy_species!(U, params)

            # Check for NaNs or Infs in heavy species solve and terminate if necessary
            if any(!isfinite, U)
                if sim.print_errors
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
                # save vector fields
                for field in _saved_fields_vector()
                    if field == :anom_variables
                        for i in 1:num_anom_variables(config.anom_model)
                            frames[save_ind][field][i] .= params.cache[field][i]
                        end
                    else
                        cached_field::Vector{Float64} = params.cache[field]
                        sv::Vector{Float64} = frames[save_ind][field]
                        sv .= cached_field
                    end
                end

                # save matrix fields
                for field in _saved_fields_matrix()
                    cached_field::Matrix{Float64} = params.cache[field]
                    sv::Matrix{Float64} = frames[save_ind][field]
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

    ind = min(save_ind, length(frames) + 1) - 1

    return Solution(
        saveat[1:ind],
		frames[1:ind],
		params,
		config,
		retcode,
		errstring,
    )
end
