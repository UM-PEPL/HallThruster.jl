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

# Perform one step of the Strong-stability-preserving RK22 algorithm
function ssprk22_step!(U, f, params, dt)
    (;k, u1) = params.cache

    # First step of SSPRK22
    f(k, U, params)
    @. u1 = U + dt * k
    stage_limiter!(u1, params)

    # Second step of SSPRK22
    f(k, u1, params)
    @. U = (U + u1 + dt * k) / 2
    stage_limiter!(U, params)

    return nothing
end


function update_ions!(U, params)
    (;index, ncells, cache) = params
    (;nn, ne, ni, ui, niui, Z_eff, ji) = cache
    mi = params.config.propellant.m

    # Update heavy species
    ssprk22_step!(U, update_heavy_species!, params, params.dt[])

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    # Update plasma quantities
    @inbounds for i in 1:ncells
        # Compute number density for each neutral fluid
        nn[i] = U[index.ρn, i] / params.config.propellant.m

        # Compute ion derived quantities
        ne[i] = 0.0
        Z_eff[i] = 0.0
        ji[i] = 0.0
        @inbounds for Z in 1:params.config.ncharge
            _ni = U[index.ρi[Z], i] / mi
            _niui = U[index.ρiui[Z], i] / mi
            ni[Z, i] = _ni
            niui[Z, i] = _niui
            ui[Z, i] = _niui / _ni
            ne[i] += Z * _ni
            Z_eff[i] += _ni
            ji[i] += Z * e * _niui
        end

        # Compute electron number density, making sure it is above floor
        ne[i] = max(params.config.min_number_density, ne[i])

        # Effective ion charge state (density-weighted average charge state)
        Z_eff[i] = max(1.0, ne[i] / Z_eff[i])
    end
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

        # update heavy species quantities
        update_ions!(U, params)

        # Update electron quantities
        update_electrons!(params, t)

        # Update plume geometry
        update_plume_geometry!(U, params)

        # Allow for system interrupts
        yield()

        # Update the current iteration
        params.iteration[1] += 1

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
