@public Solution, alternate_field_names, field_names

"""
$(TYPEDEF)

The solution of a simulation, returned by `run_simulation`.
These can be passed to any of the postprocessing functions described in [Postprocessing](@ref),
or indexed to extract specific values.

# Indexing
There are a few ways to index a solution.
First, you can extract a solution containing a single frame by indexing the `Solution` by an integer.

```jldoctest; setup = :(using HallThruster: HallThruster as het; config = het.Config(discharge_voltage = 300, thruster = het.SPT_100, anode_mass_flow_rate=5e-6, domain = (0.0, 0.08)); simparams = het.SimParams(dt = 5e-9, grid = het.EvenGrid(50), duration = 1e-3, num_save = 101, verbose=false); solution = het.run_simulation(config, simparams))
julia> solution = het.run_simulation(config, simparams)
Hall thruster solution with 101 saved frames (retcode: success, end time: 0.001 seconds)

julia> solution[51]
Hall thruster solution with 1 saved frame (retcode: success, end time: 0.0005 seconds)
```
Second, you can index a solution by a vector or vector-like object to extract a range of frames

```jldoctest; setup = :(using HallThruster: HallThruster as het; config = het.Config(discharge_voltage = 300, thruster = het.SPT_100, anode_mass_flow_rate=5e-6, domain = (0.0, 0.08)); simparams = het.SimParams(dt = 5e-9, grid = het.EvenGrid(50), duration = 1e-3, num_save = 101, verbose=false); solution = het.run_simulation(config, simparams))
julia> solution[51:end] # get last 51 frames
Hall thruster solution with 51 saved frames (retcode: success, end time: 0.001 seconds)

julia> solution[begin:2:end] # get every other frame
Hall thruster solution with 51 saved frames (retcode: success, end time: 0.001 seconds)

julia> solution[[1, 51, 101]] # get only frames 1, 51, 101
Hall thruster solution with 3 saved frames (retcode: success, end time: 0.001 seconds)
```
Lastly, you can index by a symbol or [Symbol, Integer] to get plasma data for all frames in that solution
```julia
solution[:ni, 1] 	# get ion density of first charge state for all frames
solution[:ui] 		# get ion velocity for all charge states and frames

# These return the same thing
solution[:∇pe]
solution[:grad_pe]
```

For a list of valid fields to index by, call `HallThruster.valid_fields()`
For a list of alternate names for fields containing special characters, call `HallThruster.alternate_field_names()`

See the documentation for `Base.getindex(sol::Solution, field::Symbol)` and `Base.getindex(sol::Solution, field::Symbol, charge::Integer)` for more information.

# Fields
$(TYPEDFIELDS)

"""
struct Solution{T, P, C, S}
	"""
	A vector of times (in seconds) at which simulation output has been saved
	"""
    t::T
	"""
	A vector of frames, or snapshots of the simulation state, at the times specified in `t`
	"""
    frames::S
	"""
	The solution parameters vector. Contains auxilliary information about the simulation.
	"""
    params::P
	"""
	The `Config` used to run the simulation
	"""
    config::C
	"""
	The solution return code. This can be one of three values:
	1. `:success`: the simulation completed successfully.
	2. `:failure`: the simulation failed due to a numerical issue or instability, resulting in a `NaN` or `Inf being detected somewhere in the solution`
	3. `:error`: another error occurred. Check the `error` string to see what kind of error.
	"""
    retcode::Symbol
	"""
	Holds to error text and backtrace, if an error occurred. Empty if `sol.retcode != :error`.
	"""
    error::String
end

function Base.show(io::IO, ::MIME"text/plain", sol::Solution)
	num_frames = length(sol.frames)
	retcode_str = string(sol.retcode)
	end_time = sol.t[end]
	plural = num_frames == 1 ? "" : "s"

	return print(io, "Hall thruster solution with $(num_frames) saved frame$(plural) \
			  (retcode: $(retcode_str), end time: $(end_time) seconds)")
end

Base.firstindex(sol::Solution) = 1
Base.lastindex(sol::Solution) = length(sol.t)

"""
$(TYPEDSIGNATURES)

Return a solution where `frames = [sol.frames[frame]]` and `[t = sol.t[frame]]`
All other fields remain unchanged.
"""
function Base.getindex(sol::Solution, frame::Integer)
    return Solution(
		[sol.t[frame]],
		[sol.frames[frame]],
        sol.params,
        sol.config,
        sol.retcode,
        sol.error,
    )
end

"""
$(TYPEDSIGNATURES)

Return a solution where `result.frames = sol.frames[frames] and result.t = sol.t[frames].`
This can be used to extract a contiguous slice of frames (by passing in a range like `50:end`) or a discrete sub-selection of frames (by passing in a vector like `[1, 51, 100]`)
"""
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

@inline _saved_fields_vector() = (:μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen,
    :νei, :radial_loss_frequency, :νew_momentum, :νiz, :νex,
    :νe, :Id, :ji, :nn,
    :anom_multiplier, :ohmic_heating, :wall_losses,
    :inelastic_losses, :Vs,
    :channel_area, :inner_radius, :outer_radius, :dA_dz,
    :tanδ, :anom_variables,
    :dt,)

@inline _saved_fields_matrix() = (:ni, :ui, :niui)

"""
$(SIGNATURES)
Returns a `Tuple` of symbols containing fields saved per frame.
These can be accessed by indexing a `Solution`.
Alternate names for those containing special characters can be found by calling `HallThruster.alternate_field_names()`

# Usage

```jldoctest; setup = :(using HallThruster)
julia> HallThruster.saved_fields()
$(saved_fields())
```
"""
@inline saved_fields() = (_saved_fields_vector()..., _saved_fields_matrix()...)

"""
$(SIGNATURES)
Returns a `NamedTuple` of mappings between alternate ascii field names and field names with special characters.
These can be used when indexing `Solution` objects instead of the short names.

# Usage

```jldoctest; setup = :(using HallThruster)
julia> HallThruster.alternate_field_names()
$(alternate_field_names())
```
"""
@inline alternate_field_names() = (;
	mobility = :μ,
	potential = :ϕ,
	thermal_conductivity = :κ,
	grad_pe = :∇pe,
	nu_anom = :νan,
	nu_class = :νc,
	nu_wall = :νew_momentum,	
	nu_ei = :νei,
	nu_en = :νen,
	nu_iz = :νiz,
	nu_ex = :νex,
	tan_divergence_angle = :tanδ,
)

"""
$(SIGNATURES)
Returns a `Tuple` of symbols containing fields that can be obtained by indexing a `Solution` by a `Symbol`.
This contains fields actually saved in a frame, in addition to special fields like `:z` and `:B`, as well as alternate field names for fields with special characters in their names (see `HallThruster.alternate_field_names()` for more)
```jldoctest; setup=:(using HallThruster)
julia> HallThruster.valid_fields()
$(valid_fields())
```
"""
function valid_fields()
	return (
		:z, :B, :E,
		saved_fields()...,
		keys(alternate_field_names())...,
		:E, :ωce, :cyclotron_freq,
	)
end

"""
$(TYPEDSIGNATURES)

Return plasma data indicated by the `field` for every frame in `sol`.
Type of returned data depends on the specific `field`.
A list of valid fiels can be found by calling `HallThruster.valid_fields()`.
Most of these return a vector of vectors, i.e. `[[field at time 0], [field at time 1], ...]`

For ion quantities, this method does not select a specific charge state.
Calling `sol[:ni]` returns a vector of `ncharge x ncells` matrices, each of which contains the density of ions on the grid for every charge state.
To get a specific charge, call `sol[:ni, Z]` where `1 <= Z <= ncharge` and `ncharge` is the maximum charge state of the simulation.

There are some special-cased convenience fields as well, which may return different values.
- `:B`: returns the magnetic field in each grid cell. Always returns a vector rather than vector of vectors, as the magnetic field is static.
- `:cyclotron_freq` or `ωce`: returns the electron cyclotron frequency (`e * B / m_e`)as a vector of vectors.
- `:E`: returns the electric field `-∇ϕ` as a vector of vectors.
- `:z`: returns the cell center locations for the grid as a vector.

Additionally, for values in `saved_fields` with non-ascii/special characters in their names, we provide alternate accessors, a list of which can be found by calling `HallThruster.alternate_field_names()`
"""
function Base.getindex(sol::Solution, field::Symbol)

	# Transform alternate field name, if needed
	alts = alternate_field_names()
	if field in keys(alts)
		field = alts[field]
	end

	if field in saved_fields()
		return [getproperty(frame, field) for frame in sol.frames]
	end

	# Special cases
    if field == :B
        return sol.params.cache.B
    elseif field == :ωce || field == :cyclotron_freq || field == :omega_ce
        return @. e * sol.params.cache.B / me
    elseif field == :E
        return -sol[:∇ϕ]
	elseif field == :z
		return sol.params.grid.cell_centers
	end

	throw(ArgumentError("Field :$(field) not found! Valid fields are $(valid_fields())"))
end

"""
$(TYPEDSIGNATURES)

For ion quantities (`:ni`, `:ui`, and `:niui`), indexing as `sol[field, charge]` returns a vector of vectors with the field for `charge`-charged ions.
As an example, `sol[:ui, 1]` returns the velocity of singly-charged ions for every frame in `sol.frames`.
For non-ion quantities, passing an integer as a second index causes an error.
"""
function Base.getindex(sol::Solution, field::Symbol, charge::Integer)
	is_ion_quantity = field in _saved_fields_matrix()

	if !is_ion_quantity
		throw(ArgumentError("Indexing a `solution` by `[field::Symbol, ::Integer]` is only supported for ion \
							quantities. To access a quantity at a specific frame, call `sol[field][frame]`."))
	end

    if charge <= 0 || charge > sol.config.ncharge 
		throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. \
							Maximum charge state in provided solution is $(sol.config.ncharge)."))
    end

	return [frame[field][charge, :] for frame in sol.frames]
end


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
