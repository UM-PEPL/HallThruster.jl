@public Solution, SpeciesState, alternate_field_names, field_names

"""
$(TYPEDEF)

The properties of a heavy species (neutrals or ions).

# Fields
$(TYPEDFIELDS)
"""
struct SpeciesState
    """Number density (1/m^3)"""
    n::Vector{Float64}
    """Number flux (1/m^2 s)"""
    nu::Vector{Float64}
    """Average velocity (m/s)"""
    u::Vector{Float64}
    """Molecular weight (kg)"""
    m::Float64
    """Charge number"""
    Z::Int8
    function SpeciesState(n::Int, m::Float64, Z::Integer)
        return new(zeros(n), zeros(n), zeros(n), m, Int8(Z))
    end
end

function _get_species_states(fluids_by_propellant)
    neutrals = OrderedDict{Symbol, SpeciesState}()
    ions = OrderedDict{Symbol, Vector{SpeciesState}}()

    for fluids in fluids_by_propellant
        continuity = fluids.continuity[1]
        m = continuity.species.element.m
        inv_m = 1 / m

        symbol = continuity.species.element.short_name
        neutral_state = SpeciesState(length(continuity.density), m, 0)
        @. neutral_state.n = continuity.density * inv_m
        @. neutral_state.u = continuity.const_velocity
        @. neutral_state.nu = neutral_state.n * neutral_state.u
        neutrals[symbol] = neutral_state

        ion_states = SpeciesState[]
        for ion in fluids.isothermal
            ion_state = SpeciesState(length(ion.density), m, ion.species.Z)
            @. ion_state.n = ion.density * inv_m
            @. ion_state.nu = ion.momentum * inv_m
            @. ion_state.u = ion.momentum / ion.density
            push!(ion_states, ion_state)
        end
        ions[symbol] = ion_states
    end

    return neutrals, ions
end

"""
$(TYPEDEF)

A snapshot of the simulation state at a single time, obtained by indexing the `frames` field of a `Solution` object.
Both neutral and ion species properties are stored as [`SpeciesState`](@ref) objects in a dictionary.
To access one of these objects, index by the symbol of that propellant and (if an ion species) the charge state.
For example, the number density of doubly-charged Xenon would be accessed as `frame.ions[:Xe][1].n`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct Frame
    """Dictionary containing neutral species. Indexed by that species' symbol."""
    neutrals::OrderedDict{Symbol, SpeciesState}
    """Dictionary containing ion species. Indexed by the species' symbol, followed by charge state."""
    ions::OrderedDict{Symbol, Vector{SpeciesState}}
    """Magnetic field strength (T)"""
    B::Vector{Float64}
    """Plasma density (1/m^3)"""
    ne::Vector{Float64}
    """Electron velocity (m/s)"""
    ue::Vector{Float64}
    """Ion current (A/m^2)"""
    ji::Vector{Float64}
    """Electric field (V/m)"""
    E::Vector{Float64}
    """Electron temperature (eV)"""
    Tev::Vector{Float64}
    """Electron pressure (eV/m^3)"""
    pe::Vector{Float64}
    """Electron pressure gradient (eV/m^4)"""
    grad_pe::Vector{Float64}
    """Electrostatic potential (V)"""
    potential::Vector{Float64}
    """Electron mobility"""
    mobility::Vector{Float64}
    """Anomalous collision frequency (1/s)"""
    nu_an::Vector{Float64}
    """Electron-neutral collision frequency (1/s)"""
    nu_en::Vector{Float64}
    """Electron-ion collision frequency (1/s)"""
    nu_ei::Vector{Float64}
    """Electron-wall collision frequency (1/s)"""
    nu_wall::Vector{Float64}
    """Total electron classical collision frequency (1/s)"""
    nu_class::Vector{Float64}
    """Total electron ionization collision frequency (1/s)"""
    nu_iz::Vector{Float64}
    """Total electron excitation collision frequency (1/s)"""
    nu_ex::Vector{Float64}
    """Total electron momentum transfer collision frequency (1/s)"""
    nu_e::Vector{Float64}
    """Cross-sectional area of discharge (m^2)"""
    channel_area::Vector{Float64}
    """Cross sectional area gradient (m^2 / m)"""
    dA_dz::Vector{Float64}
    """Tangent of plume divergence angle"""
    tan_div_angle::Vector{Float64}
    """Auxilliary anomalous transport variable caches"""
    anom_variables::Vector{Vector{Float64}}
    """Anomalous transport multiplier from PID controller"""
    anom_multiplier::Array{Float64, 0}
    """Discharge current (A)"""
    discharge_current::Array{Float64, 0}
    """SImulation timestep (s)"""
    dt::Array{Float64, 0}
end

function Frame(fluids_by_propellant, cache)
    neutrals, ions = _get_species_states(fluids_by_propellant)
    return Frame(;
        neutrals,
        ions,
        B = copy(cache.B),
        ne = copy(cache.ne),
        ue = copy(cache.ue),
        ji = copy(cache.ji),
        E = -copy(cache.∇ϕ),
        Tev = copy(cache.Tev),
        pe = copy(cache.pe),
        grad_pe = copy(cache.∇pe),
        potential = copy(cache.ϕ),
        mobility = copy(cache.μ),
        nu_an = copy(cache.νan),
        nu_en = copy(cache.νen),
        nu_ei = copy(cache.νei),
        nu_wall = copy(cache.νew_momentum),
        nu_class = copy(cache.νc),
        nu_iz = copy(cache.νiz),
        nu_ex = copy(cache.νex),
        nu_e = copy(cache.νe),
        channel_area = copy(cache.channel_area),
        dA_dz = copy(cache.dA_dz),
        tan_div_angle = copy(cache.tanδ),
        anom_variables = [copy(var) for var in cache.anom_variables],
        anom_multiplier = fill(cache.anom_multiplier[]),
        discharge_current = fill(cache.Id[]),
        dt = fill(cache.dt[]),
    )
end

"""
$(TYPEDEF)

The solution of a simulation, returned by `run_simulation`.
These can be passed to any of the postprocessing functions described in [Postprocessing](@ref),
or indexed to extract specific values.

# Indexing
There are a few ways to index a solution.
First, you can extract a solution containing a single `Frame` by indexing the `Solution` by an integer.

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

Lastly, you can obtain plasma properties by indexing into the frames contained in the solution

```jldoctest; setup = :(using HallThruster: HallThruster as het; config = het.Config(discharge_voltage = 300, thruster = het.SPT_100, anode_mass_flow_rate=5e-6, domain = (0.0, 0.08)); simparams = het.SimParams(dt = 5e-9, grid = het.EvenGrid(50), duration = 1e-3, num_save = 101, verbose=false); solution = het.run_simulation(config, simparams))
julia> solution.frames[4].Tev;                   # Electron temperature at fourth frame

julia> solution.frames[end].ions[:Xe][1].u;      # Ion velocity of Xenon+ at last frame

julia> solution.frames[end].neutrals[:Xe].n;     # Number density of neutral Xenon
```

For a list of valid fields in a `Frame`, call `fieldnames(HallThruster.Frame)`

# Fields
$(TYPEDFIELDS)

"""
struct Solution{T, C <: Config, CC <: CurrentController}
    """
    A vector of times (in seconds) at which simulation output has been saved
    """
    t::T
    """
    A vector of `Frame` objects representing snapshots of the simulation state, at the times specified in `t`
    """
    frames::Vector{Frame}
    """
    The grid used for the simulation
    """
    grid::Grid1D
    """
    The `Config` used to run the simulation
    """
    config::C
    """
    The simulation parameters
    """
    simulation::SimParams{CC}
    """
    The postprocessing arguments
    """
    postprocess::Postprocess
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
        sol.grid,
        sol.config,
        sol.simulation,
        sol.postprocess,
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
        sol.grid,
        sol.config,
        sol.simulation,
        sol.postprocess,
        sol.retcode,
        sol.error,
    )
end

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
@inline saved_fields() = fieldnames(Frame)


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
    μ = :mobility,
    ϕ = :potential,
    ∇pe = :grad_pe,
    νan = :nu_an,
    nu_anom = :nu_a,
    νc = :nu_class,
    νei = :nu_ei,
    νen = :nu_en,
    νiz = :nu_iz,
    νex = :nu_ex,
    tanδ = :tan_div_angle,
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
        :z,
        fieldnames(Frame)...,
        :E, :ωce, :cyclotron_freq, :ni, :ui, :niui, :nn,
        keys(alternate_field_names())...,
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

    if field in fieldnames(Frame)
        return [getproperty(frame, field) for frame in sol.frames]
    end

    # Special cases
    if field == :B
        return sol.frames[1].B
    elseif field == :ωce || field == :cyclotron_freq || field == :omega_ce
        return @. e * sol.frames[1].B / me
    elseif field == :z
        return sol.grid.cell_centers
    end

    # Heavy species properties (backwards compatibility)
    # Only allow if one ion species present
    if length(sol.config.propellants) > 1
        throw(ArgumentError("Cannot index by :nn, :ni, :niui, or :ui when more than one propellant species present. Instead, please index by the specific propellant species."))
    end
    symbol = sol.config.propellants[1].gas.short_name
    ncharge = sol.config.propellants[1].max_charge

    ncells = length(sol.grid.cell_centers)

    if field == :nn
        return [frame.neutrals[symbol].n for frame in sol.frames]
    elseif field == :ni
        return [[frame.ions[symbol][Z].n[i] for Z in 1:ncharge, i in 1:ncells] for frame in sol.frames]
    elseif field == :niui
        return [[frame.ions[symbol][Z].nu[i] for Z in 1:ncharge, i in 1:ncells] for frame in sol.frames]
    elseif field == :ui
        return [[frame.ions[symbol][Z].u[i] for Z in 1:ncharge, i in 1:ncells] for frame in sol.frames]
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
    is_ion_quantity = field in (:ni, :ui, :niui)

    if !is_ion_quantity
        throw(ArgumentError("Indexing a `solution` by `[field::Symbol, ::Integer]` is only supported for ion \
                            quantities. To access a quantity at a specific frame, call `sol.frames[frame].field`."))
    end

    if charge <= 0 || charge > sol.config.propellants[1].max_charge
        throw(ArgumentError("No ions of charge state $charge in Hall thruster solution. \
                            Maximum charge state in provided solution is $(sol.config.propellants[1].max_charge)."))
    end

    return [frame[field][charge, :] for frame in sol.frames]
end


function solve(params, config, tspan; saveat)
    # Initialize starting time and iterations
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
    frames = [Frame(params.fluids_by_propellant, params.cache)]

    # Parameters for adaptive timestep escape hatch
    small_step_count = 0
    uniform_steps = false

    sim = params.simulation

    # Extract stuff from config
    (; source_heavy_species) = config

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

            any_nan = integrate_heavy_species!(params.fluid_containers, params, source_heavy_species, params.dt[])

            # Check for NaNs or Infs in heavy species solve and terminate if necessary
            if any_nan
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

            # # Save values at designated intervals
            # # TODO interpolate these to be exact and make a bit more elegant
            if t > saveat[save_ind]
                #     # save vector fields
                #     for field in _saved_fields_vector()
                #         if field == :anom_variables
                #             for i in 1:num_anom_variables(config.anom_model)
                #                 frames[save_ind][field][i] .= params.cache[field][i]
                #             end
                #         else
                #             cached_field::Vector{Float64} = params.cache[field]
                #             sv::Vector{Float64} = frames[save_ind][field]
                #             sv .= cached_field
                #         end
                #     end

                #     # save matrix fields
                #     for field in _saved_fields_matrix()
                #         cached_field::Matrix{Float64} = params.cache[field]
                #         sv::Matrix{Float64} = frames[save_ind][field]
                #         sv .= cached_field
                #     end

                push!(frames, Frame(params.fluids_by_propellant, params.cache))

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
        params.grid,
        config,
        params.simulation,
        params.postprocess,
        retcode,
        errstring,
    )
end
