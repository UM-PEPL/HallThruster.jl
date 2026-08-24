@public write_to_json, run_simulation, SERIALIZATION_VERSION

"""
Current version of the top-level HallThruster JSON serialization schema.
Unversioned documents are treated as legacy version 0.
"""
const SERIALIZATION_VERSION = 1

function _validate_serialization_version(document, source = "serialized document")
    document isa AbstractDict || throw(
        ArgumentError(
            "$(source) must contain a JSON object at its top level."
        )
    )

    version = get(document, "serialization_version", 0)
    version isa Integer || throw(
        ArgumentError(
            "$(source) has a non-integer `serialization_version`: $(repr(version))."
        )
    )
    version in 0:SERIALIZATION_VERSION || throw(
        ArgumentError(
            "$(source) uses unsupported serialization version $(version); " *
                "this HallThruster release supports versions 0 through $(SERIALIZATION_VERSION)."
        )
    )
    return Int(version)
end

"""
    $(TYPEDSIGNATURES)
Run a simulation from a JSON input file.
If `postprocess` is set in the JSON file and `postprocess.output_file` is non-empty, output will be written
to `postprocess.output_file`.
If `restart` is a JSON file, this function will also try to restart the simulation from that file.
Returns a `Solution` object.
"""
function run_simulation(json_file::String; restart::String = "")
    json_file = if !ispath(json_file)
        joinpath(@__DIR__, json_file)
    else
        json_file
    end

    if (splitext(json_file)[2] != ".json")
        throw(ArgumentError("$json_file is not a valid JSON file"))
    end

    obj = JSON.parsefile(json_file)
    _validate_serialization_version(obj, "JSON file $(json_file)")

    # Read config and sim params from file
    input = get(obj, "input", obj)
    cfg = deserialize(Config, input["config"])
    sim = deserialize(SimParams, input["simulation"])

    postprocess::Union{Postprocess, Nothing} = nothing
    if haskey(input, "postprocess") && haskey(input["postprocess"], "output_file") &&
            !isempty(input["postprocess"]["output_file"])
        postprocess = deserialize(Postprocess, input["postprocess"])
    end

    sol = run_simulation(cfg, sim; postprocess, include_dirs = dirname(json_file), restart)

    if postprocess !== nothing
        (; average_start_time, save_time_resolved) = postprocess
        write_to_json(postprocess.output_file, sol; average_start_time, save_time_resolved)
    end

    return sol
end

"""
    $(TYPEDSIGNATURES)
Convert one frame of a `Solution` to an `OrderedDict`
"""
Base.@nospecializeinfer function frame_dict(@nospecialize(sol::Solution), frame::Integer)
    f = sol.frames[frame]
    d = OrderedDict{String, Any}()
    d["thrust"] = thrust(sol, frame)
    d["discharge_current"] = discharge_current(sol, frame)
    d["ion_current"] = ion_current(sol, frame)
    d["mass_eff"] = mass_eff(sol, frame)
    d["voltage_eff"] = voltage_eff(sol, frame)
    d["current_eff"] = current_eff(sol, frame)
    d["divergence_eff"] = divergence_eff(sol, frame)
    d["anode_eff"] = anode_eff(sol, frame)
    d["t"] = sol.t[frame]
    d["z"] = sol.grid
    d["B"] = f.B
    d["ne"] = f.ne
    d["ue"] = f.ue
    d["potential"] = f.potential
    d["E"] = f.E
    d["Tev"] = f.Tev
    d["pe"] = f.pe
    d["grad_pe"] = f.grad_pe
    d["nu_en"] = f.nu_en
    d["nu_ei"] = f.nu_ei
    d["nu_anom"] = f.nu_an
    d["nu_class"] = f.nu_class
    d["mobility"] = f.mobility
    d["channel_area"] = f.channel_area

    d["neutrals"] = OrderedDict(
        symbol => OrderedDict(
                "n" => neutral.n,
                "u" => neutral.u,
                "nu" => neutral.nu,
            ) for (symbol, neutral) in pairs(f.neutrals)
    )

    d["ions"] = OrderedDict(
        symbol => [
                OrderedDict(
                    "n" => ion.n,
                    "u" => ion.u,
                    "nu" => ion.nu,
                    "Z" => ion.Z,
                )
                for ion in ions
            ]
            for (symbol, ions) in pairs(f.ions)
    )

    d["excited_states"] = OrderedDict(
        symbol => OrderedDict(
                "n" => state.n,
                "u" => state.u,
                "nu" => state.nu,
                "m" => state.m,
                "Z" => state.Z,
                "excited_level" => state.excited_level,
                "energy_eV" => state.energy_eV,
            ) for (symbol, state) in pairs(f.excited_states)
    )

    d["photon_emissions"] = [
        OrderedDict(
                "upper" => emission.upper,
                "lower" => emission.lower,
                "frequency" => emission.frequency,
                "energy_eV" => emission.energy_eV,
                "emission_rate" => emission.emission_rate,
            ) for emission in f.photon_emissions
    ]

    return d
end

"""
    $(TYPEDSIGNATURES)
Convert `sol` to an `OrderedDict`, containing both the inputs used to run the simulation
and any requested outputs.
This function is used to convert a `Solution` to a format suitable for writing to an output file.
"""
Base.@nospecializeinfer function serialize_sol(
        @nospecialize(sol::Solution);
        average_start_time::AbstractFloat = -1, save_time_resolved::Bool = true,
    )
    output = OrderedDict{String, Any}()
    output["retcode"] = string(sol.retcode)
    output["error"] = sol.error

    if average_start_time >= 0
        first_frame = findfirst(>=(average_start_time), sol.t)
        if first_frame === nothing
            first_frame = 1
        end
        avg = time_average(sol, first_frame)
        output["average"] = frame_dict(avg, 1)
    end

    if save_time_resolved
        output["frames"] = [frame_dict(sol, i) for i in eachindex(sol.frames)]
    end

    return OrderedDict(
        "serialization_version" => SERIALIZATION_VERSION,
        "input" => OrderedDict(
            "config" => serialize(sol.config),
            "simulation" => serialize(sol.simulation),
            "postprocess" => serialize(sol.postprocess),
        ),
        "output" => output,
    )
end

"""
    $(TYPEDSIGNATURES)
Write `sol` to `file`, if `file` is a JSON file. Any NaN or Inf values in the solution will be replaced with zero.

## Mandatory arguments
- `file`: the file to which we write the solution
- `sol`: the `Solution` object to be written

## Optional keyword args
- `average_start_time` = -1: the time at which averaging begins. If < 0, no averaged output is written.
- `save_time_resolved` = true: Whether to save all frames of the simulation. If `false`, no time-resolved output is written.
"""
Base.@nospecializeinfer function write_to_json(
        file::String, @nospecialize(sol::Solution);
        average_start_time::AbstractFloat = -1.0, save_time_resolved::Bool = true,
    )

    ext = splitext(file)[2]
    if lowercase(ext) != ".json"
        throw(ArgumentError("$(file) is not a JSON file."))
    end

    output = serialize_sol(sol; average_start_time, save_time_resolved)

    # Write output dictionary to file, replacing NaN and Inf values with zeros for JSON compliance.
    open(file, "w") do f
        JSON.write_json(f, output, replace_inf = true, replace_nan = true)
    end

    return nothing
end
