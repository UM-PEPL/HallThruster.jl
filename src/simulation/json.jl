@public write_to_json, run_simulation

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

    obj = JSON3.read(read(json_file), allow_inf = true)

    # Read config and sim params from file
    if haskey(obj, "input")
        input = obj.input
    else
        input = obj
    end

    cfg = deserialize(Config, input.config)
    sim = deserialize(SimParams, input.simulation)

    postprocess::Union{Postprocess, Nothing} = nothing
    if haskey(input, "postprocess") && haskey(input.postprocess, "output_file") &&
            !isempty(input.postprocess.output_file)
        postprocess = deserialize(Postprocess, input.postprocess)
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
function frame_dict(sol::Solution, frame::Integer)
    # TODO: multiple propellants
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
    d["z"] = sol.grid.cell_centers
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

    if length(sol.config.propellants) == 1
        symbol = sol.config.propellants[1].gas.short_name
        d["nn"] = f.neutrals[symbol].n
        d["ni"] = [ion.n for ion in f.ions[symbol]]
        d["ui"] = [ion.u for ion in f.ions[symbol]]
        d["niui"] = [ion.nu for ion in f.ions[symbol]]
    end

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

    return d
end

"""
    $(TYPEDSIGNATURES)
Convert `sol` to an `OrderedDict`, containing both the inputs used to run the simulation
and any requested outputs.
This function is used to convert a `Solution` to a format suitable for writing to an output file.
"""
function serialize_sol(
        sol::Solution; average_start_time::AbstractFloat = -1, save_time_resolved::Bool = true,
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
Write `sol` to `file`, if `file` is a JSON file.

## Mandatory arguments
- `file`: the file to which we write the solution
- `sol`: the `Solution` object to be written

## Optional keyword args
- `average_start_time` = -1: the time at which averaging begins. If < 0, no averaged output is written.
- `save_time_resolved` = true: Whether to save all frames of the simulation. If `false`, no time-resolved output is written.
"""
function write_to_json(
        file::String, sol::Solution;
        average_start_time::AbstractFloat = -1.0, save_time_resolved::Bool = true,
    )

    ext = splitext(file)[2]
    if lowercase(ext) != ".json"
        throw(ArgumentError("$(file) is not a JSON file."))
    end

    output = serialize_sol(sol; average_start_time, save_time_resolved)

    open(file, "w") do f
        JSON3.write(f, output, allow_inf = true)
    end

    return nothing
end
