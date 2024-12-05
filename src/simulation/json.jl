@public write_to_json, run_simulation

"""
    $(TYPEDSIGNATURES)
Run a simulation from a JSON input file.
If `postprocess` is set and `postprocess.output_file` is non-empty, output will be written
to `postprocess.output_file`.
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

    obj = JSON3.read(read(json_file))

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
    (; grid) = sol.params
    ncharge = sol.config.ncharge
    sv = sol.savevals[frame]
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
    d["z"] = grid.cell_centers
    d["nn"] = sv.nn
    d["ni"] = [sv.ni[Z, :] for Z in 1:ncharge]
    d["ui"] = [sv.ui[Z, :] for Z in 1:ncharge]
    d["niui"] = [sv.niui[Z, :] for Z in 1:ncharge]
    d["B"] = sol.params.cache.B
    d["ne"] = sv.ne
    d["ue"] = sv.ue
    d["V"] = sv.ϕ
    d["E"] = -sv.∇ϕ
    d["Tev"] = sv.Tev
    d["pe"] = sv.pe
    d["nu_en"] = sv.νen
    d["nu_ei"] = sv.νei
    d["nu_anom"] = sv.νan
    d["nu_class"] = sv.νc
    d["mobility"] = sv.μ
    d["channel_area"] = sv.channel_area
    return d
end

"""
    $(TYPEDSIGNATURES)
Convert `sol` to an `OrderedDict`, containing both the inputs used to run the simulation 
and any requested outputs. 
This function is used to convert a `Solution` to a format suitable for writing to an output file.
"""
function serialize_sol(
        sol::Solution; average_start_time::AbstractFloat = -1, save_time_resolved::Bool = true,)
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
        output["frames"] = [frame_dict(sol, i) for i in eachindex(sol.savevals)]
    end

    return OrderedDict(
        "input" => OrderedDict(
            "config" => serialize(sol.config),
            "simulation" => serialize(sol.params.simulation),
            "postprocess" => serialize(sol.params.postprocess),
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
        average_start_time::AbstractFloat = -1.0, save_time_resolved::Bool = true,)
    ext = splitext(file)[2]
    if lowercase(ext) != ".json"
        throw(ArgumentError("$(file) is not a JSON file."))
    end

    output = serialize_sol(sol; average_start_time, save_time_resolved)

    JSON3.write(file, output)
    return nothing
end
