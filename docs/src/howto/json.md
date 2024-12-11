# Use JSON for input and output

```@meta
CurrentModule = HallThruster
```

## Running simulations from JSON input.

In addition to the usual way of running simulations using a `Config` and `SimParams`, `HallThruster` supports the use of JSON for input and output.
This can be achieved by calling `run_simulation("input.json")`, where `"input.json"` could be any json file.
The contents of this file should be laid out in one of two ways.
The first way has up to three keys -- two mandatory (`"config"` and `"simulation"`) and one optional (`"postprocess"`).
These exactly mirror the [`Config`](@ref), [`SimParams`](@ref), and [`Postprocess`](@ref) structs that are normally passed to `run_simulation` and should be relatively self-explanatory.

Below, we show an example of this type of JSON input.

```json
{
    "config": {
        "thruster": {
            "name": "SPT-100",
            "geometry": { "inner_radius": 0.0345, "outer_radius": 0.05, "channel_length": 0.025 },
            "magnetic_field": { file = "bfield_spt100.csv" }
        },
        "propellant": "Krypton",
        "discharge_voltage": 300.0,
        "anode_mass_flow_rate": 5e-6,
        "domain": [0.0, 0.08],
        "anom_model": {
            "type": "TwoZoneBohm",
            "c1": 0.00625,
            "c2": 0.0625
        }
    },
    "simulation": {
        "adaptive": true,
        "dt": 5e-9,
        "grid": {
            "type": "EvenGrid",
            "num_cells": 200,
        },
        "duration": 1e-3,
        "num_save": 1000,
    },
    "postprocess": {
        "output_file": "output.json"
        "save_time_resolved": false,
        "average_start_time": 5e-4
    }
}
```
In the second way, the JSON file has a top level field `input`, with members `config`, `simulation`, and (optionally) `postprocess`, e.g.

```json
{
    "input": {
        "config": {...},
        "simulation": {...},
        "postprocess": {...}
    }
}
```
where the contents of these keys are exactly the same as in the first method.

In both cases, the field names and types of these inputs are exactly the same as in the corresponding `HallThruster` types.
Note the `"type"` field for `config.anom_model` and `simulation.grid`, which precedes the other fields for that type.
This pattern is also used for `config.conductivity_model` and `config.wall_loss_model`, should you wish to provide those.
Note that custom anomalous transport models and propellants are not supported using the JSON interface at this time.

## Writing output files

If `postprocess` is provided and `postprocess.output_file` is not empty, `HallThruster` will write an output JSON file to that file.
The output file contains two top-level fields: `input` and `output`.
The `input` field reproduces the inputs used to run the simulation, exactly as described above.
The `output` field has at most four fields: `retcode`, `error`, `fields`, and `average`, e.g.

```json
{
    "input": {...},
    "output": {
        "retcode": "success",
        "error": "",
        "average": {...},
        "frames": [
            {...}, {...}, ...
        ]
    }
}
```

The `output_field` somewhat mirrors the [`Solution`](@ref) object.
The `retcode` gives the simulation status (one of `"success"`, `"error"`, or `"failure"`).
If `retcode` is `error`, the `error` will contain a string with the error that occurred.
If `postprocess.average_start_time` is greater than or equal to zero, the `average` field contains the output of `HallThruster.time_average(sol, average_start_time)`.
Finally, if `postprocess.save_time_resolved` is `true`, then `output.frames` contains `simulation.num_save` frames.

You can also manually write a `Solution` to a JSON file using the [`write_to_json`](@ref) function.

### Frame format

The frames in the output JSON file are laid out similarly to those in the `Solution` struct, with some additions.
In addition to plasma properties, each frame also stores its time, as well as thrust, discharge current, and component efficiencies.

## Restarts
!!! warning "Interface not finalized"
    The restart interface is not finalized and is subject to change before v1.0.
