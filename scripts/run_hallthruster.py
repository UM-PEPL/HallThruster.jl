import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any

def run_hallthruster_jl(
    input: dict[str, Any] | str | Path,
    jl_environment: str = ""
    ):
    """Run a single HallThruster.jl simulation for a given set of inputs

    :param json_input: either a dictionary containing `config`, `simulation`, and `postprocess` options for 
            HallThruster.jl, or a string/Path containing a path to a JSON file with those inputs.
            See the HallThruster.jl documentation for the required keys.
    :param jl_environment: The julia environment containing HallThruster.jl. Defaults to global Julia environment.

    :returns: `dict` of `Hallthruster.jl` outputs: `I_B0`, `I_d`, `T`, `eta_c`, `eta_m`, `eta_v`, and `u_ion` for ion
              beam current (A), discharge current (A), thrust (N), current efficiency, mass efficiency, voltage
              efficiency, and singly-charged ion velocity profile (m/s). The specific outputs depend on the settings
              provided in the `postprocess` dict in the input.
              If `postprocess['output_file']` is present, this function will also write the requested outputs and
              restart information to that file.
    """

    # Read JSON input from file if path provided
    if isinstance(input, str) or isinstance(input, Path):
        with open(input) as fp:
            json_input : dict[str, Any] = json.load(fp)
    else:
        json_input = input

    # Get output file path. If one not provided, create a temporary
    if 'output_file' in json_input['postprocess']:
        output_file = Path(json_input['postprocess'].get('output_file'))
    else:
        # Create temporary
        fd_out = tempfile.NamedTemporaryFile(suffix=".json", prefix="hallthrusterjl_")
        output_file = Path(fd_out.name)
        fd_out.close()
        json_input['postprocess']['output_file'] = str(output_file)

    # Dump input to temporary file
    fd = tempfile.NamedTemporaryFile(
        suffix=".json", prefix="hallthrusterjl_",
        mode="w", delete=False, encoding="utf-8"
    )
    input_file = fd.name
    json.dump(json_input, fd, ensure_ascii=False, indent=4)
    fd.close()

    # Run HallThruster.jl on input file
    try:
        subprocess.run([
            'julia',
            f'--project=${jl_environment}',
            '-e'
            f'using HallThruster; HallThruster.run_simulation("{input_file}")'
        ])
    finally:
        # Delete temporary input file
        os.unlink(input_file)

    # Load output data
    with open(output_file) as fp:
        output_data = json.load(fp)

    return output_data
