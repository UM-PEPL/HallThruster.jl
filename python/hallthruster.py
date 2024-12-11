import json
import os
from pathlib import Path
import subprocess
import tempfile

def run_simulation(
        input: dict | str | Path,
        jl_env: str | Path | None = None,
        jl_script: str | Path | None = None,
        **kwargs) -> dict:
    """Python wrapper for `HallThruster.run_simulation(json_input)` in Julia.

    :param json_input: either a dictionary containing `config`, `simulation`, and `postprocess` options for 
            HallThruster.jl, or a string/Path containing a path to a JSON file with those inputs.
    :param jl_env: The julia environment containing HallThruster.jl. Defaults to global Julia environment.
    :param jl_script: path to a custom Julia script to run. The script should accept the input json file path as
                      a command line argument. Defaults to just calling `HallThruster.run_simulation(input_file)`.
    :param kwargs: additional keyword arguments to pass to `subprocess.run` when calling the Julia script.

    :returns: `dict` of `Hallthruster.jl` outputs. The specific outputs depend on the settings
              provided in the `postprocess` dict in the input. If `postprocess['output_file']` is present,
              this function will also write the requested outputs and restart information to that file.
    """
    # Read JSON input from file if path provided
    if isinstance(input, str | Path):
        with open(input, 'r') as fp:
            json_input = json.load(fp)
    else:
        json_input = input

    tempfile_args = dict(suffix=".json", prefix="hallthruster_jl_", mode="w", delete=False, encoding="utf-8")

    # Get output file path. If one not provided, create a temporary
    temp_out = False
    if 'output_file' in json_input.get('postprocess', {}):
        output_file = Path(json_input['postprocess'].get('output_file'))
    elif 'output_file' in json_input.get('input', {}).get('postprocess', {}):
        output_file = Path(json_input['input']['postprocess'].get('output_file'))
    else:
        temp_out = True
        fd_out = tempfile.NamedTemporaryFile(**tempfile_args)
        output_file = Path(fd_out.name)
        fd_out.close()

        if json_input.get('input'):
            json_input['input'].setdefault('postprocess', {})
            json_input['input']['postprocess']['output_file'] = str(output_file.resolve())
        else:
            json_input.setdefault('postprocess', {})
            json_input['postprocess']['output_file'] = str(output_file.resolve())

    # Dump input to temporary file
    fd = tempfile.NamedTemporaryFile(**tempfile_args)
    input_file = fd.name
    json.dump(json_input, fd, ensure_ascii=False, indent=4)
    fd.close()

    # Run HallThruster.jl on input file
    if jl_script is None:
        cmd = ['julia', '--startup-file=no', '-e',
               f'using HallThruster; HallThruster.run_simulation(raw"{input_file}")']
    else:
        cmd = ['julia', '--startup-file=no', '--',
               str(Path(jl_script).resolve()), input_file]

    if jl_env is not None:
        if Path(jl_env).exists():
            cmd.insert(1, f'--project={Path(jl_env).resolve()}')
        else:
            raise ValueError(f"Could not find Julia environment {jl_env}. Please create it first. "
                             f"See https://github.com/JANUS-Institute/HallThrusterPEM/blob/main/scripts/install_hallthruster.py")

    try:
        subprocess.run(cmd, **kwargs)
    finally:
        # Delete temporary input file
        os.unlink(input_file)

    # Load output data
    with open(output_file, 'r') as fp:
        output_data = json.load(fp)

    if temp_out:
        os.unlink(output_file)
        if d := output_data.get('postprocess'):
            if 'output_file' in d:
                del d['output_file']
        if d := output_data.get('input'):
            if d2 := d.get('postprocess'):
                if 'output_file' in d2:
                    del d2['output_file']

    return output_data
