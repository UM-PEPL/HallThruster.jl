@public Postprocess, SimParams

"""
	$TYPEDEF

$(TYPEDFIELDS)
"""
mutable struct SimParams{C <: CurrentController}
    """
    A grid specifier, either `HallThruster.EvenGrid(n)` or `HallThruster.UnevenGrid(n)`, where n is the desired number of cells. See [Grid generation](@ref) for more information.
    """
    grid::GridSpec
    """
    The simulation base timestep in seconds. See [Timestepping](@ref) for more info. **Default:** 1e-9.
    """
    dt::Float64
    """
    The simulation duration in seconds. **Default:** 1e-3.
    """
    duration::Float64
    # Optional parameters
    """
    How many simulation frames to save in the `Solution` struct. **Default:** 1000
    """
    num_save::Int
    """
    Whether information such as the simulation run-time is printed to console. **Default:** `true`
    """
    verbose::Bool
    """
    Whether errors are printed to console in addition to being captured in the `Solution` struct. **Default:** `true`
    """
    print_errors::Bool
    """
    Whether to use adaptive timestepping. See [Timestepping](@ref) for more info. **Default:** `true`
    """
    adaptive::Bool
    """
    The CFL number used in adaptive timestepping. Maximum is 0.799. **Default:** 0.799
    """
    CFL::Float64
    """
    The minimum allowable timestep in adaptive timestepping, in seconds. **Default:** 1e-10
    """
    min_dt::Float64
    """
    The maximum allowable timestep in adaptive timestepping, in seconds. **Default:** 1e-7
    """
    max_dt::Float64
    """
    The maximum number of minimally-sized timesteps permitted in adaptive timestepping. **Default:** 100
    """
    max_small_steps::Int
    """
    Discharge current controller. **Default:** `HallThruster.NoController()`
    """
    current_control::C

    function SimParams(;
            grid::GridSpec,
            duration = 1.0e-3,
            dt = 1.0e-9,
            # Optional parameters
            num_save::Int = 1000,
            verbose::Bool = true,
            print_errors::Bool = true,
            adaptive::Bool = true,
            CFL::Float64 = 0.799,
            min_dt = 1.0e-10,
            max_dt = 1.0e-7,
            max_small_steps::Int = 100,
            current_control::C = NoController(),
        ) where {C <: CurrentController}
        return new{C}(
            grid,
            convert_to_float64(dt, units(:s)),
            convert_to_float64(duration, units(:s)),
            num_save,
            verbose,
            print_errors,
            adaptive,
            CFL,
            convert_to_float64(min_dt, units(:s)),
            convert_to_float64(max_dt, units(:s)),
            max_small_steps,
            current_control,
        )
    end
end


"""
$(TYPEDEF)
Contains postprocessing options for a given simulation.
When `run_simulation(config, sim_params; postprocess = Postprocess(...))` is called with a non-empty `output_file`, `HallThruster` will write the simulation results to a JSON file.
The results in the file will be transformed according to the fields.

## Fields
$(TYPEDFIELDS)
"""
@kwdef struct Postprocess
    """
    The file to which the output will be written. If empty, no output will be written.
    """
    output_file::String = ""
    """
    The time to begin averaging at. If less than zero, no averaged output will be written.
    """
    average_start_time::Float64 = -1
    """
    Whether time-resolved output will be saved.
    If true, each frame of the simulation will be written to the output file.
    """
    save_time_resolved::Bool = false
end
