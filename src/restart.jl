"""
    write_restart(path::AbstractString, sol)

Write a restart file to `path``.

This can be reloaded to resume a simulation. The filetype can be anything supported by FileIO, though JLD2 is preferred.
"""
function write_restart(path::AbstractString, sol)
    save(path, Dict(
        "u" =>  sol.u[end],
        "params" => sol.params
    ))
end

"""
    read_restart(path::AbstractString)

Load a restart file from `path`.

The filetype can be anything supported by FileIO, though JLD2 is preferred.
"""
function read_restart(path::AbstractString)
    dict = load(path)
    u, params = dict["u"], dict["params"]
    ncells = length(params.z_cell)-2
    grid = Grid1D(
        ncells,
        params.z_edge,
        params.z_cell,
        params.cell_volume
    )
    B = params.cache.B
    fluid_ranges = params.fluid_ranges

    return u, grid, B, fluid_ranges
end

"""
    initialize_from_restart!(U, cache, restart_file, grid, fluid_ranges)
Load a solution from `restart_file`.

If the loaded solution has the same grid as that passed into this function,
it overwrites U and cache.B with the loaded values. Otherwise, it figures out the difference between the two simulations
and interpolates the old one  onto the new one.
"""
function initialize_from_restart!(U, cache, restart_file, grid, fluid_ranges)
    U_restart, grid_restart, B_restart, fluid_ranges_restart = read_restart(restart_file)

    lf_restart = fluid_ranges_restart[end][end]

    if grid_restart == grid
        cache.B .= B_restart

        for i in 1:length(grid.cell_centers)

            # Handle differing numbers of ion charge states and different conservation laws between
            # current simulation and restart
            for (r, r_restart) in zip(fluid_ranges, fluid_ranges_restart)
                for (f, f_restart) in zip(r, r_restart)
                    U[f, i] = U_restart[f_restart, i]
                end
            end

            # Electron quantities do not need any special handling
            @views U[lf+1:end, i] .= U_restart[lf_restart+1:end]
        end
    else
        # interpolate old solution to new grid
        for (r, r_restart) in zip(fluid_ranges, fluid_ranges_restart)
            for (f, f_restart) in zip(r, r_restart)
                itp = LinearInterpolation(grid_restart.cell_centers, U_restart[f_restart, :])
                U[f, :] .= itp.(grid.cell_centers)
            end
        end
        nvars = size(U, 1)
        nvars_restart = size(U_restart, 1)
        for (i, i_restart) in zip(lf+1:nvars, lf_restart+1:nvars_restart)
            itp = LinearInterpolation(grid_restart.cell_centers, U_restart[i_restart, :])
            U[i, :] .= itp.(grid.cell_centers)
        end
        compute_bfield!(cache.B, B_restart, grid_restart)
    end
end

