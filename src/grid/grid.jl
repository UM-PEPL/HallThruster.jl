struct Grid1D
    num_cells::Int
    edges::Vector{Float64}
    cell_centers::Vector{Float64}
    dz_edge::Vector{Float64}
    dz_cell::Vector{Float64}
end

"""
    Grid1D(geometry, z_edge)
Given 1-D edge coordinates and thruster geometry, compute cell centers and cell volumes for grid
"""
function Grid1D(z_edge)
    # Compute domain
    ncells = length(z_edge) - 1

    # generate cell center coordinates
    z_cell = [0.5 * (z_edge[i + 1] + z_edge[i]) for i in 1:ncells]

    # add ghost cells on left and right boundaries
    z_ghost_left = z_cell[1] - (z_edge[2] - z_edge[1])
    z_ghost_right = z_cell[end] + (z_edge[end] - z_edge[end - 1])
    z_cell = [z_ghost_left; z_cell; z_ghost_right]

    # calculate grid spacing
    dz_edge = @. @views z_cell[2:end] - z_cell[1:(end - 1)]
    dz_cell = @. @views z_edge[2:end] - z_edge[1:(end - 1)]
    dz_cell = [dz_cell[1]; dz_cell; dz_cell[end]]

    # Add 2 to cell count to account for ghost cells
    return Grid1D(ncells + 2, z_edge, z_cell, dz_edge, dz_cell)
end

@inline interior_cells(arr) = 2:(length(arr) - 1)

@inline function face_area(grid, area, face)
    weight = (grid.edges[face] - grid.cell_centers[face]) / grid.dz_edge[face]
    return lerp(area[face], area[face + 1], weight)
end
