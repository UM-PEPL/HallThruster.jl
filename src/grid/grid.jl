struct Grid1D
    num_cells::Int
    edges::Vector{Float64}
    cell_centers::Vector{Float64}
end

"""
    Grid1D(geometry, z_edge)
Given 1-D edge coordinates and thruster geometry, compute cell centers and cell volumes for grid
"""
function Grid1D(z_edge)
    # Compute domain
    domain = (z_edge[1], z_edge[end])
    ncells = length(z_edge) - 1

    # generate cell center coordinates
    z_cell = [0.5 * (z_edge[i + 1] + z_edge[i]) for i in 1:ncells]

    # add ghost cells on left and right boundaries
    z_cell = [z_edge[1]; z_cell; z_edge[end]]

    return Grid1D(ncells, z_edge, z_cell)
end

function grid_spacing(grid)
    z_cell = grid.cell_centers
    z_edge = grid.edges

    # Fill up cell lengths and magnetic field vectors
    Δz_cell = zeros(length(z_cell))
    for i in eachindex(z_cell)
        if firstindex(z_cell) < i < lastindex(z_cell)
            Δz_cell[i] = z_edge[right_edge(i)] - z_edge[left_edge(i)]
        elseif i == firstindex(z_cell)
            Δz_cell[i] = z_edge[begin + 1] - z_edge[begin]
        elseif i == lastindex(z_cell)
            Δz_cell[i] = z_edge[end] - z_edge[end - 1]
        end
    end

    Δz_edge = @. @views z_cell[2:end] - z_cell[1:(end - 1)]

    return Δz_cell, Δz_edge
end
