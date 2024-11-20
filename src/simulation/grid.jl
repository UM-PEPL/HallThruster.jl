struct Grid1D
    ncells::Int64
    edges::Vector{Float64}
    cell_centers::Vector{Float64}
    cell_volume::Float64
end

struct HallThrusterGrid{F}
    ncells::Int
    density::F
end

"""
    EvenGrid(n)
Construct an evenly-spaced grid with n cells.
"""
EvenGrid(n) = HallThrusterGrid(n, Returns(1.0))

"""
    UnevenGrid(n, density = HallThruster.default_density)
Construct an unevenly-spaced grid according to provided density function. Defaults to twice as many grids inside of channel than outside.
Provided density functions must have a signature of (z, z0, z1, Lch) where z is the location,
(z0, z1) are the extents of the domain and Lch is the channel length
"""
UnevenGrid(n, f = default_density) = HallThrusterGrid(n, f)

function default_density(z, z0, z1, Lch)
    center = 1.5 * Lch
    width = 0.5 * Lch
    if (z < center)
        return 2.0
    else
        return 1.0 + exp(-((z - center) / (width))^2)
    end
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

"""
    Grid1D(geometry, z_edge)
Given 1-D edge coordinates and thruster geometry, compute cell centers and cell volumes for grid
"""
function Grid1D(geometry, z_edge)
    # Compute domain
    domain = (z_edge[1], z_edge[end])
    ncells = length(z_edge) - 1

    # generate cell center coordinates
    z_cell = [0.5 * (z_edge[i + 1] + z_edge[i]) for i in 1:ncells]

    # add ghost cells on left and right boundaries
    z_cell = [z_edge[1]; z_cell; z_edge[end]]

    # get cell area
    cell_volume = channel_area(geometry.outer_radius, geometry.inner_radius) *
                  abs(domain[2] - domain[1]) / ncells

    return Grid1D(ncells, z_edge, z_cell, cell_volume)
end

"""
    generate_grid(geometry, ncells)
Generate a one-dimensional uniform grid on the domain specified in the geomety. Returns number of cells, coordinates
of cell centers (plus ghost cells face coordinates), interface/edges and volume of a cell for number density calculations.
"""
function generate_grid(geometry, domain, grid::HallThrusterGrid{F}) where {F}
    # get domain
    domain = domain === nothing ? geometry.domain : domain

    # generate edge coordinates
    z_edge = nodes_from_density(
        grid.density, domain[1], domain[2], geometry.channel_length, grid.ncells + 1,
    )

    # compute centers and volumes and return
    return Grid1D(geometry, z_edge)
end

"""
    nodes_from_density(density, x0, x1, N)
Given bounds x0, x1, a number of points N, and a density function density(x), generate N nodes betweeen x0 and x1 spaced according
to the provided desity function using inverse CDF transformation.
"""
function nodes_from_density(density, x0, x1, Lch, N)
    xs = range(x0, x1; length = N)
    den = [density(x, x0, x1, Lch) for x in xs]
    cdf = cumsum(den)
    cdf_min, cdf_max = extrema(cdf)
    cdf_pts = range(cdf_min, cdf_max; length = N)
    xs_density = [HallThruster.interpolate(c, cdf, xs) for c in cdf_pts]
    return xs_density
end
