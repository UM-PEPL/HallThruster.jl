Geometry1D = @NamedTuple begin
    domain::Tuple{Float64,Float64}
    channel_length::Float64
    inner_radius::Float64
    outer_radius::Float64
end

struct Grid1D
    ncells::Int64
    edges::Vector{Float64}
    cell_centers::Vector{Float64}
    cell_volume::Float64
end

"""
    channel_area(geometry::Geometry1D)
Compute the area of the Hall thruster channel from the given Geometry1D object
"""
function channel_area(geometry::Geometry1D)
    return channel_area(geometry.outer_radius, geometry.inner_radius)
end

"""
    channel_area(outer_radius, inner_radius, length)
Compute the area of a Hall thruster channel from its dimensions
"""
channel_area(outer_radius, inner_radius) = π * (outer_radius^2 - inner_radius^2)

"""
    generate_grid(geometry, ncells)
Generate a one-dimensional uniform grid on the domain specified in the geomety. Returns number of cells, coordinates
of cell centers (plus ghost cells face coordinates), interface/edges and volume of a cell for number density calculations. 
"""
function generate_grid(geometry, ncells, domain = nothing)
    domain = domain === nothing ? geometry.domain : domain

    # generate cell interface coordinates
    z_edge = LinRange(domain[1], domain[2], ncells+1)

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
    SPT_100
Geometry of the SPT_100 thruster
"""

const SPT_100 = (domain=(0.0, 0.05), channel_length=0.025, inner_radius=0.0345,
                 outer_radius=0.05)

const SPT_100_1 = (domain=(0.0, 0.08), channel_length=0.025, inner_radius=0.0345,
                 outer_radius=0.05)