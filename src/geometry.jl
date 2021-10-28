Geometry1D = @NamedTuple begin
    domain::Tuple{Float64, Float64}
    channel_length::Float64
    inner_radius::Float64
    outer_radius::Float64
end

"""
    channel_area(geometry::Geometry1D)
Compute the area of the Hall thruster channel from the given Geometry1D object
"""
channel_area(geometry::Geometry1D) =
    channel_area(geometry.outer_radius, geometry.inner_radius)

"""
    channel_area(outer_radius, inner_radius, length)
Compute the area of a Hall thruster channel from its dimensions
"""
channel_area(outer_radius, inner_radius) = π * (outer_radius^2 - inner_radius^2)


"""
    generate_grid(geometry, ncells)
Generate a one-dimensional uniform grid on the domain specified in the geomety. Returns coordinates
of cell centers (plus ghost cells) as well as cell interface/edges
"""

function generate_grid(geometry, ncells)
    # generate cell interface coordinates
    z_edge = LinRange(geometry.domain[1], geometry.domain[2], ncells+1)

    # generate cell center coordinates
    z_cell = [0.5 * (z_edge[i+1] + z_edge[i]) for i in 1:ncells]

    # add ghost cells on left and right boundaries
    z_cell = [z_edge[1]; z_cell; z_edge[end]]
    return z_cell, z_edge
end

"""
    SPT_100
Geometry of the SPT_100 thruster
"""

const SPT_100 = (
    domain = (0.0, 0.05),
    channel_length = 0.025,
    inner_radius = 0.0345,
    outer_radius = 0.05
)