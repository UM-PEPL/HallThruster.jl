function B_field_SPT_100(B_max, L_ch, z) #same in Landmark and in FFM model Hara
    B = if z < L_ch
        B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2) #for SPT_100
    else
        B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
    end
    return B
end

struct Geometry1D
    channel_length::Float64
    inner_radius::Float64
    outer_radius::Float64
    channel_area::Float64
    function Geometry1D(;channel_length, inner_radius, outer_radius)
        A_ch = channel_area(outer_radius, inner_radius)
        return new(channel_length, inner_radius, outer_radius, A_ch)
    end
end

Base.@kwdef struct Thruster{B}
    name::String = "noname"
    geometry::HallThruster.Geometry1D
    magnetic_field::B
    shielded::Bool
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
@inline channel_area(geometry::Geometry1D) = geometry.channel_area

@inline channel_area(thruster::Thruster) = thruster.geometry.channel_area

"""
    channel_area(outer_radius, inner_radius, length)
Compute the area of a Hall thruster channel from its dimensions
"""
@inline channel_area(outer_radius, inner_radius) = π * (outer_radius^2 - inner_radius^2)

"""
    generate_grid(geometry, ncells)
Generate a one-dimensional uniform grid on the domain specified in the geomety. Returns number of cells, coordinates
of cell centers (plus ghost cells face coordinates), interface/edges and volume of a cell for number density calculations.
"""
function generate_grid(geometry, ncells, domain)
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

const geometry_SPT_100 = Geometry1D(
    inner_radius = 0.0345,
    outer_radius = 0.05,
    channel_length = 0.025
)

const SPT_100 = Thruster(
    name = "SPT-100",
    geometry = geometry_SPT_100,
    magnetic_field = B_field_SPT_100 $ (0.015, geometry_SPT_100.channel_length),
    shielded = false
)
