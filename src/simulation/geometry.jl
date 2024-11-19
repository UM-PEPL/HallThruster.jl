"""
    Geometry1D
Struct containing information about Hall thruster geometry.
Required fields are `channel_length`, `inner_radius`, and `outer_radius`, all in meters.
Contains a fourth field, `channel_area`, which is computed from the above three.

# Fields
$(TYPEDFIELDS)
"""

export EvenGrid, UnevenGrid

struct Geometry1D
    channel_length::Float64
    inner_radius::Float64
    outer_radius::Float64
    channel_area::Float64
    function Geometry1D(; channel_length, inner_radius, outer_radius)
        if channel_length isa Quantity
            channel_length = ustrip(uconvert(u"m", channel_length))
        end

        if inner_radius isa Quantity
            inner_radius = ustrip(uconvert(u"m", inner_radius))
        end

        if outer_radius isa Quantity
            outer_radius = ustrip(uconvert(u"m", outer_radius))
        end

        A_ch = channel_area(outer_radius, inner_radius)

        return new(channel_length, inner_radius, outer_radius, A_ch)
    end
end

function Geometry1D(channel_length, inner_radius, outer_radius, channel_area = nothing)
    return Geometry1D(; channel_length, inner_radius, outer_radius)
end

"""
    Magnetic field

# Fields
$(TYPEDFIELDS)
"""
struct MagneticField
    file::String
    z::Vector{Float64}
    B::Vector{Float64}
end

function MagneticField(file::String)
    data = readdlm(file, ',')
    return MagneticField(file, data[:, 1], data[:, 2])
end
@inline MagneticField(file::String, ::Nothing, ::Nothing) = MagneticField(file)
@inline MagneticField(z, B) = MagneticField("", z, B)
@inline MagneticField(::Nothing, z, B) = MagneticField(z, B)

"""
    Thruster
Struct containing information about a Hall thruster. This includes a `name`, `geometry` (a `Geometry1D` object), `magnetic_field` (radial magnetic field along
centerline, a function which takes z in meters and outputs B in Tesla), and a `shielded` (a flag indicating whether the thruster is magnetically-shielded).

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Thruster
    name::String = "noname"
    geometry::Geometry1D   # A Geometry1D object containing geometric information about the thruster
    magnetic_field::MagneticField      # A function which takes z in meters and outputs B in Tesla
    shielded::Bool         # Whether the thruster is magnetically-shielded
end

"""
    channel_area(outer_radius, inner_radius)
Compute the cross-sectional area of a Hall thruster channel from its dimensions
"""
@inline channel_area(outer_radius, inner_radius) = π * (outer_radius^2 - inner_radius^2)
@inline channel_area(geometry::Geometry1D) = geometry.channel_area
@inline channel_area(thruster::Thruster) = thruster.geometry.channel_area

"""
    channel_perimeter(outer_radius, inner_radius)
Compute the perimeteter of the thruster channel, equal to the sum of the inner and outer circumferences
"""
@inline channel_perimeter(outer_radius, inner_radius) = 2π * (outer_radius + inner_radius)
@inline channel_perimeter(geometry::Geometry1D) = channel_perimeter(
    geometry.outer_radius, geometry.inner_radius)
@inline channel_perimeter(thruster::Thruster) = channel_perimeter(thruster.geometry)

"""
    channel_width(outer_radius, inner_radius)
Compute the thruster channel width
"""
@inline channel_width(outer_radius, inner_radius) = outer_radius - inner_radius
@inline channel_width(geometry::Geometry1D) = channel_width(
    geometry.outer_radius, geometry.inner_radius)
@inline channel_width(thruster::Thruster) = channel_width(thruster.geometry)

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
        grid.density, domain[1], domain[2], geometry.channel_length, grid.ncells + 1
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

const geometry_SPT_100 = Geometry1D(;
    inner_radius = 0.0345, outer_radius = 0.05, channel_length = 0.025
)

function B_field_SPT_100() # same in Landmark and in FFM model Hara
    L_ch = geometry_SPT_100.channel_length
    B_max = 0.015 # T
    zs = LinRange(0, 4 * L_ch, 256)
    Bs = zeros(length(zs))
    for (i, z) in enumerate(zs)
        Bs[i] = if z < L_ch
            B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2)
        else
            B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
        end
    end
    return MagneticField("SPT-100 Default", zs, Bs)
end

const SPT_100 = Thruster(;
    name = "SPT-100",
    geometry = geometry_SPT_100,
    magnetic_field = B_field_SPT_100(),
    shielded = false
)
