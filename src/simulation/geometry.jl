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

"""
    Thruster
Struct containing information about a Hall thruster. This includes a `name`, `geometry` (a `Geometry1D` object), `magnetic_field` (radial magnetic field along
centerline, a function which takes z in meters and outputs B in Tesla), and a `shielded` (a flag indicating whether the thruster is magnetically-shielded).

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Thruster{B}
    name::String = "noname"
    geometry::Geometry1D   # A Geometry1D object containing geometric information about the thruster
    magnetic_field::B      # A function which takes z in meters and outputs B in Tesla
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
    geometry.outer_radius, geometry.inner_radius,)
@inline channel_perimeter(thruster::Thruster) = channel_perimeter(thruster.geometry)

"""
    channel_width(outer_radius, inner_radius)
Compute the thruster channel width
"""
@inline channel_width(outer_radius, inner_radius) = outer_radius - inner_radius
@inline channel_width(geometry::Geometry1D) = channel_width(
    geometry.outer_radius, geometry.inner_radius,)
@inline channel_width(thruster::Thruster) = channel_width(thruster.geometry)

const geometry_SPT_100 = Geometry1D(
    inner_radius = 0.0345,
    outer_radius = 0.05,
    channel_length = 0.025,
)

function B_field_SPT_100(B_max, L_ch, z) #same in Landmark and in FFM model Hara
    B = if z < L_ch
        B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2) #for SPT_100
    else
        B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
    end
    return B
end

const SPT_100 = Thruster(
    name = "SPT-100",
    geometry = geometry_SPT_100,
    magnetic_field = z -> B_field_SPT_100(0.015, geometry_SPT_100.channel_length, z),
    shielded = false,
)
