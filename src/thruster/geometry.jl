@public Geometry1D

"""
$(TYPEDEF)

Describes the geometry of a Hall thruster discharge channel.

# Fields
$(TYPEDFIELDS)
"""
struct Geometry1D
    """
    The discharge channel length, in meters
    """
    channel_length::Float64
    """
    The inner radius of the discharge channel, in meters
    """
    inner_radius::Float64
    """
    The outer radius of the discharge channel, in meters
    """
    outer_radius::Float64
    """
    The discharge channel cross-sectional area, computed from `inner_radius` and `outer_radius`
    """
    channel_area::Float64
    function Geometry1D(; channel_length, inner_radius, outer_radius)
        channel_length = convert_to_float64(channel_length, units(:m))
        inner_radius = convert_to_float64(inner_radius, units(:m))
        outer_radius = convert_to_float64(outer_radius, units(:m))
        A_ch = channel_area(outer_radius, inner_radius)
        return new(channel_length, inner_radius, outer_radius, A_ch)
    end
end

"""
    channel_area(outer_radius, inner_radius)
Compute the cross-sectional area of a Hall thruster channel from its dimensions
"""
@inline channel_area(outer_radius, inner_radius) = π * (outer_radius^2 - inner_radius^2)
@inline channel_area(geometry::Geometry1D) = geometry.channel_area

"""
    channel_perimeter(outer_radius, inner_radius)
Compute the perimeteter of the thruster channel, equal to the sum of the inner and outer circumferences
"""
@inline channel_perimeter(outer_radius, inner_radius) = 2π * (outer_radius + inner_radius)
@inline channel_perimeter(geometry::Geometry1D) = channel_perimeter(
    geometry.outer_radius, geometry.inner_radius,
)

"""
    channel_width(outer_radius, inner_radius)
Compute the thruster channel width
"""
@inline channel_width(outer_radius, inner_radius) = outer_radius - inner_radius
@inline channel_width(geometry::Geometry1D) = channel_width(
    geometry.outer_radius, geometry.inner_radius,
)

#=============================================================================
 Serialization
==============================================================================#
Serialization.exclude(::Type{Geometry1D}) = [:channel_area]
