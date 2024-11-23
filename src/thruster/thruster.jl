"""
    Thruster
Struct containing information about a Hall thruster. This includes a `name`, `geometry` (a `Geometry1D` object), `magnetic_field` (radial magnetic field along
centerline, a function which takes z in meters and outputs B in Tesla), and a `shielded` (a flag indicating whether the thruster is magnetically-shielded).

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct Thruster
    name::String = "noname"
    geometry::Geometry1D   # A Geometry1D object containing geometric information about the thruster
    magnetic_field::MagneticField      # A function which takes z in meters and outputs B in Tesla
    shielded::Bool = false         # Whether the thruster is magnetically-shielded
end

@inline channel_area(thruster::Thruster) = thruster.geometry.channel_area

@inline channel_perimeter(thruster::Thruster) = channel_perimeter(thruster.geometry)

@inline channel_width(thruster::Thruster) = channel_width(thruster.geometry)
