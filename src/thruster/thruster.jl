@public Thruster

"""
    $(TYPEDEF)

Defines a Hall thruster from a name, geometry, magnetic field.
The thruster may also be shielded, in which case the wall electron temperature is assumed to be equal to the anode electron temperature inside the discharge channel.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct Thruster
    """
    Name of the thruster. Not used during the simulation, but useful for certain cases.
    """
    name::String = "noname"
    """
    The thruster geometry
    """
    geometry::Geometry1D   
    """
    Contains a magnetic field at discrete axial locations
    """
    magnetic_field::MagneticField
    """
    Whether the thruster is magnetically-shielded
    """
    shielded::Bool = false
end

@inline channel_area(thruster::Thruster) = thruster.geometry.channel_area

@inline channel_perimeter(thruster::Thruster) = channel_perimeter(thruster.geometry)

@inline channel_width(thruster::Thruster) = channel_width(thruster.geometry)
