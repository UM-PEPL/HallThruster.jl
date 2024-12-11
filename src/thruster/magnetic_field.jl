@public MagneticField, load_magnetic_field, load_magnetic_field!

"""
    $(TYPEDEF)

Specifies the radial magnetic field of a Hall thruster, measured along channel centerline.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MagneticField
    """
    The file where magnetic field information can be found. Can be left empty if `z` and `B` are explicitly specified.
    """
    file::String
    """
    The axial coordinates at which the magnetic field is known, in meters.
    """
    z::Vector{Float64}
    """
    The magnetic field at each point in `z`, measured in Teslas.
    """
    B::Vector{Float64}
    function MagneticField(file, z = Float64[], B = Float64[])
        return new(file, z, B)
    end
end

function MagneticField(; file = "", z = Float64[], B = Float64[])
    return MagneticField(file, z, B)
end

"""
	$(TYPEDSIGNATURES)
Given a magnetic field `field` with empty `z` or `B` fields, looks for a magnetic field at `field`.file.
If one is found in the present directory or in any of the provided `include_dirs`, sets `z` and `B` accordingly.
Throws an ArgumentError if one is not found.
"""
function load_magnetic_field!(magnetic_field::MagneticField; include_dirs = String[])::Nothing
    if isempty(magnetic_field.z) || isempty(magnetic_field.B)
        include_dirs = [include_dirs; pwd()]
        for dir in include_dirs
            file = joinpath(dir, magnetic_field.file)
            if ispath(file)
                data = readdlm(file, ',')
                magnetic_field.z = data[:, 1]
                magnetic_field.B = data[:, 2]
                return
            end
        end
        throw(ArgumentError("Magnetic field field $(magnetic_field.file) not found in dirs $(include_dirs)"))
    end
end

"""
	$(TYPEDSIGNATURES)
Given a path to a file, loads a `MagneticField` from the data in that file.
Looks in the present working directory and any additional directories passed to `include_dirs`.
"""
function load_magnetic_field(file::String; include_dirs = String[])::MagneticField
	field = MagneticField(;file);
	load_magnetic_field!(field; include_dirs)
	return field
end
