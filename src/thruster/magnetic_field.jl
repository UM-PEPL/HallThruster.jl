"""
    $(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct MagneticField
    file::String
    z::Vector{Float64}
    B::Vector{Float64}
    function MagneticField(file, z = Float64[], B = Float64[])
        return new(file, z, B)
    end
end

function MagneticField(; file = "", z = Float64[], B = Float64[])
    return MagneticField(file, z, B)
end

function load_magnetic_field!(magnetic_field::MagneticField; include_dirs)
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
