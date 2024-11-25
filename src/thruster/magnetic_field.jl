"""
    $(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
struct MagneticField
    file::String
    z::Vector{Float64}
    B::Vector{Float64}
    function MagneticField(file, z = Float64[], B = Float64[])
        if isempty(z) || isempty(B)
            data = readdlm(file, ',')
            return new(file, data[:, 1], data[:, 2])
        end
        return new(file, z, B)
    end
end

function MagneticField(; file = "", z = Float64[], B = Float64[])
    return MagneticField(file, z, B)
end
