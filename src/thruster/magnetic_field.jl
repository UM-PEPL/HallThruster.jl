"""
    $(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct MagneticField
    file::String = ""
    z::Vector{Float64} = []
    B::Vector{Float64} = []
end

function MagneticField(file::String)
    data = readdlm(file, ',')
    return MagneticField(file, data[:, 1], data[:, 2])
end
