serialize(x::Any) = x
deserialize(::Type{T}, arg) where {T} = T(arg)
deserialize(::Type{T}, d::AbstractDict) where {T} = T(d)
