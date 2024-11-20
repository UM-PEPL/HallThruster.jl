serialize(x::Number) = x
deserialize(::Type{T}, x) where {T <: Number} = T(x)

serialize(x::AbstractString) = x
deserialize(::Type{T}, x) where {T <: AbstractString} = T(x)

serialize(x::AbstractVector) = serialize.(x)
deserialize(::Type{T}, x) where {T <: AbstractVector} = deserialize.(eltype(T), x)
