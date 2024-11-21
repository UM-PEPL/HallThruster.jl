# Trait-based approach to serialization, inspired by StructTypes.jl
module Serialization

using Base.Iterators: Iterators
using OrderedCollections: OrderedDict

abstract type SType end
struct Null <: SType end
struct Boolean <: SType end
struct Number <: SType end
struct String <: SType end
struct ArrayType <: SType end
struct Enum <: SType end

abstract type Composite <: SType end
struct TaggedUnion <: Composite end
struct Struct <: Composite end

SType(::Type{T}) where {T}                   = Struct()
SType(::Type{Bool})                          = Boolean()
SType(::Type{Nothing})                       = Null()
SType(::Type{T}) where {T <: Base.Number}    = Number()
SType(::Type{T}) where {T <: AbstractString} = String()
SType(::Type{T}) where {T <: AbstractVector} = ArrayType()

exclude(::Any) = Symbol[]
options(::Any) = NamedTuple()
typetag(::Any) = "type"
iterate_fields(::Type{T}) where {T} = iterate_fields(SType(T), T)

function iterate_fields(::Composite, ::Type{T}) where {T}
    excluded = exclude(T)
    names = fieldnames(T)
    types = fieldtypes(T)
    return ((name, type) for (name, type) in zip(names, types) if !(name in excluded))
end

serialize(x::T) where {T} = serialize(SType(T), x)
deserialize(::Type{T}, x) where {T} = deserialize(SType(T), T, x)

serialize(::Null, x) = :null
deserialize(::Null, ::Type{Null}, x) = nothing

serialize(::Boolean, x) = x ? :true : :false
deserialize(::Boolean, ::Type{Bool}, x) = x == :true ? true : false

serialize(::Number, x) = x
deserialize(::Number, ::Type{T}, x) where {T} = T(x)

serialize(::String, x) = string(x)
deserialize(::String, ::Type{T}, x) where {T} = T(x)

serialize(::ArrayType, x) = serialize.(x)
deserialize(::ArrayType, ::Type{T}, x) where {T} = deserialize.(eltype(T), x)

function serialize(::Struct, x::T) where {T}
    return OrderedDict(
        string(field) => serialize(getfield(x, field))
    for (field, _) in iterate_fields(T)
    )
end

function deserialize(::Struct, ::Type{T}, dict::AbstractDict) where {T}
    args = NamedTuple(
        Symbol(field) => deserialize(fieldtype(T, Symbol(field)), dict[field])
    for field in keys(dict)
    )
    return T(; args...)
end

function serialize(::TaggedUnion, x::T) where {T}
    opts = options(T)
    for k in keys(opts)
        if T <: opts[k]
            pairs = (
                string(field) => serialize(getfield(x, field))
            for (field, _) in iterate_fields(T)
            )
            return OrderedDict(typetag(T) => k, pairs...)
        end
    end
    throw(ArgumentError("Invalid type $(T). Valid options are $(keys(opts))"))
end

function deserialize(::TaggedUnion, ::Type{T}, dict::AbstractDict) where {T}
    tag = dict[typetag(T)]
    subtype = getfield(options(T), Symbol(tag))

    pairs = NamedTuple(
        field => deserialize(type, dict[field])
    for (field, type) in iterate_fields(subtype)
    )
    return subtype(; pairs...)
end

function serialize(::Enum, x::T) where {T}
    opts = options(T)
    for k in keys(opts)
        if opts[k] == x
            return k
        end
    end
    throw(ArgumentError("Invalid value $(x) for type $(T). \
                        Valid options are $(keys(opts))"))
end

function deserialize(::Enum, ::Type{T}, x) where {T}
    name = Symbol(x)
    return getfield(options(T), name)
end

end
