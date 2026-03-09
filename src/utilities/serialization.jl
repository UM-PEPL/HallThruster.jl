# Trait-based approach to serialization, inspired by StructTypes.jl
module Serialization

using Base.Iterators: Iterators
using OrderedCollections: OrderedDict

abstract type SType end
struct Null <: SType end
struct Boolean <: SType end
struct NumberType <: SType end
struct StringType <: SType end
struct ArrayType <: SType end
struct TupleType <: SType end
struct EnumType <: SType end
struct DictType <: SType end

abstract type Composite <: SType end
struct TaggedUnion <: Composite end
struct Struct <: Composite end

SType(::Type{T}) where {T} = Struct()
SType(::Type{T}) where {T <: Base.Bool} = Boolean()
SType(::Type{T}) where {T <: Base.Nothing} = Null()
SType(::Type{T}) where {T <: Base.Number} = NumberType()
SType(::Type{T}) where {T <: AbstractString} = StringType()
SType(::Type{T}) where {T <: AbstractVector} = ArrayType()
SType(::Type{T}) where {V, T <: AbstractArray{V, 0}} = ArrayType()
SType(::Type{T}) where {T <: Tuple} = TupleType()
SType(::Type{T}) where {T <: AbstractDict} = DictType()

SType(::Type{Symbol}) = StringType()

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

serialize(::Null, x) = nothing
deserialize(::Null, ::Type{T}, x) where {T} = nothing

serialize(::Boolean, x) = x ? true : false
deserialize(::Boolean, ::Type{T}, x) where {T} = T(x)

serialize(::NumberType, x) = x
deserialize(::NumberType, ::Type{T}, x) where {T} = T(x)

serialize(::StringType, x) = string(x)
deserialize(::StringType, ::Type{T}, x) where {T} = T(x)

serialize(::ArrayType, x) = serialize.(x)
deserialize(::ArrayType, ::Type{T}, x) where {T} = [deserialize(eltype(T), _x) for _x in x]

serialize(::TupleType, x) = serialize.(x)
deserialize(::TupleType, ::Type{T}, x) where {T} = T(deserialize(eltype(T), _x) for _x in x)

serialize(::DictType, x) = OrderedDict(serialize(k) => serialize(v) for (k, v) in x)

function deserialize(::DictType, ::Type{T}, x) where {K, V, T <: AbstractDict{K, V}}
    return T(deserialize(K, k) => deserialize(V, v) for (k, v) in x)
end

# Fallback for Any
deserialize(::S, ::Type{Any}, x::T) where {S <: SType, T} = deserialize(T, x)

function serialize(::Struct, x::T) where {T}
    return OrderedDict(
        string(field) => serialize(getfield(x, field))
            for (field, _) in iterate_fields(T)
    )
end

function deserialize(::Struct, ::Type{T}, dict::AbstractDict) where {T}
    valid_fields = Set(fieldnames(T))

    typed_pairs = (
        let key = Symbol(field)
            key => deserialize(fieldtype(T, key), dict[field])
        end
        for field in keys(dict) if Symbol(field) in valid_fields
    )

    extra_pairs = (
        let key = Symbol(field)
            key => dict[field]
        end
        for field in keys(dict) if !(Symbol(field) in valid_fields)
    )

    args = merge((;), NamedTuple(typed_pairs), NamedTuple(extra_pairs))
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
            return OrderedDict(typetag(T) => string(k), pairs...)
        end
    end
    throw(ArgumentError("Invalid type $(T). Valid options are $(keys(opts))"))
end

function deserialize(::TaggedUnion, ::Type{T}, dict::AbstractDict) where {T}
    tag = string(typetag(T))
    type = dict[tag]
    subtype = getfield(options(T), Symbol(type))

    args = NamedTuple(
        let key = Symbol(field)
                key => deserialize(fieldtype(subtype, key), dict[field])
        end
            for field in keys(dict) if field != tag
    )
    return subtype(; args...)
end

function serialize(::EnumType, x::T) where {T}
    opts = options(T)
    for k in keys(opts)
        if opts[k] == x
            return k
        end
    end
    throw(ArgumentError("Invalid value $(x) for type $(T). \
                        Valid options are $(keys(opts))"))
end

function deserialize(::EnumType, ::Type{T}, x) where {T}
    name = Symbol(x)
    return getfield(options(T), name)
end

end
