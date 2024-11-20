"""
    @__register_stringtype(type, options)
Registers a type as serializable to JSON as a StringType
This means that that the type should be read as a string, and written as a string,
but may be something else internally.

# Arguments 
---
- type: the type to register 
- options: a NamedTuple mapping symbols to values for this type.

# Example
---
```julia
struct MyStruct
    x::Int
end

const struct_1 = MyStruct(1)
const struct_2 = MyStruct(2)
const struct_3 = MyStruct(3)
const structs = (;struct_1, struct_2, struct_3)

@ __register_stringtype(MyStruct, structs)

s = JSON3.write(struct_1) // "struct1"
s = JSON3.read("struct2") // struct2
```
"""
macro __register_stringtype(type, options)
    @eval function serialize(x::$(type))
        for (k, v) in pairs($(options))
            if v == x
                return string(k)
            end
        end
        throw(ArgumentError("Invalid option $(x) for type $($(type))"))
    end

    @eval function deserialize(::Type{$(type)}, s)
        return $(options)[Symbol(s)]
    end
end

"""
    @__register_abstracttype(type, options)
Registers a type as serializable to an AbstractType
Works the same as the normal StructTypes.AbstractType(), but preserves the "type: " key
when serialized.

# Arguments 
---
- type: the type to register 
- options: a NamedTuple mapping symbols to values for this type.

#Example
---
```julia
abstract type MyAbstract end

struct Concrete1 <: MyAbstract
    x::Int
end

struct Concrete2 <: MyAbstract
    s::String
end

@__register_abstracttype(MyAbstract, (;Concrete1, Concrete2))
```
"""
macro __register_abstracttype(type, options)
    @eval begin
        function serialize(model::T) where {T <: $(type)}
            type_pair = "type" => string(nameof(T))
            kv_pairs = (
                string(field) => serialize(getfield(model, field))
            for field in fieldnames(T)
            )

            return OrderedDict(type_pair, kv_pairs...)
        end

        function deserialize(
                ::Type{T}, dict::AbstractDict,) where {T <: $(type)}
            # Get the specific concrete type
            model = $(options)[Symbol(dict["type"])]

            # Map keys to fieldnames of this type
            args = NamedTuple(
                field => deserialize(field_type, dict[string(field)])
            for (field, field_type) in zip(fieldnames(model), fieldtypes(model))
            )

            # Pass args into keyword constructor of this type
            return model(args...)
        end
    end
end

"""
    @NTuple(ex)
Allows the definition of NTuples from compile-time-known generator expressions

Example:
- x = @NTuple([i for i in 1:3]) // (1,2,3)
"""
macro NTuple(ex)
    if !isa(ex, Expr)
        error("Bad input for @NTuple")
    end
    head = ex.head
    if head === :comprehension
        if length(ex.args) != 1 || !isa(ex.args[1], Expr) || ex.args[1].head != :generator
            error("Expected generator in comprehension, e.g. [f(i) for i = 1:3]")
        end
        ex = ex.args[1]
        if length(ex.args) != 2
            error("Use a one-dimensional comprehension for @NamedTuple")
        end
        rng = eval(ex.args[2].args[2])
        exprs = (:(f($j)) for j in rng)
        return quote
            let
                f($(esc(ex.args[2].args[1]))) = $(esc(ex.args[1]))
                NTuple{$(length(rng))}($tuple($(exprs...)))
            end
        end
    elseif head === :typed_comprehension
        if length(ex.args) != 2 || !isa(ex.args[2], Expr) || ex.args[2].head != :generator
            error(
                "Expected generator in typed comprehension, e.g. Float64[f(i) for i = 1:3]"
            )
        end
        T = esc(ex.args[1])
        ex = ex.args[2]
        if length(ex.args) != 2
            error("Use a one-dimensional comprehension for @NamedTuple")
        end
        rng = eval(ex.args[2].args[2])
        exprs = (:(f($j)) for j in rng)
        return quote
            let
                f($(esc(ex.args[2].args[1]))) = $(esc(ex.args[1]))
                NTuple{$(length(rng)), $T}($tuple($(exprs...)))
            end
        end
    else
        error("Expected comprehension")
    end
end
