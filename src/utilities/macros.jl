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
    return quote
        StructTypes.StructType(::Type{$type}) = StructTypes.StringType()

        function StructTypes.construct(::Type{$type}, x::String; kw...)
            sym = Symbol(x)
            if haskey($options, sym)
                return $options[sym]
            else
                throw(
                    ArgumentError("Invalid $(nameof($type)), select one of $(keys($options)).")
                )
            end
        end

        function Base.string(x::$type)
            for k in keys($options)
                if $options[k] == x
                    return string(k)
                end
            end
            throw(
                ArgumentError("Invalid $(nameof($type)), select one of $(keys($options)).")
            )
        end
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

json3.write(Concrete1(1))                                       // {"type": "Concrete1", "x": 1}
json3.read(\"\"\"{"type": "Concrete2", "s": "test" }\"\"\")     // Concrete2("test")
```

"""
macro __register_abstracttype(type, options)
    return quote
        StructTypes.StructType(::Type{$type}) = StructTypes.AbstractType()
        StructTypes.subtypes(::Type{$type}) = $options
        StructTypes.subtypekey(::Type{$type}) = :type

        StructTypes.StructType(::Type{T}) where {T <: $type} = StructTypes.DictType()

        function StructTypes.construct(::Type{T}, d::Dict; kw...) where {T <: $type}
            return T((d[name] for name in fieldnames(T))...)
        end

        function StructTypes.keyvaluepairs(x::T) where {T <: $type}
            p1 = :type => nameof(T)
            return [p1; [name => getfield(x, name) for name in fieldnames(T)]]
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
