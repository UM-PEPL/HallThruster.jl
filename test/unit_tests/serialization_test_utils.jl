using Test
using HallThruster: HallThruster as ht, JSON as js

function struct_eq(x::T, y::T) where {T}
    if x != y
        for f in fieldnames(T)
            xv = getfield(x, f)
            yv = getfield(y, f)
            if !(xv == yv || struct_eq(xv, yv))
                return false
            end
        end
    end
    return true
end

function test_roundtrip(::Type{T}, x::X) where {T, X <: T}
    dict = ht.serialize(x)
    json = js.json(dict)
    dict2 = js.parse(json)
    x2 = ht.deserialize(T, dict2)
    return @test struct_eq(x, x2)
end

function test_roundtrip(::Type{T}, x::AbstractDict) where {T}
    obj1 = ht.deserialize(T, x)
    return test_roundtrip(T, obj1)
end

function test_instances(::Type{T}, instances::NamedTuple) where {T}
    for instance in instances
        test_roundtrip(T, instance)
    end
    return
end

function test_subtype(::Type{T}, subtype::S; show_js = false) where {T, S <: T}
    dict = ht.serialize(subtype)

    @test dict["type"] == string(nameof(S))
    @test length(keys(dict)) == fieldcount(S) + 1
    if show_js
        js.json(dict)
        println()
    end

    return test_roundtrip(T, subtype)
end
