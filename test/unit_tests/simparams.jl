using HallThruster: HallThruster as ht
using Test

include("serialization_test_utils.jl")

function test_simparams_serialization()
    p = ht.SimParams()
    test_roundtrip(p)
end

function test_simparams()
end
