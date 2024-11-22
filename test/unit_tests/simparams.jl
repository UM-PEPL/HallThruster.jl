using HallThruster: HallThruster as ht
using Test

include("serialization_test_utils.jl")

function test_simparams_serialization()
    @testset "Serialization" begin
        p = ht.SimParams()
        test_roundtrip(p)

        p_dict = ht.OrderedDict(
            :dt => 1e-8,
            :current_control => ht.OrderedDict(
                :type => "PIDController",
                :target_value => 15.0,
            ),
        )
        test_roundtrip(ht.SimParams, p_dict)
    end
end

function test_simparams()
    test_simparams_serialization()
end
