using HallThruster: HallThruster as ht
using Test

include("serialization_test_utils.jl")

function test_controller_serialization()
    @testset "Serialization" begin
        none = ht.NoController()
        @show ht.serialize(none)

        test_roundtrip(ht.NoController())
        test_roundtrip(ht.PIDController(target_value = 15.0))

        input = ht.Serialization.OrderedDict(
            :type => "PIDController",
            :target_value => 15.0,
            :integral_constant => 1.0,
        )

        control = ht.deserialize(ht.CurrentController, input)
        @test control isa PIDController
        @test control.target_value == 15.0
        @test control.integral_constant == 1.0
        test_roundtrip(ht.CurrentController, input)
    end
end

function test_current_controller()
end
