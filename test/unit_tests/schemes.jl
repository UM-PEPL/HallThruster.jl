using HallThruster: HallThruster as het

include("serialization_test_utils.jl")

function test_scheme_serialization()
    @testset "Serialization" begin
        scheme1 = het.HyperbolicScheme()
        test_roundtrip(het.HyperbolicScheme, scheme1)

        dict1 = het.Serialization.OrderedDict(
            :flux_function => "rusanov"
        )
        test_roundtrip(het.HyperbolicScheme, dict1)

        dict2 = het.Serialization.OrderedDict(
            :reconstruct => false
        )
        test_roundtrip(het.HyperbolicScheme, dict2)
        scheme = het.deserialize(het.HyperbolicScheme, dict2)
        @test scheme.flux_function == scheme1.flux_function
        @test scheme.limiter == scheme1.limiter
        @test scheme.reconstruct == false
    end
end

function test_schemes()
    @testset "Hyperbolic schemes" begin
        test_scheme_serialization()
    end
end
