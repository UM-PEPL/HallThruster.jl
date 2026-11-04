using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

@testset "Serialization" begin
    @testset "Braginskii" begin
        test_subtype(het.ThermalConductivityModel, het.Braginskii())
    end

    @testset "Mitchner" begin
        test_subtype(het.ThermalConductivityModel, het.Mitchner())
    end

    @testset "Landmark" begin
        test_subtype(het.ThermalConductivityModel, het.LANDMARK_conductivity())
    end
end
