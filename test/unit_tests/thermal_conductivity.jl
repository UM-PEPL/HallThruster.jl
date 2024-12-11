using HallThruster: HallThruster as het

function test_conductivity_serialization()
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
end

function test_thermal_conductivity()
    @testset "Thermal conductivity" begin
        test_conductivity_serialization()
    end
end
