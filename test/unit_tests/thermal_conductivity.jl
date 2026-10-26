using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

@testset "Braginskii charge dependence" begin
    cache = (;
        μ = fill(0.7, 8), ne = fill(2.0e17, 8), Tev = fill(8.0, 8),
        Z_eff = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 39.0, Inf],
    )
    κ = zeros(8)
    het.Braginskii()(κ, (; cache))
    coefficients = κ ./ (cache.μ .* cache.ne .* cache.Tev)
    # Braginskii Table II values, plus fractional and high-charge cases.
    expected = [4.664, 4.192809, 3.957, 3.721, 3.604, 3.532843, 3.286262, 3.25]
    for i in eachindex(expected)
        @test coefficients[i] ≈ expected[i] atol = 5.0e-4 rtol = 0.0
    end
end

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
